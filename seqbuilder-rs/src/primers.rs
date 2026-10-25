// Primer analysis as SeqBuilder's Primer Design view shows it: Tm and dG, dimers, hairpins and binding
// sites.
//
// The formulas and constants here were read out of SeqBuilder Pro 19's own code (Ghidra):
// PCRConditionData holds _deltaH_DNA and _deltaS_DNA, 4x4 tables of doubles indexed
// [second base][first base] over A,C,G,T, whose values are SantaLucia (1998) unified; makeDeltaG()
// fills a third table with ΔS·T/1000 − ΔH, the magnitude of ΔG at the ΔG temperature; and
// ComputePrimingSites::CalcFalsePrimerStats sums those tables and writes Tm and dG.
// _deltaS_DNAini is 0.0108 and the salt and primer concentrations come from the Create Primer Pairs
// dialog (50 mM, 250 nM).
const SALT: f64 = 0.05;      // [Na+] in M (PCRConditionData salt field, mM)
const PRIMER_CT: f64 = 250e-9; // primer concentration in M; the Tm uses Ct/2
const DG_TEMP_C: f64 = 25.0; // ΔG temperature
const DS_INI: f64 = 0.0108;  // _deltaS_DNAini

fn complement(b: u8) -> u8 {
    match b { b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G', _ => b'N' }
}

/// ΣΔH (kcal/mol) and ΣΔS (cal/K·mol) over the nearest neighbours; None for non-ACGT bases.
fn nn_sums(p: &[u8]) -> Option<(f64, f64)> {
    if p.len() < 2 || p.iter().any(|b| !b"ACGT".contains(b)) { return None; }
    let nn = |d: &[u8]| -> (f64, f64) {
        match d {
            b"AA" | b"TT" => (-7.9, -22.2), b"AT" => (-7.2, -20.4), b"TA" => (-7.2, -21.3),
            b"CA" | b"TG" => (-8.5, -22.7), b"GT" | b"AC" => (-8.4, -22.4), b"CT" | b"AG" => (-7.8, -21.0),
            b"GA" | b"TC" => (-8.2, -22.2), b"CG" => (-10.6, -27.2), b"GC" => (-9.8, -24.4), _ => (-8.0, -19.9), // GG, CC
        }
    };
    Some(p.windows(2).map(nn).fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1)))
}

fn upper(p: &[u8]) -> Vec<u8> {
    p.iter().map(|b| b.to_ascii_uppercase()).collect()
}

/// ΔG of a perfectly paired stretch at the ΔG temperature (kcal/mol), without the initiation term.
fn duplex_dg(s: &[u8]) -> f64 {
    nn_sums(s).map_or(0.0, |(dh, ds)| dh - (DG_TEMP_C + 273.15) * ds / 1000.0)
}

/// Tm: ΣΔH·1000 / (0.368·(N−1)·ln[Na+] + ΣΔS + ΔS_ini + R·ln(Ct/2)) − 273.15, as
/// CalcFalsePrimerStats computes it for a perfectly matched primer.
pub fn tm(p: &[u8]) -> f64 {
    let p = upper(p);
    let Some((dh, ds)) = nn_sums(&p) else { return f64::NAN };
    let denom = 0.368 * (p.len() - 1) as f64 * SALT.ln() + ds + DS_INI + 1.987 * (PRIMER_CT / 2.0).ln();
    dh * 1000.0 / denom - 273.15
}

/// The whole primer's ΔG, SeqBuilder's dG column: the nearest-neighbour sum at the ΔG temperature plus
/// the initiation term ΔS_ini × T(°C) = 0.27 kcal/mol at 25 °C.
pub fn dg(p: &[u8]) -> f64 {
    let p = upper(p);
    if nn_sums(&p).is_none() { return f64::NAN; }
    duplex_dg(&p) + DS_INI * DG_TEMP_C
}

/// Two primers paired antiparallel: `offset` is where `b`'s 3' end sits under `a` (index in `a`).
pub struct Dimer {
    pub offset: isize,
    pub pairs: Vec<bool>, // per base of `a`: paired with `b` at this offset
    pub len: usize,       // longest paired stretch
    pub dg: f64,          // its ΔG, kcal/mol
}

/// The most stable dimer of `a` with `b` that a polymerase could extend: the alignment must hold the
/// 3' end of `a` or of `b`. Stability is that of the best contiguous paired stretch, as SeqBuilder
/// reports it ("Most stable dimer 3 bp, -3.3 kc/m" for ALDH2's reverse primer with itself).
pub fn dimer(a: &[u8], b: &[u8]) -> Option<Dimer> {
    let (a, b) = (upper(a), upper(b));
    let (na, nb) = (a.len() as isize, b.len() as isize);
    let mut best: Option<Dimer> = None;
    for offset in -(nb - 1)..na {
        // b runs 3'→5' under a: b's 3' end (j = nb-1) sits under a[offset], then j falls as i rises
        let j_of = |i: isize| offset - i + nb - 1;
        let overlap: Vec<isize> = (0..na).filter(|&i| (0..nb).contains(&j_of(i))).collect();
        let holds_3prime = overlap.contains(&(na - 1)) || overlap.iter().any(|&i| j_of(i) == nb - 1);
        if !holds_3prime { continue; }
        let pairs: Vec<bool> = (0..na).map(|i| {
            let j = j_of(i);
            (0..nb).contains(&j) && complement(a[i as usize]) == b[j as usize]
        }).collect();
        // best contiguous stretch of at least 2 pairs
        let mut i = 0;
        while i < pairs.len() {
            if !pairs[i] { i += 1; continue; }
            let start = i;
            while i < pairs.len() && pairs[i] { i += 1; }
            if i - start < 2 { continue; }
            let g = duplex_dg(&a[start..i]);
            // ties go to the later alignment, the stretch nearer the 3' end (SeqBuilder picks CCGG over CCCG
            // for pBR322's reverse 1483..1506)
            if best.as_ref().is_none_or(|d| g <= d.dg) {
                best = Some(Dimer { offset, pairs: pairs.clone(), len: i - start, dg: g });
            }
        }
    }
    best
}

/// A hairpin: stem of `stem` pairs closing a loop between `start + stem` and `end - stem`.
pub struct Hairpin {
    pub start: usize,
    pub end: usize, // inclusive, the base pairing with `start`
    pub stem: usize,
    pub dg: f64,
}

// ponytail: one loop penalty for every loop size, fitted to SeqBuilder (ALDH2 "3 bp, -0.2", pBR322 reverse
// 1483..1506 "4 bp, -2.3"); the two allow 4.31-4.37.
const HAIRPIN_LOOP: f64 = 4.34;

/// The most stable hairpin (stem of at least 2 pairs, loop of at least 3 bases) with ΔG below zero;
/// SeqBuilder says "No hairpins > 1 bp" otherwise.
pub fn hairpin(p: &[u8]) -> Option<Hairpin> {
    let p = upper(p);
    let mut best: Option<Hairpin> = None;
    for start in 0..p.len() {
        for end in (start + 1..p.len()).rev() {
            let mut stem = 0;
            while start + stem < end - stem && complement(p[start + stem]) == p[end - stem] { stem += 1; }
            if stem < 2 || (end - stem) < (start + stem) + 2 { continue; } // loop of 3 or more
            let g = duplex_dg(&p[start..start + stem]) + HAIRPIN_LOOP;
            if g < 0.0 && best.as_ref().is_none_or(|h| g < h.dg) { best = Some(Hairpin { start, end, stem, dg: g }); }
        }
    }
    best
}

/// A hairpin is marked "(BAD!)" when its ΔG reaches −1.0 (DisplayObjMisprimingHairpin::DrawHairpin
/// tests the stored magnitude against 1.0).
pub fn hairpin_is_bad(dg: f64) -> bool {
    dg <= -1.0
}

// ponytail: DrawDimer also has a softer rule that looks at the 3'-terminal pair; only the plain
// threshold is implemented, so a dimer between -2.07 and -2.53 may be marked differently.
/// A dimer is marked "(BAD!)" when it is more stable than a single GC step, which is the threshold
/// DrawDimer reads out of the ΔG table (ΔG(GC) at the ΔG temperature).
pub fn dimer_is_bad(dg: f64) -> bool {
    dg < duplex_dg(b"GC")
}

/// Where the primer's 3'-terminal 12 bases (SeqBuilder's "Unique at 3' for 12 bp") match the sequence, on
/// the top strand or with `reverse` on the bottom strand: the 0-based position of the primer's 3' end.
pub fn binding_sites(seq: &[u8], primer: &[u8], reverse: bool) -> Vec<usize> {
    let k = 12.min(primer.len());
    let tail = upper(&primer[primer.len() - k..]);
    let target: Vec<u8> = if reverse { tail.iter().rev().map(|&b| complement(b)).collect() } else { tail };
    let seq = upper(seq);
    let three_prime = if reverse { 0 } else { k - 1 };
    seq.windows(k).enumerate().filter(|(_, w)| *w == target.as_slice()).map(|(i, _)| i + three_prime).collect()
}

/// Product Tm: 81.5 + 16.6·log10[Na+] + 41·GC fraction − 675/N, the GC-based formula
/// CalcFalsePrimerStats uses.
pub fn product_tm(product: &[u8]) -> f64 {
    let n = product.len() as f64;
    let gc = product.iter().filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C')).count() as f64 / n;
    81.5 + 16.6 * SALT.log10() + 41.0 * gc - 675.0 / n
}

/// The pair table of the Primer Design view: (Prod Tm, Pri dTm, Ta OPT, Prod dTm), where
/// SeqDataHolder::ComputeProductAnnealingTemperature gives Ta OPT = 0.7·ProdTm + 0.3·min(Tm) − 14.9
/// and Prod dTm = ProdTm − min(Tm).
pub fn pair_values(forward: &[u8], reverse: &[u8], product: &[u8]) -> (f64, f64, f64, f64) {
    let (tf, tr, tp) = (tm(forward), tm(reverse), product_tm(product));
    let low = tf.min(tr);
    (tp, (tf - tr).abs(), 0.3 * low + 0.7 * tp - 14.9, tp - low)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aldh2_reverse_primer_self_dimer() {
        let rev = b"CCGCCTCCCGGGTTCATGCCATTCTCCTGCCTCAGCCT";
        let d = dimer(rev, rev).unwrap();
        // SeqBuilder: "Most stable dimer 3 bp, -3.3 kc/m"
        assert_eq!((d.len, format!("{:.1}", d.dg)), (3, "-3.3".to_string()));
        assert_eq!(format!("{:.1}", tm(rev)), "77.7");
    }

    /// The nearest-neighbour tables read out of SeqBuilder's _deltaH_DNA and _deltaS_DNA.
    #[test]
    fn nn_tables_are_santalucia() {
        assert_eq!(nn_sums(b"CG"), Some((-10.6, -27.2)));
        assert_eq!(nn_sums(b"GC"), Some((-9.8, -24.4)));
        assert_eq!(nn_sums(b"AT"), Some((-7.2, -20.4)));
        assert_eq!(nn_sums(b"TA"), Some((-7.2, -21.3)));
        // the dimer threshold DrawDimer reads from the table
        assert_eq!(format!("{:.3}", duplex_dg(b"GC")), "-2.525");
    }

    /// Values SeqBuilder's Primer Design view showed for pBR322 pairs and ALDH2.
    #[test]
    fn matches_seqbuilder_primer_design() {
        let f704 = b"CGGGGCATGACTATCG";
        let r1483 = b"GCAACCCCGCCAGCCTAGCCGGGT";
        let d = dimer(r1483, r1483).unwrap();
        assert_eq!((d.len, format!("{:.1}", d.dg)), (4, "-6.6".to_string())); // "Most stable dimer 4 bp, -6.6 kc/m (BAD!)"
        assert_eq!(d.offset, 16); // the 3'-end CCGG, as SeqBuilder draws it
        let h = hairpin(r1483).unwrap();
        assert_eq!((h.stem, format!("{:.1}", h.dg)), (4, "-2.3".to_string())); // "Hairpin 4 bp, -2.3 kc/m (BAD!)"
        assert!(dimer_is_bad(d.dg) && hairpin_is_bad(h.dg));
        assert!(hairpin(f704).is_none()); // "No hairpins > 1 bp"
        let h = hairpin(b"CCGCCTCCCGGGTTCATGCCATTCTCCTGCCTCAGCCT").unwrap();
        assert_eq!((h.stem, format!("{:.1}", h.dg)), (3, "-0.2".to_string())); // ALDH2 "Hairpin 3 bp, -0.2 kc/m"
        assert!(!hairpin_is_bad(h.dg)); // -0.2 is above the -1.0 hairpin threshold

        let seq = std::fs::read("testdata/pBR322.gb").unwrap();
        let seq = gb_io::reader::parse_slice(&seq).unwrap().remove(0).seq;
        let fmt = |v: (f64, f64, f64, f64)| [v.0, v.1, v.2, v.3].map(|x| format!("{x:.1}"));
        // Set 2: Forward 704..719 + Reverse 1483..1506, 803 bp: Prod Tm 84.3, Pri dTm 17.8, Ta OPT 61.3, Prod dTm 27.0
        assert_eq!(fmt(pair_values(f704, r1483, &seq[703..1506])), ["84.3", "17.8", "61.3", "27.0"]);
        // Set 1: Forward 635..659 + Reverse 1267..1284, 650 bp: Tm 84.2, Ta OPT 61.3, dTm 26.4
        let v = fmt(pair_values(b"TCGCATAAGGGAGAGCGTCGACCGA", b"GCTTCCATTCAGGTCGAG", &seq[634..1284]));
        assert_eq!((v[0].as_str(), v[2].as_str(), v[3].as_str()), ("84.2", "61.3", "26.4"));
        // only each primer's own site, as SeqBuilder marks them on pBR322
        assert_eq!(binding_sites(&seq, f704, false), vec![718]); // 3' end of 704..719
        assert_eq!(binding_sites(&seq, r1483, true), vec![1482]); // 3' end of 1483..1506 on the bottom strand
    }
}
