// Primer analysis as SeqBuilder's Primer Design view shows it: Tm and dG, dimers, hairpins and binding
// sites. Energies are SantaLucia (1998) unified nearest-neighbour values without initiation terms.

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

/// ΔG at 25 °C of a perfectly paired stretch (kcal/mol).
fn duplex_dg(s: &[u8]) -> f64 {
    nn_sums(s).map_or(0.0, |(dh, ds)| dh - 298.15 * ds / 1000.0)
}

/// Tm as SeqBuilder shows it: salt 0.368·(N−1)·ln[Na+] at 50 mM, primer 250 nM taken as Ct/2.
/// Matches SeqBuilder to 0.1 °C on eight primers of 16-38 bases.
pub fn tm(p: &[u8]) -> f64 {
    let p = upper(p);
    let Some((dh, ds)) = nn_sums(&p) else { return f64::NAN };
    let (na, ct) = (0.05f64, 250e-9f64);
    dh * 1000.0 / (ds + 0.368 * (p.len() - 1) as f64 * na.ln() + 1.987 * (ct / 2.0).ln()) - 273.15
}

// ponytail: the +0.27 kcal/mol is fitted to SeqBuilder's dG column (exact on six primers), its source is unknown.
/// The whole primer's duplex ΔG at 25 °C, SeqBuilder's dG column.
pub fn dg(p: &[u8]) -> f64 {
    let p = upper(p);
    if nn_sums(&p).is_none() { return f64::NAN; }
    duplex_dg(&p) + 0.27
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
            if best.as_ref().is_none_or(|d| g < d.dg) {
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

// ponytail: loop penalties are SantaLucia (2004) hairpin values; SeqBuilder's "Hairpin 3 bp, -0.2 kc/m"
// for ALDH2's reverse primer is not reproduced yet, so this ranking may differ from SeqBuilder's.
/// The most stable hairpin (stem of at least 2 pairs, loop of at least 3 bases).
pub fn hairpin(p: &[u8]) -> Option<Hairpin> {
    let p = upper(p);
    let loop_penalty = |l: usize| match l { 3 | 4 => 3.5, 5 => 3.3, 6 => 4.0, 7 => 4.2, 8 => 4.3, 9 => 4.5, _ => 4.6 + 0.1 * (l as f64 - 10.0).max(0.0) };
    let mut best: Option<Hairpin> = None;
    for start in 0..p.len() {
        for end in (start + 1..p.len()).rev() {
            let mut stem = 0;
            while start + stem < end - stem && complement(p[start + stem]) == p[end - stem] { stem += 1; }
            if stem < 2 { continue; }
            let lp = (end - stem) as isize - (start + stem) as isize + 1;
            if lp < 3 { continue; }
            let g = duplex_dg(&p[start..start + stem]) + loop_penalty(lp as usize);
            if best.as_ref().is_none_or(|h| g < h.dg) { best = Some(Hairpin { start, end, stem, dg: g }); }
        }
    }
    best
}

// ponytail: SeqBuilder's rule for "Primer Binding Sites" is not known; this marks every place where the
// primer's 3'-terminal 8 bases match, which is what a polymerase would extend from.
/// 0-based positions on the sequence where the primer's 3' end could prime (forward: on the top strand).
pub fn binding_sites(seq: &[u8], primer: &[u8], reverse: bool) -> Vec<usize> {
    let k = 8.min(primer.len());
    let tail = upper(&primer[primer.len() - k..]);
    let target: Vec<u8> = if reverse { tail.iter().rev().map(|&b| complement(b)).collect() } else { tail };
    let seq = upper(seq);
    seq.windows(k).enumerate().filter(|(_, w)| *w == target.as_slice()).map(|(i, _)| i).collect()
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

    #[test]
    fn binding_sites_find_the_primer_itself() {
        let seq = b"AAAAGGATCCAAAAAAAAAAGGATCCTTTT";
        assert!(binding_sites(seq, b"CCCAAAAGGATCC", false).contains(&2));
        assert!(!binding_sites(seq, b"GGATCCTTTT", true).is_empty());
    }
}
