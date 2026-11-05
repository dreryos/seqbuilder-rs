// ORF detection matching SeqBuilder Pro's Sequence ▸ ORF Detection Options
// (orfOptions in SeqBuilder.xrc): a minimum ORF size in amino acids and a
// three-way start-codon filter — "Show only ORFs with start codon" /
// "without start codon" / "both". ORFs carry a reading frame (1–6, like the
// frames in the Sequence view) and a strand, are named "ORF n" at display
// time, and can cross the origin of a circular sequence (start > end then,
// the same convention as features).
//
// Criteria assumed from the dialog labels and the binary's ORF classes
// (ORFSearcher/ORFCriteriaDialog/DisplayObjORF, source paths visible in the
// binary): standard genetic code (stops TAA/TAG/TGA, start ATG) until
// selectable genetic codes exist; ORFs run to the sequence edges when no
// stop codon closes them; only whole codons count. The minus-strand span
// uses plus-strand coordinates with the start codon at the high end, like
// complement features.

use crate::complement;
use crate::genetic::GeneticCode;

#[derive(Clone, Copy, PartialEq, Default)]
pub enum StartFilter {
    With,    // Show only ORFs with start codon
    Without, // Show only ORFs without start codon
    #[default]
    Both, // Show both ORFs with and without start codon
}

#[derive(Clone)]
pub struct OrfOpts {
    pub min_aa: usize, // Minimum ORF size, in amino acids
    pub starts: StartFilter,
}

impl Default for OrfOpts {
    fn default() -> Self { Self { min_aa: 30, starts: StartFilter::Both } }
}

#[derive(Clone, Debug)]
pub struct Orf {
    /// Plus-strand span, half-open; start > end wraps the origin of a circular sequence.
    pub start: usize,
    pub end: usize,
    pub strand: i8, // +1 forward, -1 reverse
    pub frame: u8,  // 1..=3 forward, 4..=6 reverse, as displayed
    pub start_codon: bool,
    pub len_aa: usize,
}

fn is_stop(c: &[u8], code: &GeneticCode) -> bool { code.is_stop(c) }
fn is_start(c: &[u8], code: &GeneticCode) -> bool { code.is_start(c) }

/// True when `i` lies inside the (possibly origin-wrapping) span of `o`.
#[allow(dead_code)] // exercised by the orf unit tests; UI hit-testing uses screen rects
pub fn covers(o: &Orf, i: usize) -> bool {
    if o.start <= o.end { i >= o.start && i < o.end } else { i >= o.start || i < o.end }
}

/// Span length in bases, wrap-aware.
pub fn span_len(o: &Orf, n: usize) -> usize {
    if o.start <= o.end { o.end - o.start } else { n - o.start + o.end }
}

/// All ORFs of `seq` under `opts`, sorted for display: strand (+ then -), frame, start.
pub fn find_orfs(seq: &[u8], circular: bool, opts: &OrfOpts, code: &GeneticCode) -> Vec<Orf> {
    let n = seq.len();
    if n < 3 { return vec![] }
    let mut out: Vec<Orf> = Vec::new();
    // Forward strand: raw spans are already plus-strand spans.
    for (s, e, frame, start_codon, aa) in scan(seq, circular, opts, code) {
        let (start, end) = circle_span(n, s, e);
        out.push(Orf { start, end, strand: 1, frame, start_codon, len_aa: aa });
    }
    // Reverse strand: scan the reverse complement, then map spans back to plus coordinates.
    let rc: Vec<u8> = seq.iter().rev().map(|&b| complement(b)).collect();
    for (s, e, f, start_codon, aa) in scan(&rc, circular, opts, code) {
        let (start, end) = reverse_span(n, s, e);
        out.push(Orf { start, end, strand: -1, frame: f + 3, start_codon, len_aa: aa });
    }
    out.sort_by_key(|o| (o.strand, o.frame, o.start));
    out
}

/// Convert a raw span (start < n, end <= start + n) to circle representation:
/// normal [start, end) with end <= n, or start > end when it wraps the origin.
fn circle_span(n: usize, s: usize, e: usize) -> (usize, usize) {
    let len = e - s;
    if len >= n { (0, n) } else if e <= n { (s, e) } else { (s, e - n) }
}

/// Map a reverse-complement span [s, e) (s < n, e <= s + n) to plus-strand
/// coordinates. rc index x is plus index n - 1 - (x mod n), so the span maps
/// to [n - e, n - s) on the circle.
fn reverse_span(n: usize, s: usize, e: usize) -> (usize, usize) {
    let len = e - s;
    if len >= n { return (0, n); }
    let p = (2 * n - e) % n; // (n - e) mod n
    let end = p + len;
    if end <= n { (p, end) } else { (p, end - n) }
}

/// Scan one strand (the caller passes the reverse complement for the minus
/// strand) collecting raw spans (start < n, end <= start + n): frame 1..=3
/// within the strand, start-codon flag, amino-acid length.
fn scan(seq: &[u8], circular: bool, opts: &OrfOpts, code: &GeneticCode) -> Vec<(usize, usize, u8, bool, usize)> {
    let n = seq.len();
    // uppercase so lowercase GenBank sequence blocks parse the same (codons are case-insensitive)
    let up: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    // Circular: the standard concat-and-cap convention — codons tile the doubled
    // sequence from offset f, ORFs start in the first copy and are capped at one
    // circle; the codon phase after the origin simply continues (this is what
    // reading frames mean on a circle whose length is not a multiple of 3).
    let (ext, walk_limit) = if circular {
        let mut e: Vec<u8> = up.iter().chain(up.iter()).copied().collect();
        e.extend_from_slice(&up[..3.min(n)]); // so stop codons at the walk edge are readable
        (e, 2 * n)
    } else {
        (up, n)
    };
    let m = ext.len();
    let mut out = Vec::new();
    for f in 0..3usize {
        // stop-codon starts in this frame
        let mut stops: Vec<usize> = Vec::new();
        let mut i = f;
        while i + 3 <= walk_limit {
            if is_stop(&ext[i..i + 3], code) { stops.push(i); }
            i += 3;
        }
        let last_codon_end = f + (walk_limit - f) / 3 * 3;
        // regions between stops, stop-exclusive: (start, end, ends_at_stop)
        let mut regions: Vec<(usize, usize, bool)> = Vec::new();
        if stops.is_empty() {
            if last_codon_end > f { regions.push((f, last_codon_end, false)); }
        } else {
            if stops[0] > f { regions.push((f, stops[0], true)); }
            for w in stops.windows(2) { regions.push((w[0] + 3, w[1], true)); }
            let last = *stops.last().unwrap();
            if last_codon_end > last + 3 { regions.push((last + 3, last_codon_end, false)); }
        }
        if circular {
            // keep only regions starting in the first copy (drops the duplicates
            // the doubled scan produces), and cap each to one circle
            regions.retain(|&(rs, _, _)| rs < n);
            for r in &mut regions {
                r.1 = r.1.min(r.0 + n);
                r.2 = r.2 && r.1 + 3 <= m && is_stop(&ext[r.1..r.1 + 3], code);
            }
        }
        for (rs, re, closed) in regions {
            // ORFs with a start codon: every in-frame ATG in the region that
            // itself starts in the first copy
            if !matches!(opts.starts, StartFilter::Without) {
                let mut a = rs;
                while a + 3 <= re && a < n {
                    if is_start(&ext[a..a + 3], code) {
                        let e = if closed { re + 3 } else { re };
                        let aa = if closed { (e - a) / 3 - 1 } else { (e - a).div_ceil(3) };
                        if aa >= opts.min_aa { out.push((a, e, f as u8 + 1, true, aa)); }
                    }
                    a += 3;
                }
            }
            // ORFs without a start codon: the inter-stop region itself
            if !matches!(opts.starts, StartFilter::With) {
                let begins_start = re - rs >= 3 && is_start(&ext[rs..rs + 3], code);
                if !begins_start {
                    let aa = (re - rs) / 3;
                    if aa >= opts.min_aa { out.push((rs, re, f as u8 + 1, false, aa)); }
                }
            }
        }
    }
    // cap spans to one full circle (spans already start below n)
    out.iter()
        .map(|&(s, e, fr, sc, aa)| {
            let e = if e - s > n { s + n } else { e };
            (s, e, fr, sc, aa)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn orfs(seq: &str, starts: StartFilter, min_aa: usize) -> Vec<Orf> {
        find_orfs(seq.as_bytes(), false, &OrfOpts { min_aa, starts }, &GeneticCode::standard())
    }

    #[test]
    fn forward_orf_with_start() {
        // ATG TAA: one 1-aa ORF from the start codon through the stop codon
        let v = orfs("ATGTAA", StartFilter::With, 1);
        assert_eq!(v.len(), 1);
        assert_eq!((v[0].start, v[0].end), (0, 6));
        assert_eq!((v[0].strand, v[0].frame, v[0].len_aa), (1, 1, 1));
    }

    #[test]
    fn without_start_regions() {
        // TAA AAA CCC TAA: the frame-1 region between the stops, no start codon
        let v = orfs("TAAAAACCCTAA", StartFilter::Without, 1);
        let f1: Vec<&Orf> = v.iter().filter(|o| o.frame == 1).collect();
        assert_eq!(f1.len(), 1);
        assert_eq!((f1[0].start, f1[0].end), (3, 9));
        assert_eq!(f1[0].len_aa, 2);
    }

    #[test]
    fn region_beginning_with_start_is_not_a_without_orf() {
        let v = orfs("ATGTAA", StartFilter::Without, 1);
        assert!(!v.iter().any(|o| o.frame == 1), "the frame-1 region begins with ATG");
    }

    #[test]
    fn min_size_filters() {
        assert_eq!(orfs("ATGTAA", StartFilter::With, 1).len(), 1);
        assert!(orfs("ATGTAA", StartFilter::With, 2).is_empty());
    }

    #[test]
    fn open_orf_runs_to_edge() {
        // no stop codon: the ORF runs to the sequence end, whole codons only
        let v = orfs("ATGAAAAAA", StartFilter::With, 1);
        assert_eq!(v.len(), 1);
        assert_eq!((v[0].start, v[0].end), (0, 9));
        assert_eq!(v[0].len_aa, 3);
    }

    #[test]
    fn reverse_strand_maps_to_plus_coordinates() {
        // plus = TTAAAACAT, whose reverse complement is ATGTTTTAA:
        // rc frame 1 reads ATG TTT TAA, a 2-aa ORF covering the whole circle
        let v = find_orfs(b"TTAAAACAT", true, &OrfOpts { min_aa: 1, starts: StartFilter::With }, &GeneticCode::standard());
        assert!(v.iter().any(|o| o.strand < 0 && o.len_aa == 2), "reverse ORF over the whole circle");
        // rc [3, 12) of a 12-mer maps to plus [0, 9)
        let v = find_orfs(b"TTAAAACATTTT", true, &OrfOpts { min_aa: 1, starts: StartFilter::With }, &GeneticCode::standard());
        let rev: Vec<&Orf> = v.iter().filter(|o| o.strand < 0).collect();
        assert!(rev.iter().any(|o| o.start == 0 && o.end == 9), "rc [3,12) -> plus [0,9), got {:?}", rev);
    }

    #[test]
    fn circular_wrap() {
        // frame 1 reads TAA ... ATG with no stop between the ATG and the origin:
        // the ATG's ORF is closed by the stop after the origin, so it wraps (start > end).
        // TAA GGG CCC TTT TTT CCC CCC ATG | TAA...
        let seq = b"TAAGGGCCCTTTTTTCCCCCATG";
        let v = find_orfs(seq, true, &OrfOpts { min_aa: 1, starts: StartFilter::With }, &GeneticCode::standard());
        assert!(v.iter().any(|o| o.start > o.end), "ORF crossing the origin wraps, got {v:?}");
        // and the wrapped ORF covers both its parts
        let w = v.iter().find(|o| o.start > o.end).unwrap();
        assert!(covers(w, w.start) && covers(w, w.end - 1));
    }

    #[test]
    fn circular_no_duplicate_regions() {
        // the same circle scanned twice (doubled sequence) must not duplicate ORFs
        let seq = b"ATGAAATAACCCGGGTTT";
        let v = find_orfs(seq, true, &OrfOpts { min_aa: 1, starts: StartFilter::Both }, &GeneticCode::standard());
        let mut keys: Vec<(i8, u8, usize, usize, bool)> =
            v.iter().map(|o| (o.strand, o.frame, o.start, o.end, o.start_codon)).collect();
        keys.sort();
        keys.dedup();
        assert_eq!(keys.len(), v.len(), "duplicates: {v:?}");
    }

    #[test]
    fn frames_are_consistent() {
        // every ORF's frame matches its start offset from the strand's origin
        let seq = b"ATGAAATAACCCGGGTTTAAAGGGCCCATG";
        let n = seq.len();
        for o in find_orfs(seq, false, &OrfOpts { min_aa: 1, starts: StartFilter::Both }, &GeneticCode::standard()) {
            if o.strand > 0 {
                assert_eq!(((o.start % 3) + 1) as u8, o.frame, "{o:?}");
            } else if o.end - o.start < n {
                // for a reverse ORF the rc start codon sits at plus n - end
                let rc_offset = (n - o.end) % 3;
                assert_eq!((rc_offset + 1) as u8, o.frame - 3, "{o:?}");
            }
        }
    }

    #[test]
    fn pbr322_finds_bla() {
        // pBR322's ampR (bla) CDS is complement(3293..4155), 286 aa with a start
        // codon: a with-start ORF of >= 286 aa must exist on the reverse strand.
        let gb = std::fs::read_to_string("testdata/pBR322.gb").unwrap();
        let seq: Vec<u8> = gb.split('\n').skip_while(|l| !l.starts_with("ORIGIN")).skip(1)
            .flat_map(|l| l.chars().filter(|c| c.is_ascii_alphabetic()))
            .map(|c| c as u8).collect();
        let v = find_orfs(&seq, true, &OrfOpts { min_aa: 100, starts: StartFilter::With }, &GeneticCode::from_ncbi(11).unwrap());
        assert!(v.iter().any(|o| o.strand < 0 && o.len_aa >= 286), "bla ORF (286 aa, reverse) not found in {:?}", v.iter().map(|o| (o.strand, o.len_aa)).collect::<Vec<_>>());
    }
}
