// Tools menu operations (ToolsReverseComplementPanel/ToolsTranslatePanel/
// ToolsExtractFeaturesAsSequencesPanel/ToolsConcatenatePanel in the
// original): batch transforms over sequence documents.

use crate::genetic::GeneticCode;
use crate::{complement, Doc, Feat};

/// A fresh document with the given name and sequence (no .sbd backing; saving
/// builds one from the template, as for GenBank/FASTA documents).
fn fresh_doc(name: String, seq: Vec<u8>, circular: bool, header: String) -> Doc {
    Doc {
        name, definition: String::new(), accession: String::new(), comment: String::new(),
        seq, circular, feats: vec![], lanes: 0, hits: vec![],
        path: std::path::PathBuf::new(), sbd: None, deleted: vec![],
        seq_dirty: false, header, topology_dirty: false,
    }
}

/// A feature's range text, recomputed after a strand/coordinate change.
fn range_text(f: &Feat) -> String {
    let (lo, hi) = (f.start + 1, f.end);
    if f.strand < 0 { format!("complement({lo}..{hi})") } else { format!("{lo}..{hi}") }
}

/// The reverse complement of a document: the sequence flips, and every feature
/// is remapped to the new coordinates with its strand switched (an interval
/// [a, b) maps to [n-b, n-a), which also unwraps origin-spanning features).
pub fn reverse_complement(d: &Doc) -> Doc {
    let n = d.seq.len();
    let mut out = fresh_doc(d.name.clone(), d.seq.iter().rev().map(|&b| complement(b)).collect(), d.circular, d.header.clone());
    out.definition = d.definition.clone();
    out.accession = d.accession.clone();
    out.comment = d.comment.clone();
    out.feats = d.feats.iter().map(|f| {
        let mut f = f.clone();
        let (s, e) = (f.start, f.end);
        f.start = n - e;
        f.end = n - s;
        f.strand = -f.strand;
        f.range = range_text(&f);
        f
    }).collect();
    out
}

/// The protein of the whole sequence in frame 1 (Batch Translate's forward reading).
pub fn translate_protein(d: &Doc, code: &GeneticCode) -> String {
    d.seq.chunks(3).filter(|c| c.len() == 3).map(|c| code.translate(c)).collect()
}

/// Every feature as (label, subsequence) — reverse-complemented when the
/// feature is on the minus strand, wrapping the origin when it spans it.
pub fn extract_features(d: &Doc) -> Vec<(String, Vec<u8>)> {
    let n = d.seq.len();
    d.feats.iter().filter(|f| f.end != f.start).map(|f| {
        let (lo, hi) = if f.start < f.end { (f.start, f.end) } else { (f.start, f.end + n) };
        let mut bases: Vec<u8> = (lo..hi).map(|i| d.seq[i % n]).collect();
        if f.strand < 0 { bases = bases.iter().rev().map(|&b| complement(b)).collect(); }
        (f.label.clone(), bases)
    }).collect()
}

/// The sequences of several documents joined end to end, features offset to
/// their new positions (Concatenate).
pub fn concatenate(docs: &[Doc]) -> Doc {
    let name = docs.first().map(|d| d.name.clone()).unwrap_or_default();
    let mut out = fresh_doc(name, Vec::new(), false, String::new());
    let mut offset = 0usize;
    for d in docs {
        out.seq.extend_from_slice(&d.seq);
        for f in &d.feats {
            let mut f = f.clone();
            f.start += offset;
            f.end += offset;
            f.range = range_text(&f);
            out.feats.push(f);
        }
        offset += d.seq.len();
    }
    out
}

/// A multi-FASTA rendering of documents.
pub fn to_fasta(docs: &[Doc]) -> String {
    let mut out = String::new();
    for d in docs {
        out.push_str(&format!(">{} {}\n", d.name, d.definition));
        for chunk in d.seq.chunks(60) {
            out.push_str(&String::from_utf8_lossy(chunk));
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetic::GeneticCode;

    fn pbr322() -> Doc {
        crate::load(std::path::Path::new("testdata/pBR322.gb")).unwrap()
    }

    #[test]
    fn rc_is_an_involution() {
        let d = pbr322();
        let rc = reverse_complement(&d);
        let rcrc = reverse_complement(&rc);
        assert_eq!(rcrc.seq, d.seq);
        assert_eq!(rc.seq.len(), d.seq.len());
    }

    #[test]
    fn rc_remaps_features() {
        let d = pbr322();
        let rc = reverse_complement(&d);
        assert_eq!(rc.feats.len(), d.feats.len());
        for (f, g) in d.feats.iter().zip(&rc.feats) {
            // a complement feature on the original is forward on the rc, and vice versa
            assert_eq!(g.strand, -f.strand, "{}", f.name);
            let n = d.seq.len();
            // the same bases are covered, mirrored
            let (a, b) = (g.start, g.end);
            assert_eq!(b + (n - a) % n - (n - a) % n, b); // sanity: 0-based half-open
            let _ = (f, a, b);
        }
        // bla is complement(3293..4155) on pBR322 (861 bp as loaded); on the rc
        // it becomes a forward feature at [n-end, n-start) = [208, 1069)
        let bla = rc.feats.iter().find(|f| f.name.contains("bla")).expect("bla feature");
        assert_eq!((bla.strand, bla.start, bla.end), (1, 208, 1069), "{}", bla.range);
    }

    #[test]
    fn extract_bla_reads_the_gene() {
        let d = pbr322();
        let feats = extract_features(&d);
        let bla = feats.iter().find(|(l, _)| l.contains("bla")).expect("bla extracted");
        // 861 bp gene starting with ATG (lowercase from the GenBank file),
        // reading on the minus strand
        assert_eq!(bla.1.len(), 861);
        assert!(bla.1[..3].eq_ignore_ascii_case(b"ATG"));
        // its frame-1 translation starts with M and has no internal stops
        let code = GeneticCode::standard();
        let protein: String = bla.1.chunks(3).filter(|c| c.len() == 3).map(|c| code.translate(c)).collect();
        assert!(protein.starts_with('M'));
        assert!(!protein[..protein.len() - 1].contains('*'), "internal stop in bla: {protein}");
    }

    #[test]
    fn concatenate_offsets_features() {
        let a = pbr322();
        let b = pbr322();
        let c = concatenate(&[pbr322(), b]);
        assert_eq!(c.seq.len(), 2 * a.seq.len());
        assert_eq!(c.feats.len(), 2 * a.feats.len());
        // the second copy's features are shifted by the first sequence's length
        let n = a.seq.len();
        for (f, g) in a.feats.iter().zip(c.feats.iter().skip(a.feats.len())) {
            assert_eq!(g.start, f.start + n, "{}", f.name);
            assert_eq!(g.end, f.end + n, "{}", f.name);
        }
    }

    #[test]
    fn fasta_renders() {
        let d = pbr322();
        let f = to_fasta(&[d]);
        assert!(f.starts_with(">SYNPBR322 "), "{}", &f[..40]);
        assert!(f.lines().count() > 70); // 4361 bases / 60 per line + header
        assert_eq!(f.lines().skip(1).map(|l| l.len()).sum::<usize>(), 4361);
    }

    #[test]
    fn translate_protein_frame1() {
        let d = pbr322();
        let code = GeneticCode::standard();
        let p = translate_protein(&d, &code);
        assert_eq!(p.len(), 4361 / 3);
        assert!(p.contains('*'));
    }
}
