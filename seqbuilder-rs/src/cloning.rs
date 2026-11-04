// Cloning suite, matching SeqBuilder Pro's Cloning menu (CloningCommands/
// CloneFragmentDialog in the original): restriction inserts with real
// overhang geometry, Gibson/In-Fusion-style homology assembly, and fragment
// extraction.
//
// Restriction fragment geometry: a fragment cut out between two enzymes'
// top-strand cuts at [a, b) carries the top strand seq[a..b); the bottom
// strand starts left_d bases after a and ends right_d bases after b
// (d = bot_cut - top_cut of the respective enzyme). With d > 0 the left end
// carries a 5' top overhang (top[a..a+d)) and the right end a 5' bottom
// overhang (the complement of the source's top[b..b+d), read 5'→3' on the
// bottom strand). Blunt ends have d = 0; 3' overhangs (d < 0, KpnI-like)
// are not supported yet.

use crate::complement;

/// A restriction fragment held for pasting (Cloning ▸ Copy Restriction Insert).
#[derive(Clone)]
pub struct Insert {
    pub name: String,
    pub top: Vec<u8>,         // top strand between the cuts
    pub left_d: i32,          // bottom strand starts this far after the top start
    pub right_d: i32,         // bottom strand ends this far after the top end
    pub left_enz: String,
    pub right_enz: String,
    right_over: Vec<u8>,      // right-end bottom overhang, 5'→3'
}

fn revcomp(s: &[u8]) -> Vec<u8> { s.iter().rev().map(|&b| complement(b)).collect() }

/// `k` bases of `seq` starting at `start`, wrapping a circular sequence.
fn take(seq: &[u8], start: usize, k: usize) -> Vec<u8> {
    let n = seq.len();
    (0..k).map(|i| seq[(start + i) % n]).collect()
}

/// Cut an insert between top-strand cut positions [a, b) using the enzymes
/// that made the cuts; `d` values are each enzyme's bot_cut - top_cut.
pub fn copy_insert(seq: &[u8], a: usize, b: usize, left_d: i32, right_d: i32, left_enz: &str, right_enz: &str, name: &str) -> Result<Insert, String> {
    if left_d < 0 || right_d < 0 { return Err("3' overhang ends (KpnI-like cuts) are not supported yet".into()); }
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    if b > seq.len() { return Err("selection past the end of the sequence".into()); }
    let top = seq[a..b].to_vec();
    let right_over = if right_d > 0 { revcomp(&take(seq, b, right_d as usize)) } else { vec![] };
    Ok(Insert { name: name.to_string(), top, left_d, right_d, left_enz: left_enz.to_string(), right_enz: right_enz.to_string(), right_over })
}

impl Insert {
    /// The 5' overhang sequences the fragment carries: (left top-strand,
    /// right bottom-strand), both 5'→3'.
    pub fn overhangs(&self) -> (Vec<u8>, Vec<u8>) {
        let left = if self.left_d > 0 { self.top[..self.left_d as usize].to_vec() } else { vec![] };
        (left, self.right_over.clone())
    }

    pub fn len(&self) -> usize { self.top.len() }
}

/// Whether an insert ligates into a target whose region [x, y) was cut out
/// (left enzyme d = dl, right enzyme d = dr). The insert's left top overhang
/// must anneal to the target's bottom overhang past y, and its right bottom
/// overhang to the target's top overhang past x.
pub fn compatible(insert: &Insert, target: &[u8], x: usize, y: usize, dl: i32, dr: i32) -> Result<(), String> {
    if dl < 0 || dr < 0 { return Err("3' overhang target ends (KpnI-like cuts) are not supported yet".into()); }
    if insert.right_d != dl || insert.left_d != dr {
        return Err(format!("overhang lengths differ: insert {}+{}, target {}+{}",
            insert.left_d, insert.right_d, dl, dr));
    }
    if insert.left_d > 0 {
        let t = take(target, y, insert.left_d as usize);
        if insert.top[..insert.left_d as usize] != t[..] {
            return Err("the left overhang bases do not match (AATT vs target past the right cut)".into());
        }
    }
    if insert.right_d > 0 {
        let t = take(target, x, insert.right_d as usize);
        if insert.right_over != revcomp(&t) {
            return Err("the right overhang bases do not match".into());
        }
    }
    Ok(())
}

/// The ligated product: the target with [x, y) replaced by the insert.
pub fn paste_insert(target: &[u8], insert: &Insert, x: usize, y: usize) -> Vec<u8> {
    let mut out = target[..x].to_vec();
    out.extend_from_slice(&insert.top);
    out.extend_from_slice(&target[y..]);
    out
}

/// Gibson / In-Fusion assembly: join fragments sharing terminal homology of
/// at least `min_overlap` bases (the longest overlap up to 60 bases is used);
/// when the last fragment also overlaps the first, the circle is closed.
pub fn gibson_assemble(fragments: &[Vec<u8>], min_overlap: usize) -> Result<Vec<u8>, String> {
    if fragments.is_empty() { return Err("no fragments".into()); }
    if fragments.len() == 1 { return Ok(fragments[0].clone()); }
    let mut out = fragments[0].clone();
    for w in fragments.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let ov = longest_suffix_prefix(a, b).ok_or("two neighbouring fragments do not overlap")?;
        if ov < min_overlap { return Err(format!("overlap of {ov} bp is below the minimum of {min_overlap}")); }
        out.extend_from_slice(&b[ov..]);
    }
    let (last, first) = (fragments.last().unwrap(), &fragments[0]);
    if let Some(close) = longest_suffix_prefix(last, first).filter(|k| *k >= min_overlap) {
        let n = out.len();
        out.truncate(n - close);
    }
    Ok(out)
}

/// Longest overlap where `a` ends with what `b` starts with (up to 60 bases,
/// shorter than half of either fragment).
pub fn overlap(a: &[u8], b: &[u8]) -> Option<usize> {
    let max = 60.min(a.len() / 2).min(b.len() / 2);
    (1..=max).rev().find(|&k| a.ends_with(&b[..k]))
}

fn longest_suffix_prefix(a: &[u8], b: &[u8]) -> Option<usize> { overlap(a, b) }

/// Clone Selected Fragment: the subsequence of a span, wrapping the origin.
pub fn clone_fragment(seq: &[u8], a: usize, b: usize) -> Vec<u8> {
    if a <= b { seq[a..b].to_vec() } else { seq[a..].iter().chain(seq[..b].iter()).copied().collect() }
}

#[cfg(test)]
mod tests {
    use super::*;

    // EcoRI G^AATTC: d = +4, 5' AATT overhang on both ends.
    const VECTOR: &[u8] = b"AAAAGAATTCCGGGGGGGGGGGAATTCGGTT";

    #[test]
    fn ecori_geometry() {
        let ins = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "fragment").unwrap();
        assert_eq!(ins.top, VECTOR[5..22]);
        let (l, r) = ins.overhangs();
        assert_eq!(l, b"AATT"); // top bases right after the left cut
        assert_eq!(r, b"AATT"); // complement of the source's AATT past the right cut
    }

    #[test]
    fn paste_into_same_enzyme_cuts() {
        let ins = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "insert").unwrap();
        // a target opened with EcoRI on both sides: TTGAATTC...GAATTCCT
        let target = b"TTTTGAATTCAAAAAAAAGAATTCCTTTT";
        // top cuts at 5 and 19
        compatible(&ins, target, 5, 19, 4, 4).unwrap();
        let out = paste_insert(target, &ins, 5, 19);
        assert_eq!(&out[..5], &target[..5]);
        assert_eq!(&out[5..5 + ins.top.len()], &ins.top);
        assert_eq!(&out[5 + ins.top.len()..], &target[19..]);
        // the left junction reconstitutes the EcoRI site: ...G + AATTC...
        assert_eq!(&out[4..10], b"GAATTC");
    }

    #[test]
    fn blunt_and_incompatible() {
        let blunt = Insert { name: "b".into(), top: b"ACGTACGTAC".to_vec(), left_d: 0, right_d: 0, left_enz: "SmaI".into(), right_enz: "SmaI".into(), right_over: vec![] };
        compatible(&blunt, b"AAAAAAAAAA", 3, 7, 0, 0).unwrap();
        // an AATT insert does not ligate into a blunt target
        let sticky = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "s").unwrap();
        assert!(compatible(&sticky, b"AAAAAAAAAAA", 3, 7, 0, 0).is_err());
        // and not into a target whose overhang bases differ
        let t = b"TTTTGACTAGAAAAAAAAGACTAGCTTTT";
        assert!(compatible(&sticky, t, 5, 17, 4, 4).is_err());
    }

    #[test]
    fn circular_target_wraps() {
        // cut at the origin: target[y..] then [x..y) spans the wrap
        let ins = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "i").unwrap();
        let target = b"GAATTCAAAAAAAGAATTCGG";
        // left cut EcoRI at 1, right cut EcoRI at 14: target[..1] ++ ins ++ target[14..]
        compatible(&ins, target, 1, 14, 4, 4).unwrap();
        let out = paste_insert(target, &ins, 1, 14);
        assert_eq!(out.len(), 1 + ins.top.len() + (target.len() - 14));
    }

    #[test]
    fn gibson_two_fragments() {
        let a = b"ACGTACGTACGGGGGGGGGGTTTTTTTTTT".to_vec();
        let b = b"TTTTTTTTTTCCCCCCCCCC".to_vec();
        let out = gibson_assemble(&[a, b], 10).unwrap();
        assert_eq!(out, b"ACGTACGTACGGGGGGGGGGTTTTTTTTTTCCCCCCCCCC");
    }

    #[test]
    fn gibson_three_fragments_circular() {
        let a = b"AAAAAAAAAAGGGGGG".to_vec();
        let b = b"GGGGGGCCCCCC".to_vec();
        let c = b"CCCCCCTTTTTTAAAA".to_vec();
        let out = gibson_assemble(&[a, b, c], 4).unwrap();
        // the last fragment's AAAA tail overlaps the first fragment's head
        assert_eq!(out, b"AAAAAAAAAAGGGGGGCCCCCCTTTTTT");
    }

    #[test]
    fn gibson_rejects_short_overlap() {
        let a = b"ACGTACGTACGTACGTACGTA".to_vec();
        let b = b"GTACGTA".to_vec();
        assert!(gibson_assemble(&[a, b], 10).is_err());
    }

    #[test]
    fn clone_fragment_wraps_origin() {
        let seq = b"AAAACCCCGGGGTTTT";
        assert_eq!(clone_fragment(seq, 12, 4), b"TTTTAAAA");
        assert_eq!(clone_fragment(seq, 4, 8), b"CCCC");
    }

    #[test]
    fn three_prime_overhang_rejected() {
        assert!(copy_insert(VECTOR, 5, 22, -4, 4, "KpnI", "EcoRI", "x").is_err());
    }
}
