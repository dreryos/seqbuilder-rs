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
// bottom strand). Blunt ends have d = 0; with d < 0 the end carries a 3'
// overhang instead (KpnI-like) — the bottom strand starts |left_d| before a
// (its bases pair with the target's top bases before the left cut) or the
// top strand ends |right_d| before the bottom strand (its tail pairs with
// the target's bottom bases before the right cut).

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
    right_over: Vec<u8>,      // right-end 5' bottom overhang, 5'→3'
    left_under: Vec<u8>,      // left-end 3' overhang bases: the source's top [a+left_d..a), empty unless left_d < 0
}

fn revcomp(s: &[u8]) -> Vec<u8> { s.iter().rev().map(|&b| complement(b)).collect() }

/// `k` bases of `seq` starting at `start`, wrapping a circular sequence.
fn take(seq: &[u8], start: usize, k: usize) -> Vec<u8> {
    let n = seq.len();
    (0..k).map(|i| seq[(start + i) % n]).collect()
}

/// Cut an insert between top-strand cut positions [a, b) using the enzymes
/// that made the cuts; `d` values are each enzyme's bot_cut - top_cut. On a
/// circular source a left 3' overhang reaching before position zero wraps
/// around the origin.
pub fn copy_insert(seq: &[u8], a: usize, b: usize, left_d: i32, right_d: i32, left_enz: &str, right_enz: &str, name: &str, circular: bool) -> Result<Insert, String> {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    if b > seq.len() { return Err("selection past the end of the sequence".into()); }
    // every overhang carried inside the cut range must fit in it: the 5'
    // left overhang (top[a..a+left_d)) and the 3' right overhang (the last
    // |right_d| bases of top)
    if (b - a) < left_d.max(0) as usize { return Err(format!("the {} bp selection is shorter than its {} bp left overhang", b - a, left_d)); }
    if (b - a) < right_d.min(0).unsigned_abs() as usize { return Err(format!("the {} bp selection is shorter than its {} bp right overhang", b - a, right_d.min(0).abs())); }
    // a 5' left overhang with a 3' right one can overlap each other: the
    // bottom strand would start after it ends
    if ((b - a) as i64 + right_d as i64 - left_d as i64) < 0 {
        return Err("the overhangs overlap: the bottom strand would start after it ends".into());
    }
    // a left 3' overhang reaches before a: the bottom strand starts that early
    // (wrapping the origin on a circular source)
    let n = seq.len().max(1);
    let left_under = if left_d < 0 {
        let k = left_d.unsigned_abs() as usize;
        if k > n { return Err("the 3' overhang is longer than the whole sequence".into()); }
        if k > a && !circular { return Err("the left 3' overhang runs past the start of the sequence".into()); }
        (0..k).map(|i| seq[(a + n - k + i) % n]).collect()
    } else { vec![] };
    let top = seq[a..b].to_vec();
    let right_over = if right_d > 0 { revcomp(&take(seq, b, right_d as usize)) } else { vec![] };
    Ok(Insert { name: name.to_string(), top, left_d, right_d, left_enz: left_enz.to_string(), right_enz: right_enz.to_string(), right_over, left_under })
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

fn eq_ci(a: &[u8], b: &[u8]) -> bool { a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.eq_ignore_ascii_case(y)) }

/// Whether an insert ligates into a target whose region [x, y) was cut out:
/// the left enzyme (d = dl) cut at x, the right enzyme (d = dr) at y.
/// 5' junctions: the insert's left top overhang anneals to the target's
/// bases past the left cut, its right bottom overhang to the target's top
/// bases past the right cut — direct pairings. 3' junctions (d < 0): the
/// insert's left bottom overhang anneals to the target's top bases before
/// the left cut, its right top tail to the target's top bases before the
/// right cut — also direct pairings, compared case-insensitively.
pub fn compatible(insert: &Insert, target: &[u8], x: usize, y: usize, dl: i32, dr: i32, circular: bool) -> Result<(), String> {
    let n = target.len().max(1);
    if insert.left_d != dl || insert.right_d != dr {
        return Err(format!("overhang lengths differ: insert {}+{}, target {}+{}",
            insert.left_d, insert.right_d, dl, dr));
    }
    // the target's bottom strand runs [x+dl, y+dr): its span must not invert
    if (y as i64 + dr as i64) - (x as i64 + dl as i64) < 0 {
        return Err("the target's overhangs overlap: its bottom strand would start after it ends".into());
    }
    // 3' overhangs must fit inside the target before wrapping
    if dl < 0 && dl.unsigned_abs() as usize > n { return Err("the 3' overhang is longer than the whole target".into()); }
    if dr < 0 && dr.unsigned_abs() as usize > n { return Err("the 3' overhang is longer than the whole target".into()); }
    // left junction
    if dl > 0 {
        let t = take(target, x, dl as usize);
        if !eq_ci(&insert.top[..dl as usize], &t) {
            return Err("the left overhang bases do not match the target bases past the left cut".into());
        }
    } else if dl < 0 {
        // the insert's bottom overhang pairs with the target's top bases
        // immediately before the cut (wrapping the origin on a circle)
        let k = dl.unsigned_abs() as usize;
        if k > x && !circular { return Err("the target's left 3' overhang runs before its start".into()); }
        let t: Vec<u8> = (0..k).map(|i| target[(x + n - k + i) % n]).collect();
        if !eq_ci(&insert.left_under, &t) {
            return Err("the left overhang bases do not match the target bases before the left cut".into());
        }
    }
    // right junction
    if dr > 0 {
        let t = take(target, y, dr as usize);
        if !eq_ci(&insert.right_over, &revcomp(&t)) {
            return Err("the right overhang bases do not match the target bases past the right cut".into());
        }
    } else if dr < 0 {
        // the insert's top tail pairs with the target's bottom bases —
        // equivalently its top bases — immediately before the cut
        let k = dr.unsigned_abs() as usize;
        if k > y && !circular { return Err("the target's right 3' overhang runs before its start".into()); }
        let t: Vec<u8> = (0..k).map(|i| target[(y + n - k + i) % n]).collect();
        if !eq_ci(&insert.top[insert.top.len() - k..], &t) {
            return Err("the right overhang bases do not match the target bases before the right cut".into());
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
/// Returns the assembled sequence, whether the product is circular, and for
/// every fragment its (out_start, in_start) — where its contributed bases
/// begin in the product and in the fragment — so callers can remap feature
/// locations.
pub fn gibson_assemble(fragments: &[Vec<u8>], min_overlap: usize) -> Result<(Vec<u8>, bool, Vec<(usize, usize)>), String> {
    if fragments.is_empty() { return Err("no fragments".into()); }
    if fragments.len() == 1 { return Ok((fragments[0].clone(), false, vec![(0, 0)])); }
    let mut out = fragments[0].clone();
    let mut spans = vec![(0usize, 0usize)];
    for w in fragments.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        let ov = overlap(a, b).ok_or("two neighbouring fragments do not overlap")?;
        if ov < min_overlap { return Err(format!("overlap of {ov} bp is below the minimum of {min_overlap}")); }
        spans.push((out.len(), ov));
        out.extend_from_slice(&b[ov..]);
    }
    let (last, first) = (fragments.last().unwrap(), &fragments[0]);
    if let Some(close) = overlap(last, first).filter(|k| *k >= min_overlap) {
        let n = out.len();
        out.truncate(n - close);
        return Ok((out, true, spans));
    }
    Ok((out, false, spans))
}

/// Longest overlap where `a` ends with what `b` starts with (case-insensitive,
/// up to 60 bases — the maximum Gibson homology the assembler looks for).
pub fn overlap(a: &[u8], b: &[u8]) -> Option<usize> {
    let max = 60.min(a.len()).min(b.len());
    (1..=max).rev().find(|&k| a.ends_with_ci(&b[..k]))
}

trait EndsCi { fn ends_with_ci(&self, suffix: &[u8]) -> bool; }
impl EndsCi for [u8] {
    fn ends_with_ci(&self, suffix: &[u8]) -> bool {
        self.len() >= suffix.len() && self[self.len() - suffix.len()..].iter().zip(suffix).all(|(x, y)| x.eq_ignore_ascii_case(y))
    }
}

/// Clone Selected Fragment: the subsequence of a span, wrapping the origin.
pub fn clone_fragment(seq: &[u8], a: usize, b: usize) -> Vec<u8> {
    if a <= b { seq[a..b].to_vec() } else { seq[a..].iter().chain(seq[..b].iter()).copied().collect() }
}

/// Shift every coordinate in a GenBank location text ("join(101..110,201..210)",
/// "complement(3..50)", "<101..>110") by `delta`, preserving the location
/// structure. Only pure coordinate numbers shift: range operands move both
/// sides, while words carrying letters — remote accessions such as
/// J12345.1 — are left untouched, as are fuzzy markers (<, >).
pub fn offset_range(range: &str, delta: i64) -> String {
    let mut out = String::new();
    let chars: Vec<char> = range.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        if chars[i].is_ascii_alphanumeric() || chars[i] == '.' {
            // a word runs over alphanumerics and single dots, but the ".."
            // range operator always ends it, so fuzzy bounds like <101..>110
            // parse as two numbers around the operator
            let start = i;
            while i < chars.len() {
                let c = chars[i];
                // a dot never starts a word, and the ".." operator ends one
                if c == '.' && (i == start || (i + 1 < chars.len() && chars[i + 1] == '.')) { break; }
                if !(c.is_ascii_alphanumeric() || c == '.') { break; }
                i += 1;
            }
            if i == start {
                // a lone dot or the first dot of the ".." operator
                out.push(chars[start]);
                i += 1;
                continue;
            }
            let word: String = chars[start..i].iter().collect();
            if let Ok(v) = word.parse::<i64>() {
                out.push_str(&v.saturating_add(delta).to_string());
            } else {
                out.push_str(&word);
            }
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}

/// The segments of a GenBank location as (start, end, fuzzy_start, fuzzy_end)
/// tuples, 1-based inclusive as written; the fuzzy fields carry the exact
/// marker written ("<" or ">", or none for an exact bound). Understands
/// "a..b", fuzzy bounds like "<a" and ">b", join(...) and complement(...)
/// wrappers. Returns nothing when the location holds no "a..b" pair
/// (callers fall back to the feature's numeric bounds).
pub fn parse_location(loc: &str) -> Vec<(i64, i64, Option<char>, Option<char>)> {
    let c: Vec<char> = loc.chars().collect();
    let mut out = vec![];
    let mut i = 0;
    while i < c.len() {
        let fuzzy_s = (c[i] == '<' || c[i] == '>').then_some(c[i]);
        let ds = i + usize::from(fuzzy_s.is_some());
        if ds < c.len() && c[ds].is_ascii_digit() {
            let mut j = ds;
            while j < c.len() && c[j].is_ascii_digit() { j += 1; }
            let a: i64 = c[ds..j].iter().collect::<String>().parse().unwrap_or(0);
            if j + 1 < c.len() && c[j] == '.' && c[j + 1] == '.' {
                let mut k = j + 2;
                let fuzzy_e = if k < c.len() && (c[k] == '<' || c[k] == '>') { Some(c[k]) } else { None };
                let de = k + usize::from(fuzzy_e.is_some());
                if de < c.len() && c[de].is_ascii_digit() {
                    let mut m = de;
                    while m < c.len() && c[m].is_ascii_digit() { m += 1; }
                    let b: i64 = c[de..m].iter().collect::<String>().parse().unwrap_or(0);
                    out.push((a, b, fuzzy_s, fuzzy_e));
                    i = m;
                    continue;
                }
            }
        }
        i += 1;
    }
    out
}

/// The top-strand cut position of a recognition-site hit starting at `p`.
/// A forward-oriented site cuts at `p + top_cut` (bottom at `p + bot_cut`);
/// a reverse-complement hit cuts at `p + len - bot_cut` (bottom at
/// `p + len - top_cut`). On a circular sequence cuts wrap modulo the length;
/// on a linear one both strand cuts must lie within [0, n] for the boundary
/// to exist. `matches` is the IUPAC pattern match (pattern byte, sequence
/// byte) -> bool.
pub fn hit_top_cut(seq: &[u8], p: usize, site: &[u8], top_cut: i32, bot_cut: i32, circular: bool, matches: impl Fn(u8, u8) -> bool) -> Option<usize> {
    let (n, l) = (seq.len(), site.len());
    if l == 0 || n == 0 { return None; }
    if p + l > n && !circular { return None; }
    let forward = (0..l).all(|k| matches(site[k], seq[(p + k) % n]));
    let (top, bot) = if forward {
        (p as i64 + top_cut as i64, p as i64 + bot_cut as i64)
    } else {
        (p as i64 + l as i64 - bot_cut as i64, p as i64 + l as i64 - top_cut as i64)
    };
    if circular { Some(top.rem_euclid(n as i64) as usize) }
    else if (0..=n as i64).contains(&top) && (0..=n as i64).contains(&bot) { Some(top as usize) }
    else { None }
}

#[cfg(test)]
mod tests {
    use super::*;

    // EcoRI G^AATTC: d = +4, 5' AATT overhang on both ends.
    const VECTOR: &[u8] = b"AAAAGAATTCCGGGGGGGGGGGAATTCGGTT";

    #[test]
    fn ecori_geometry() {
        let ins = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "fragment", false).unwrap();
        assert_eq!(ins.top, VECTOR[5..22]);
        let (l, r) = ins.overhangs();
        assert_eq!(l, b"AATT"); // top bases right after the left cut
        assert_eq!(r, b"AATT"); // complement of the source's AATT past the right cut
    }

    #[test]
    fn paste_into_same_enzyme_cuts() {
        let ins = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "insert", false).unwrap();
        // a target opened with EcoRI on both sides: TTGAATTC...GAATTCCT
        let target = b"TTTTGAATTCAAAAAAAAGAATTCCTTTT";
        // top cuts at 5 and 19
        compatible(&ins, target, 5, 19, 4, 4, false).unwrap();
        let out = paste_insert(target, &ins, 5, 19);
        assert_eq!(&out[..5], &target[..5]);
        assert_eq!(&out[5..5 + ins.top.len()], &ins.top);
        assert_eq!(&out[5 + ins.top.len()..], &target[19..]);
        // the left junction reconstitutes the EcoRI site: ...G + AATTC...
        assert_eq!(&out[4..10], b"GAATTC");
    }

    #[test]
    fn mixed_enzyme_overhangs_are_rejected() {
        // an EcoRI-left / BamHI-right insert (overhang AATT ... GATC on the
        // bottom strand) does not ligate into an EcoRI / EcoRI target
        let src = b"AAAAGAATTCCGGGGGGGGGGGGATCCTT";
        let ins = copy_insert(src, 5, 22, 4, 4, "EcoRI", "BamHI", "mixed", false).unwrap();
        let (_, r) = ins.overhangs();
        assert_eq!(r, b"GATC"); // revcomp of GGATCC's tail past the cut
        let target = b"TTTTGAATTCAAAAAAAAGAATTCCTTTT";
        assert!(compatible(&ins, target, 5, 19, 4, 4, false).is_err());
        // and into a matching EcoRI / BamHI target it does
        let t2 = b"TTTTGAATTCAAAAAAAAGGATCCCTTTT";
        compatible(&ins, t2, 5, 19, 4, 4, false).unwrap();
    }

    #[test]
    fn overhang_comparison_is_case_insensitive() {
        let ins = Insert { name: "s".into(), top: b"aattCGCGCG".to_vec(), left_d: 4, right_d: 0, left_enz: "EcoRI".into(), right_enz: "SmaI".into(), left_under: vec![], right_over: vec![] };
        let target = b"TTTTGAATTCAAAAAAAAAAAAAACTTTT";
        compatible(&ins, target, 5, 9, 4, 0, false).unwrap();
    }

    #[test]
    fn selection_shorter_than_overhang_rejected() {
        // MspI C^CGG (d=2) and SmaI CCC^GGG (d=0) bound a 1-base top strand
        assert!(copy_insert(b"CCCGGG", 2, 3, 2, 0, "MspI", "SmaI", "t", false).is_err());
    }

    #[test]
    fn blunt_and_incompatible() {
        let blunt = Insert { name: "b".into(), top: b"ACGTACGTAC".to_vec(), left_d: 0, right_d: 0, left_enz: "SmaI".into(), right_enz: "SmaI".into(), left_under: vec![], right_over: vec![] };
        compatible(&blunt, b"AAAAAAAAAA", 3, 7, 0, 0, false).unwrap();
        // an AATT insert does not ligate into a blunt target
        let sticky = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "s", false).unwrap();
        assert!(compatible(&sticky, b"AAAAAAAAAAA", 3, 7, 0, 0, false).is_err());
        // and not into a target whose overhang bases differ
        let t = b"TTTTGACTAGAAAAAAAAGACTAGCTTTT";
        assert!(compatible(&sticky, t, 5, 17, 4, 4, false).is_err());
    }

    #[test]
    fn circular_target_wraps() {
        // cut at the origin: target[y..] then [x..y) spans the wrap
        let ins = copy_insert(VECTOR, 5, 22, 4, 4, "EcoRI", "EcoRI", "i", false).unwrap();
        let target = b"GAATTCAAAAAAAGAATTCGG";
        // left cut EcoRI at 1, right cut EcoRI at 14: target[..1] ++ ins ++ target[14..]
        compatible(&ins, target, 1, 14, 4, 4, false).unwrap();
        let out = paste_insert(target, &ins, 1, 14);
        assert_eq!(out.len(), 1 + ins.top.len() + (target.len() - 14));
    }

    #[test]
    fn gibson_two_fragments() {
        let a = b"ACGTACGTACGGGGGGGGGGTTTTTTTTTT".to_vec();
        let b = b"TTTTTTTTTTCCCCCCCCCC".to_vec();
        let (out, circular, _) = gibson_assemble(&[a, b], 10).unwrap();
        assert_eq!(out, b"ACGTACGTACGGGGGGGGGGTTTTTTTTTTCCCCCCCCCC");
        assert!(!circular, "no last-to-first overlap: the product is linear");
    }

    #[test]
    fn gibson_overlap_is_case_insensitive() {
        let a = b"ACGTACGTACGGGGGGGGGGtttttttttt".to_vec(); // GenBank-style lowercase
        let b = b"TTTTTTTTTTCCCCCCCCCC".to_vec();
        let (out, _, _) = gibson_assemble(&[a, b], 10).unwrap();
        assert_eq!(out, b"ACGTACGTACGGGGGGGGGGttttttttttCCCCCCCCCC");
    }

    #[test]
    fn gibson_three_fragments_circular() {
        let a = b"AAAAAAAAAAGGGGGG".to_vec();
        let b = b"GGGGGGCCCCCC".to_vec();
        let c = b"CCCCCCTTTTTTAAAA".to_vec();
        let (out, circular, _) = gibson_assemble(&[a, b, c], 4).unwrap();
        // the last fragment's AAAA tail overlaps the first fragment's head
        assert_eq!(out, b"AAAAAAAAAAGGGGGGCCCCCCTTTTTT");
        assert!(circular);
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
    fn offset_range_shifts_numbers_and_keeps_structure() {
        assert_eq!(offset_range("101..110", -100), "1..10");
        assert_eq!(offset_range("complement(3..50)", -2), "complement(1..48)");
        assert_eq!(offset_range("join(101..110,201..210)", -100), "join(1..10,101..110)");
    }

    #[test]
    fn offset_range_keeps_remote_accessions() {
        // remote references carry their own accession numbers, which must not shift
        assert_eq!(offset_range("J12345.1:20..30", 10), "J12345.1:30..40");
        assert_eq!(offset_range("join(5..9,U00004.5:12..19)", 100), "join(105..109,U00004.5:112..119)");
        assert_eq!(offset_range("U00004.5", 1000), "U00004.5");
    }

    #[test]
    fn offset_range_shifts_fuzzy_bounds() {
        assert_eq!(offset_range("<101..>110", 10), "<111..>120");
        assert_eq!(offset_range("join(<1..5,>900..950)", -1), "join(<0..4,>899..949)");
    }

    #[test]
    fn parse_location_reads_segments_and_markers() {
        assert_eq!(parse_location("101..110"), vec![(101, 110, None, None)]);
        assert_eq!(parse_location("<101..>110"), vec![(101, 110, Some('<'), Some('>'))]);
        assert_eq!(parse_location(">101..<110"), vec![(101, 110, Some('>'), Some('<'))]);
        assert_eq!(parse_location("complement(join(5..9,201..210))"), vec![(5, 9, None, None), (201, 210, None, None)]);
        assert_eq!(parse_location("J12345.1:20..30"), vec![(20, 30, None, None)]);
        assert!(parse_location("complement(56)").is_empty());
    }

    #[test]
    fn gibson_spans_report_fragment_positions() {
        let (out, circular, spans) = gibson_assemble(&[b"AAAACCCC".to_vec(), b"CCCCGGGG".to_vec()], 4).unwrap();
        assert_eq!(out, b"AAAACCCCGGGG");
        assert!(!circular);
        // fragment 0 contributes all of itself from product position 0;
        // fragment 1 contributes from its own base 4 (the overlap) at
        // product position 8
        assert_eq!(spans, vec![(0, 0), (8, 4)]);

        // a closed circle: the second fragment's tail is the product's start
        let (out, circular, spans) = gibson_assemble(&[b"TTTTAAAA".to_vec(), b"AAAATTTT".to_vec()], 4).unwrap();
        assert_eq!(out, b"TTTTAAAA");
        assert!(circular);
        assert_eq!(spans, vec![(0, 0), (8, 4)]);
    }

    #[test]
    fn hit_top_cut_forward_and_reverse() {
        // BsaI GGTCTC(1/5): top_cut 7, bot_cut 11
        let eq = |a: u8, b: u8| a == b;
        // forward site GGTCTC at 4 cuts the top strand at 4 + 7 = 11
        let fwd = b"AAAAGGTCTCAAAAAAAAAA";
        assert_eq!(hit_top_cut(fwd, 4, b"GGTCTC", 7, 11, false, eq), Some(11));
        // reverse hit: revcomp(GGTCTC) = GAGACC at 5 cuts at 5 + 6 - 11 = 0
        let rev = b"AAAAAGAGACCAAAAAAAA";
        assert_eq!(hit_top_cut(rev, 5, b"GGTCTC", 7, 11, false, eq), Some(0));
    }

    #[test]
    fn hit_top_cut_circular_wraps_and_linear_rejects() {
        let eq = |a: u8, b: u8| a == b;
        // circular: a cut landing on len is boundary 0
        let seq = b"TTGAATTC";
        // EcoRI GAATTC at 2: forward, top cut 3; on the 8-base circle a site at
        // the end whose cut runs past len wraps
        assert_eq!(hit_top_cut(seq, 2, b"GAATTC", 1, 5, true, eq), Some(3));
        // a BsaI site late in a circular sequence cuts before the origin
        let late = b"AAAAAAAAAGGTCTC";
        assert_eq!(hit_top_cut(late, 9, b"GGTCTC", 7, 11, true, eq), Some(16 % 15));
        // linear: the same cut does not exist
        assert_eq!(hit_top_cut(late, 9, b"GGTCTC", 7, 11, false, eq), None);
    }

    #[test]
    fn hit_top_cut_linear_requires_both_strand_cuts() {
        let eq = |a: u8, b: u8| a == b;
        // forward BsaI at 4 in a 14 bp molecule: top cut 11 is inside, but the
        // bottom cut 15 is past the end — the restriction end does not exist
        let seq = b"AAAAGGTCTCAAAA";
        assert_eq!(seq.len(), 14);
        assert_eq!(hit_top_cut(seq, 4, b"GGTCTC", 7, 11, false, eq), None);
        // with two more bases both cuts exist and the top cut is 11
        let seq = b"AAAAGGTCTCAAAAAA";
        assert_eq!(hit_top_cut(seq, 4, b"GGTCTC", 7, 11, false, eq), Some(11));
    }

    #[test]
    fn gibson_finds_homology_longer_than_half_a_fragment() {
        // 16 bp of terminal homology in a 20 bp fragment
        let a = b"ACGTAAAAAAAAAAAAAAAA".to_vec();
        let b = b"AAAAAAAAAAAAAAAAGGG".to_vec();
        let (out, _, _) = gibson_assemble(&[a, b], 12).unwrap();
        assert_eq!(out, b"ACGTAAAAAAAAAAAAAAAAGGG");
    }

    #[test]
    fn three_prime_overhangs_clone() {
        // KpnI-like left end (d = -4): the fragment's bottom strand starts
        // four bases before the top cut, so it carries a left 3' overhang
        let ins = copy_insert(VECTOR, 7, 20, -4, 0, "KpnI", "SmaI", "x", false).unwrap();
        assert_eq!(ins.top, VECTOR[7..20]);
        assert_eq!(ins.left_under, VECTOR[3..7]);
        // matching target: same bases before its left cut
        assert!(compatible(&ins, VECTOR, 7, 20, -4, 0, false).is_ok());
        // a target whose bases before the cut differ does not ligate
        let other = b"TTTTGAATTCCGGGGGGGGGGGAATTCGGTT";
        assert!(compatible(&ins, other, 7, 20, -4, 0, false).is_err());
        // the pasted product keeps the target's flanks and the insert
        assert_eq!(paste_insert(VECTOR, &ins, 7, 20), [VECTOR[..7].to_vec(), ins.top.clone(), VECTOR[20..].to_vec()].concat());

        // PstI-like right end (d = -4): the fragment's top tail pairs with
        // the target's bases before its right cut
        let ins = copy_insert(VECTOR, 5, 20, 0, -4, "SmaI", "PstI", "y", false).unwrap();
        assert!(compatible(&ins, VECTOR, 5, 20, 0, -4, false).is_ok());
        let mismatched = [VECTOR[..16].to_vec(), b"TTTT".to_vec(), VECTOR[20..].to_vec()].concat();
        assert!(compatible(&ins, &mismatched, 5, 20, 0, -4, false).is_err());

        // a left 3' overhang reaching before the sequence start is refused
        assert!(copy_insert(VECTOR, 2, 10, -4, 0, "KpnI", "SmaI", "z", false).is_err());
    }
}
