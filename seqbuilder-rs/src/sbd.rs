// SeqBuilder Pro document (.sbd) reader.
//
// Container (big-endian): u32 at 0x2C -> object table (0x40-byte records: idx, 1, type4, subtype4, 0,
// group offset, Pascal name). Group: u32 0xCC, count, count, 0, first entry offset. Entries are 32 bytes:
// type4, id, flags, len, len, timestamp, 0, offset. Chunks used here:
//   UVSQ  sequence, ASCII
//   UVTX  GenBank header text (LOCUS..COMMENT), CRLF; "^^" when the sequence was made in SeqBuilder
//   UVSI  sequence info; its trailing property block holds SeqShape (1 linear, 2 circular)
//   UVFT  one feature per id: 416-byte record + CR-separated "/key=value" qualifiers
//   UVTL  layout as a CommonBlockFmt property tree (8-byte key, u32 LE length, value), holding styles:
//         ALLDSTAT.CMNSTATE.VSSTYLE.OBJINFO.STYLE (per object) and .TYPEINFO.TDATA{TINFO,STYLE} (per feature type)
use std::collections::HashMap;

pub struct SbdFeat {
    pub id: u32,
    pub name: String,
    pub kind: String,
    pub start: usize, // 1-based, as stored; start > end on the complement strand
    pub end: usize,
    pub complement: bool,
    pub quals: Vec<(String, String)>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Style {
    pub line: Option<[u8; 3]>,
    pub weight: f64,
    pub fill: Option<[u8; 3]>,
}

pub struct Sbd {
    pub seq: Vec<u8>,
    pub header: String,
    pub circular: bool,
    pub feats: Vec<SbdFeat>,
    pub type_styles: HashMap<String, Style>,
    pub obj_styles: HashMap<u32, Style>, // by UVFT id
}

const FEATURE_OBJECT: u32 = 9; // object type of features in STYLE records

fn u32be(d: &[u8], o: usize) -> Option<u32> {
    Some(u32::from_be_bytes(d.get(o..o + 4)?.try_into().ok()?))
}

fn u32le(d: &[u8], o: usize) -> Option<u32> {
    Some(u32::from_le_bytes(d.get(o..o + 4)?.try_into().ok()?))
}

fn pascal(d: &[u8], o: usize) -> String {
    let n = d.get(o).copied().unwrap_or(0) as usize;
    d.get(o + 1..o + 1 + n).map(|s| s.iter().map(|&b| b as char).collect()).unwrap_or_default()
}

/// All chunks as (type, id, bytes).
fn chunks(d: &[u8]) -> Option<Vec<([u8; 4], u32, &[u8])>> {
    let mut out = vec![];
    let mut o = u32be(d, 0x2c)? as usize;
    while o + 0x40 <= d.len() {
        let grp = u32be(d, o + 20)? as usize;
        let (n, first) = (u32be(d, grp + 4)? as usize, u32be(d, grp + 16)? as usize);
        for k in 0..n {
            let e = first + 32 * k;
            let typ: [u8; 4] = d.get(e..e + 4)?.try_into().ok()?;
            let (id, len, off) = (u32be(d, e + 4)?, u32be(d, e + 12)? as usize, u32be(d, e + 28)? as usize);
            out.push((typ, id, d.get(off..off + len)?));
        }
        o += 0x40;
    }
    Some(out)
}

/// Children of a CommonBlockFmt property block.
fn props(b: &[u8]) -> Vec<(&str, &[u8])> {
    let mut i = if b.starts_with(b"CommonBlockFmt\0") { 17 } else { 0 };
    let mut out = vec![];
    while let (Some(key), Some(len)) = (b.get(i..i + 8), u32le(b, i + 8)) {
        let Some(val) = b.get(i + 12..i + 12 + len as usize) else { break };
        let key = std::str::from_utf8(key).unwrap_or("").trim_end_matches('\0');
        out.push((key, val));
        i += 12 + len as usize;
    }
    out
}

fn child<'a>(b: &'a [u8], key: &str) -> Option<&'a [u8]> {
    props(b).into_iter().find(|(k, _)| *k == key).map(|(_, v)| v)
}

// STYLE (35 bytes, LE): u32 object index, u32 object type, u8 line on, rgb line, f64 weight,
// u8 fill on, rgb fill, u8 shadow on, rgb shadow, 7 bytes not decoded yet (shape, linkage, label position).
fn style(s: &[u8]) -> Option<(u32, u32, Style)> {
    let rgb = |o: usize| -> [u8; 3] { [s[o], s[o + 1], s[o + 2]] };
    if s.len() < 24 { return None; }
    Some((u32le(s, 0)?, u32le(s, 4)?, Style {
        line: (s[8] != 0).then(|| rgb(9)),
        weight: f64::from_le_bytes(s[12..20].try_into().ok()?),
        fill: (s[20] != 0).then(|| rgb(21)),
    }))
}

impl Sbd {
    /// A GenBank header field (DEFINITION, ACCESSION, COMMENT...) from UVTX, continuation lines included.
    pub fn header_field(&self, key: &str, sep: &str) -> String {
        let mut lines = self.header.lines().skip_while(|l| l.split_whitespace().next() != Some(key));
        let Some(first) = lines.next() else { return String::new() };
        std::iter::once(first[key.len()..].trim())
            .chain(lines.take_while(|l| l.starts_with("            ")).map(str::trim))
            .collect::<Vec<_>>().join(sep)
    }
}

pub fn read(d: &[u8]) -> Result<Sbd, String> {
    let chunks = chunks(d).ok_or("not a SeqBuilder document (bad chunk directory)")?;
    let get = |t: [u8; 4]| chunks.iter().filter(move |c| c.0 == t);
    let text = |b: &[u8]| -> String { b.iter().map(|&c| c as char).collect() }; // Latin-1
    let seq = get(*b"UVSQ").next().ok_or("no sequence")?.2.to_vec();
    let header = get(*b"UVTX").next().map(|c| text(c.2).replace("\r\n", "\n")).unwrap_or_default();
    // UVSI ends with a property block whose SeqShape is 1 (linear) or 2 (circular); fall back to LOCUS.
    let shape = get(*b"UVSI").next().and_then(|c| {
        let at = c.2.windows(15).position(|w| w == b"CommonBlockFmt\0")?;
        child(&c.2[at..], "SeqShape").and_then(|v| u32le(v, 0))
    });
    let circular = shape.map_or_else(|| header.lines().next().is_some_and(|l| l.starts_with("LOCUS") && l.contains("circular")), |s| s == 2);

    let feats = get(*b"UVFT").filter_map(|&(_, id, r)| {
        let (start, end) = (u32be(r, 394)? as usize, u32be(r, 398)? as usize);
        let quals = text(r.get(416..)?).split('\r').filter_map(|l| {
            let (k, v) = l.strip_prefix('/')?.split_once('=').unwrap_or((l.trim_start_matches('/'), ""));
            Some((k.to_string(), v.trim_matches('"').to_string()))
        }).collect();
        Some(SbdFeat {
            id, name: pascal(r, 6), kind: pascal(r, 266), start, end,
            complement: u16::from_be_bytes([r[402], r[403]]) != 0, quals,
        })
    }).collect();

    let (mut type_styles, mut obj_styles) = (HashMap::new(), HashMap::new());
    let vs = get(*b"UVTL").next()
        .and_then(|c| child(c.2, "ALLDSTAT")).and_then(|b| child(b, "CMNSTATE")).and_then(|b| child(b, "VSSTYLE"));
    if let Some(vs) = vs {
        for (k, v) in props(vs) {
            match k {
                "OBJINFO" => for (_, s) in props(v).into_iter().filter(|(k, _)| *k == "STYLE") {
                    if let Some((idx, FEATURE_OBJECT, st)) = style(s) { obj_styles.insert(idx, st); }
                },
                "TYPEINFO" => for (_, td) in props(v).into_iter().filter(|(k, _)| *k == "TDATA") {
                    let (Some(ti), Some(st)) = (child(td, "TINFO"), child(td, "STYLE").and_then(style)) else { continue };
                    let n = u32le(ti, 0).unwrap_or(0) as usize;
                    if let Some(name) = ti.get(4..4 + n) { type_styles.insert(text(name), st.2); }
                },
                _ => {}
            }
        }
    }
    Ok(Sbd { seq, header, circular, feats, type_styles, obj_styles })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_pbr322_and_custom_colour() {
        let d = read(&std::fs::read("testdata/pBR322.sbd").unwrap()).unwrap();
        assert_eq!(d.seq.len(), 4361);
        assert!(d.circular);
        assert_eq!(d.feats.len(), 50);
        let p1 = &d.feats[3];
        assert_eq!((p1.kind.as_str(), p1.start, p1.end, p1.complement), ("regulatory", 33, 27, true));
        assert!(p1.quals.contains(&("note".into(), "promoter P1 (6)".into())));
        assert_eq!(d.type_styles["CDS"].fill, Some([0xe8, 0xfa, 0xfd]));
        assert!(!d.obj_styles.contains_key(&38));

        let lin = read(&std::fs::read("testdata/ALDH2.sbd").unwrap()).unwrap();
        assert_eq!((lin.seq.len(), lin.circular, lin.feats.len()), (114842, false, 2));
        let (fwd, rev) = (&lin.feats[0], &lin.feats[1]);
        assert_eq!((fwd.kind.as_str(), fwd.start, fwd.end, fwd.complement), ("PCR_primer", 5087, 5124, false));
        assert_eq!((rev.name.as_str(), rev.start, rev.end, rev.complement), ("Reverse (5087..5124)", 5124, 5087, true));
        assert_eq!(lin.header_field("DEFINITION", " "), "");

        let red = read(&std::fs::read("testdata/pBR322_blaRED.sbd").unwrap()).unwrap();
        assert_eq!((red.feats[38].name.as_str(), red.feats[38].kind.as_str()), ("bla", "CDS"));
        assert_eq!(red.obj_styles[&38].fill, Some([0xff, 0, 0]));
    }
}
