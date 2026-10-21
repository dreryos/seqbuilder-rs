// SeqBuilder Pro document (.sbd) reader and writer.
//
// Container (big-endian): u32 at 0x2C -> object table (0x40-byte records: idx, 1, type4, subtype4, 0,
// group offset, Pascal name). Group: u32 0xCC, count, count, 0, first entry offset. Entries are 32 bytes:
// type4, id, flags, len, len, timestamp, 0, offset. Chunks used here:
//   UVSQ  sequence, ASCII
//   UVTX  GenBank header text (LOCUS..COMMENT), CRLF; "^^" when the sequence was made in SeqBuilder
//   UVSI  sequence info: u32 BE sequence length at 16; its trailing property block holds SeqShape
//         (1 linear, 2 circular)
//   UVFT  one feature per id: 416-byte record + CR-separated "/key=value" qualifiers. Record: u16 0x65,
//         Pascal name at 6 (256 bytes), u32 416 at 262, Pascal type at 266, start/end (u32, 1-based,
//         reversed on the complement strand) at 394/398, u16 complement at 402, u16 1 at 404; the
//         remaining bytes are the same in every record seen. FXIF id N+1 belongs to UVFT id N.
//   UVTL  layout as a CommonBlockFmt property tree (8-byte key, u32 LE length, value; some values start
//         with a few bytes of their own before their items). Styles live in VSSTYLE blocks: OBJINFO.STYLE
//         per object, TYPEINFO.TDATA{TINFO,STYLE} per feature type. There is one VSSTYLE per view state:
//         ALLDSTAT, the open view (InfSpltF...SeqView) and every default layout in DefLayLs; SeqBuilder
//         draws from the view copies, so an edit has to reach all of them.
// Writing keeps every chunk as read and replaces only the ones that were edited, then lays the file
// out again: header, chunk data, each object's entries followed by its group header, object table.
use std::collections::HashMap;

pub struct SbdFeat {
    pub id: u32,
    pub name: String,
    pub kind: String,
    pub start: usize, // 1-based, as stored; start > end on the complement strand
    pub end: usize,
    pub complement: bool,
    pub quals: Vec<String>, // "key=value" lines as stored, value quoting included
}

impl SbdFeat {
    /// A qualifier's value without its quotes.
    pub fn qual(&self, key: &str) -> Option<String> {
        self.quals.iter().find_map(|l| l.strip_prefix(key)?.strip_prefix('=').map(|v| v.trim_matches('"').to_string()))
    }

    fn encode(&self) -> Vec<u8> {
        let latin1 = |s: &str, max: usize| -> Vec<u8> {
            let b: Vec<u8> = s.chars().map(|c| u8::try_from(c as u32).unwrap_or(b'?')).take(max).collect();
            std::iter::once(b.len() as u8).chain(b).collect()
        };
        let mut r = vec![0u8; 416];
        r[1] = 0x65;
        let name = latin1(&self.name, 255);
        r[6..6 + name.len()].copy_from_slice(&name);
        r[262..266].copy_from_slice(&416u32.to_be_bytes());
        let kind = latin1(&self.kind, 60);
        r[266..266 + kind.len()].copy_from_slice(&kind);
        (r[327], r[331], r[353]) = (1, 1, 1);
        r[394..398].copy_from_slice(&(self.start as u32).to_be_bytes());
        r[398..402].copy_from_slice(&(self.end as u32).to_be_bytes());
        r[402..404].copy_from_slice(&(self.complement as u16).to_be_bytes());
        r[404..406].copy_from_slice(&1u16.to_be_bytes());
        for q in &self.quals {
            r.push(b'/');
            r.extend(latin1(q, usize::MAX).into_iter().skip(1));
            r.push(b'\r');
        }
        r
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Style {
    pub line: Option<[u8; 3]>,
    pub weight: f64,
    pub fill: Option<[u8; 3]>,
}

pub struct Sbd {
    container: Container,
    pub seq: Vec<u8>,
    pub header: String,
    pub circular: bool,
    pub feats: Vec<SbdFeat>,
    pub type_styles: HashMap<String, Style>,
    pub obj_styles: HashMap<u32, Style>, // by UVFT id
    pub shown: HashMap<u32, bool>,       // Show checkbox by UVFT id (VSFEAT.STATE)
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

struct Entry {
    meta: [u8; 32], // as read; length and offset are rewritten on save
    data: Vec<u8>,
    rank: (u64, u32), // place of the data in the file: (original offset, n-th chunk added after it)
}

impl Entry {
    fn is(&self, typ: &[u8; 4]) -> bool {
        &self.meta[..4] == typ
    }

    fn id(&self) -> u32 {
        u32::from_be_bytes(self.meta[4..8].try_into().unwrap())
    }
}

struct Object {
    record: [u8; 0x40],
    group: [u8; 20],
    entries: Vec<Entry>,
}

struct Container {
    head: [u8; 0x30],
    objects: Vec<Object>,
}

impl Container {
    fn parse(d: &[u8]) -> Option<Container> {
        let mut objects = vec![];
        let mut o = u32be(d, 0x2c)? as usize;
        while o + 0x40 <= d.len() {
            let grp = u32be(d, o + 20)? as usize;
            let (n, first) = (u32be(d, grp + 4)? as usize, u32be(d, grp + 16)? as usize);
            let mut entries = vec![];
            for k in 0..n {
                let e = first + 32 * k;
                let meta: [u8; 32] = d.get(e..e + 32)?.try_into().ok()?;
                let (len, off) = (u32be(&meta, 12)? as usize, u32be(&meta, 28)? as usize);
                entries.push(Entry { meta, data: d.get(off..off + len)?.to_vec(), rank: (off as u64, 0) });
            }
            objects.push(Object { record: d[o..o + 0x40].try_into().ok()?, group: d.get(grp..grp + 20)?.try_into().ok()?, entries });
            o += 0x40;
        }
        Some(Container { head: d.get(..0x30)?.try_into().ok()?, objects })
    }

    fn write(&self) -> Vec<u8> {
        let put = |b: &mut [u8], at: usize, v: usize| b[at..at + 4].copy_from_slice(&(v as u32).to_be_bytes());
        let mut out = self.head.to_vec();
        let mut offsets: Vec<Vec<usize>> = self.objects.iter().map(|o| vec![0; o.entries.len()]).collect();
        let mut order: Vec<((u64, u32), usize, usize)> = self.objects.iter().enumerate()
            .flat_map(|(o, obj)| obj.entries.iter().enumerate().map(move |(e, en)| (en.rank, o, e))).collect();
        order.sort();
        for (_, o, e) in order {
            offsets[o][e] = out.len();
            out.extend(&self.objects[o].entries[e].data);
        }
        let mut groups = vec![];
        for (o, obj) in self.objects.iter().enumerate() {
            let first = out.len();
            for (e, entry) in obj.entries.iter().enumerate() {
                let mut meta = entry.meta;
                put(&mut meta, 12, entry.data.len());
                put(&mut meta, 16, entry.data.len());
                put(&mut meta, 28, offsets[o][e]);
                out.extend(meta);
            }
            let mut group = obj.group;
            put(&mut group, 4, obj.entries.len());
            put(&mut group, 8, obj.entries.len());
            put(&mut group, 16, first);
            groups.push(out.len());
            out.extend(group);
        }
        let table = out.len();
        for (obj, grp) in self.objects.iter().zip(groups) {
            let mut record = obj.record;
            put(&mut record, 20, grp);
            out.extend(record);
        }
        put(&mut out, 0x2c, table);
        out
    }

    fn chunks(&self) -> impl Iterator<Item = ([u8; 4], u32, &[u8])> {
        self.objects.iter().flat_map(|o| &o.entries).map(|e| {
            (e.meta[..4].try_into().unwrap(), u32::from_be_bytes(e.meta[4..8].try_into().unwrap()), e.data.as_slice())
        })
    }

    fn chunk_mut(&mut self, typ: &[u8; 4]) -> Option<&mut Vec<u8>> {
        self.chunk_id_mut(typ, None)
    }

    fn chunk_id_mut(&mut self, typ: &[u8; 4], id: Option<u32>) -> Option<&mut Vec<u8>> {
        self.objects.iter_mut().flat_map(|o| &mut o.entries)
            .find(|e| e.is(typ) && id.is_none_or(|id| e.id() == id))
            .map(|e| &mut e.data)
    }

    /// Add a chunk right after the last one of its type, in the directory and in the data (as SeqBuilder does).
    fn add_after_last(&mut self, typ: &[u8; 4], id: u32, data: Vec<u8>) -> Option<()> {
        let obj = self.objects.iter_mut().find(|o| o.entries.iter().any(|e| e.is(typ)))?;
        let at = obj.entries.iter().rposition(|e| e.is(typ))?;
        let mut meta = obj.entries[at].meta;
        meta[4..8].copy_from_slice(&id.to_be_bytes());
        let rank = (obj.entries[at].rank.0, obj.entries[at].rank.1 + 1);
        obj.entries.insert(at + 1, Entry { meta, data, rank });
        Some(())
    }

    /// Replace every chunk of a type with `chunks`, where the first old one was.
    fn replace_all(&mut self, typ: &[u8; 4], chunks: Vec<(u32, Vec<u8>)>) -> Option<()> {
        let obj = self.objects.iter_mut().find(|o| o.entries.iter().any(|e| e.is(typ)))?;
        let at = obj.entries.iter().position(|e| e.is(typ))?;
        let (meta, rank) = (obj.entries[at].meta, obj.entries[at].rank);
        obj.entries.retain(|e| !e.is(typ));
        let new = chunks.into_iter().enumerate().map(|(k, (id, data))| {
            let mut meta = meta;
            meta[4..8].copy_from_slice(&id.to_be_bytes());
            Entry { meta, data, rank: (rank.0, rank.1 + k as u32) }
        });
        obj.entries.splice(at..at, new);
        Some(())
    }

    fn remove(&mut self, typ: &[u8; 4], id: u32) {
        for o in &mut self.objects { o.entries.retain(|e| !(e.is(typ) && e.id() == id)); }
    }

    fn ids(&self, typ: &[u8; 4]) -> Vec<u32> {
        self.objects.iter().flat_map(|o| &o.entries).filter(|e| e.is(typ)).map(Entry::id).collect()
    }
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

/// A property block rebuilt from (key, value) items, keeping the CommonBlockFmt prefix if `like` has one.
fn block<'a>(like: &[u8], items: impl IntoIterator<Item = (&'a str, Vec<u8>)>) -> Vec<u8> {
    let mut out = if like.starts_with(b"CommonBlockFmt\0") { like[..17].to_vec() } else { vec![] };
    for (k, v) in items {
        let mut key = [0u8; 8];
        key[..k.len()].copy_from_slice(k.as_bytes());
        out.extend(key);
        out.extend((v.len() as u32).to_le_bytes());
        out.extend(v);
    }
    out
}

fn valid_key(k: &[u8]) -> bool {
    let n = k.iter().position(|&b| b == 0).unwrap_or(k.len());
    n >= 3 && k[0].is_ascii_alphabetic() && k[..n].iter().all(|b| b.is_ascii_alphanumeric() || *b == b'_') && k[n..].iter().all(|&b| b == 0)
}

/// Headers (offsets of their 8-byte keys) of the items enclosing position `p`, outermost first.
/// A value may start with up to 16 bytes of its own before its items.
fn enclosing(b: &[u8], p: usize) -> Vec<usize> {
    let (mut s, mut e, mut chain) = (17, b.len(), vec![]);
    'down: loop {
        for skip in 0..=16 {
            let mut i = s + skip;
            while let (Some(key), Some(len)) = (b.get(i..i + 8), u32le(b, i + 8)) {
                let end = i + 12 + len as usize;
                if !valid_key(key) || end > e { break; }
                if (i + 12..end).contains(&p) {
                    chain.push(i);
                    (s, e) = (i + 12, end);
                    continue 'down;
                }
                i = end;
            }
        }
        return chain;
    }
}

/// Offsets of every item header named `key` in a UVTL block.
fn find_items(b: &[u8], key: &str) -> Vec<usize> {
    let mut k = [0u8; 8];
    k[..key.len()].copy_from_slice(key.as_bytes());
    b.windows(8).enumerate().filter(|(_, w)| *w == k).map(|(i, _)| i)
        .filter(|&i| enclosing(b, i + 12).last() == Some(&i)).collect()
}

/// Replace the value of every `key` item in a UVTL block, fixing the lengths of all enclosing items.
fn edit_all(tl: &mut Vec<u8>, key: &str, mut f: impl FnMut(&[u8]) -> Vec<u8>) {
    for at in find_items(tl, key).into_iter().rev() {
        let chain = enclosing(tl, at + 12);
        let len = u32le(tl, at + 8).unwrap_or(0) as usize;
        let new = f(&tl[at + 12..at + 12 + len]);
        let delta = new.len() as i64 - len as i64;
        tl.splice(at + 12..at + 12 + len, new);
        for h in chain {
            let l = u32le(tl, h + 8).unwrap_or(0) as i64 + delta;
            tl[h + 8..h + 12].copy_from_slice(&(l as u32).to_le_bytes());
        }
    }
}

// VSFEAT.STATE: u32 count, then per feature id ever made (u32 id, u8 shown). Deleting a feature leaves its entry.
fn feature_state(st: &[u8]) -> Vec<(u32, bool)> {
    let n = u32le(st, 0).unwrap_or(0) as usize;
    (0..n).filter_map(|k| Some((u32le(st, 4 + 5 * k)?, *st.get(8 + 5 * k)? != 0))).collect()
}

fn encode_state(rows: &[(u32, bool)]) -> Vec<u8> {
    let mut v = (rows.len() as u32).to_le_bytes().to_vec();
    for &(id, shown) in rows {
        v.extend(id.to_le_bytes());
        v.push(shown as u8);
    }
    v
}

/// Apply `f` to the feature state of every VSFEAT copy.
fn edit_feature_state(tl: &mut Vec<u8>, mut f: impl FnMut(&mut Vec<(u32, bool)>)) {
    edit_all(tl, "VSFEAT", |vf| {
        block(vf, props(vf).into_iter().map(|(k, v)| {
            if k != "STATE" { return (k, v.to_vec()); }
            let mut rows = feature_state(v);
            f(&mut rows);
            (k, encode_state(&rows))
        }))
    });
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

/// One VSSTYLE block with feature `id`'s object style set or removed.
fn restyle(vs: &[u8], id: u32, kind: &str, st: Option<Style>) -> Vec<u8> {
    let is_kind = |td: &[u8]| child(td, "TINFO").is_some_and(|ti| {
        let n = u32le(ti, 0).unwrap_or(0) as usize;
        ti.get(4..4 + n) == Some(kind.as_bytes())
    });
    let type_style = props(vs).into_iter().find(|(k, _)| *k == "TYPEINFO")
        .and_then(|(_, ti)| props(ti).into_iter().find(|(k, td)| *k == "TDATA" && is_kind(td)))
        .and_then(|(_, td)| child(td, "STYLE").map(<[u8]>::to_vec));
    let template = type_style.unwrap_or_else(|| {
        let mut r = vec![0u8; 35];
        r[12..20].copy_from_slice(&1f64.to_le_bytes());
        r
    });
    let objinfo = |oi: &[u8]| -> Vec<u8> {
        let mine = |v: &[u8]| matches!(style(v), Some((i, FEATURE_OBJECT, _)) if i == id);
        let items = props(oi);
        let old = items.iter().find(|(k, v)| *k == "STYLE" && mine(v)).map(|(_, v)| v.to_vec());
        let mut rest: Vec<(&str, Vec<u8>)> = items.iter().filter(|(k, v)| !(*k == "STYLE" && mine(v))).map(|&(k, v)| (k, v.to_vec())).collect();
        if let Some(st) = st {
            let mut r = old.unwrap_or_else(|| template.clone());
            r[0..4].copy_from_slice(&id.to_le_bytes());
            r[4..8].copy_from_slice(&FEATURE_OBJECT.to_le_bytes());
            r[8] = st.line.is_some() as u8;
            if let Some(c) = st.line { r[9..12].copy_from_slice(&c); }
            r[12..20].copy_from_slice(&st.weight.to_le_bytes());
            r[20] = st.fill.is_some() as u8;
            if let Some(c) = st.fill { r[21..24].copy_from_slice(&c); }
            rest.insert(0, ("STYLE", r)); // SeqBuilder puts the newest object style first
        }
        block(oi, rest)
    };
    // The type style's first field names the object last styled with it (0xFFFFFFFF = none).
    // ponytail: resetting it on removal is a guess; SeqBuilder's own removal is not sampled yet.
    let mark = if st.is_some() { id } else { u32::MAX };
    let typeinfo = |ti: &[u8]| -> Vec<u8> {
        block(ti, props(ti).into_iter().map(|(k, td)| {
            if k != "TDATA" || !is_kind(td) { return (k, td.to_vec()); }
            (k, block(td, props(td).into_iter().map(|(k2, v)| {
                let mut v = v.to_vec();
                if k2 == "STYLE" && (st.is_some() || u32le(&v, 0) == Some(id)) { v[0..4].copy_from_slice(&mark.to_le_bytes()); }
                (k2, v)
            })))
        }))
    };
    block(vs, props(vs).into_iter().map(|(k, v)| match k {
        "OBJINFO" => (k, objinfo(v)),
        "TYPEINFO" => (k, typeinfo(v)),
        _ => (k, v.to_vec()),
    }))
}

impl Sbd {
    /// Add a feature as Features ▸ New Feature does: next unused id, a UVFT after the last one, one more
    /// FXIF, and a shown entry in every feature state. Returns the id.
    pub fn add_feature(&mut self, mut f: SbdFeat) -> u32 {
        let used = self.shown.keys().copied().chain(self.container.ids(b"UVFT"));
        f.id = used.max().map_or(0, |m| m + 1);
        self.container.add_after_last(b"UVFT", f.id, f.encode());
        let fx = self.container.ids(b"FXIF");
        let template = fx.last().and_then(|&i| self.container.chunk_id_mut(b"FXIF", Some(i)).map(|d| d.clone()));
        if let (Some(&last), Some(data)) = (fx.last(), template) { self.container.add_after_last(b"FXIF", last + 1, data); }
        if let Some(tl) = self.container.chunk_mut(b"UVTL") { edit_feature_state(tl, |rows| rows.push((f.id, true))); }
        self.shown.insert(f.id, true);
        let id = f.id;
        self.feats.push(f);
        id
    }

    /// Delete a feature as SeqBuilder does: its UVFT and the last FXIF go, ids are not renumbered.
    pub fn remove_feature(&mut self, id: u32) {
        self.container.remove(b"UVFT", id);
        if let Some(&last) = self.container.ids(b"FXIF").last() { self.container.remove(b"FXIF", last); }
        self.feats.retain(|f| f.id != id);
    }

    /// The Show checkbox of the feature table.
    pub fn set_feature_shown(&mut self, id: u32, shown: bool) {
        if let Some(tl) = self.container.chunk_mut(b"UVTL") {
            edit_feature_state(tl, |rows| rows.iter_mut().filter(|r| r.0 == id).for_each(|r| r.1 = shown));
        }
        self.shown.insert(id, shown);
    }

    /// Replace the sequence (features are moved by the caller and saved with `update_feature`).
    pub fn set_sequence(&mut self, seq: Vec<u8>) {
        if let Some(sq) = self.container.chunk_mut(b"UVSQ") { *sq = seq.clone(); }
        if let Some(si) = self.container.chunk_mut(b"UVSI") { si[16..20].copy_from_slice(&(seq.len() as u32).to_be_bytes()); }
        self.seq = seq;
    }

    /// Replace an existing feature (matched by id) with edited name, type, range or qualifiers.
    pub fn update_feature(&mut self, f: SbdFeat) -> Result<(), String> {
        let data = self.container.chunk_id_mut(b"UVFT", Some(f.id)).ok_or(format!("no feature {}", f.id))?;
        *data = f.encode();
        if let Some(old) = self.feats.iter_mut().find(|o| o.id == f.id) { *old = f; }
        Ok(())
    }

    /// The document as SeqBuilder Pro would read it; unchanged chunks are written back byte for byte.
    pub fn to_bytes(&self) -> Vec<u8> {
        self.container.write()
    }

    /// Set (or with `None` remove) the custom style of feature `id`, as Style ▸ Fill/Line does in SeqBuilder:
    /// every VSSTYLE copy gets the object style first in OBJINFO, and its type style records the object id.
    pub fn set_feature_style(&mut self, id: u32, st: Option<Style>) {
        let kind = self.feats.iter().find(|f| f.id == id).map(|f| f.kind.clone()).unwrap_or_default();
        let Some(tl) = self.container.chunk_mut(b"UVTL") else { return };
        edit_all(tl, "VSSTYLE", |vs| restyle(vs, id, &kind, st));
        match st {
            Some(st) => { self.obj_styles.insert(id, st); }
            None => { self.obj_styles.remove(&id); }
        }
    }

    /// A GenBank header field (DEFINITION, ACCESSION, COMMENT...) from UVTX, continuation lines included.
    pub fn header_field(&self, key: &str, sep: &str) -> String {
        let mut lines = self.header.lines().skip_while(|l| l.split_whitespace().next() != Some(key));
        let Some(first) = lines.next() else { return String::new() };
        std::iter::once(first[key.len()..].trim())
            .chain(lines.take_while(|l| l.starts_with("            ")).map(str::trim))
            .collect::<Vec<_>>().join(sep)
    }
}

/// The name SeqBuilder shows for a feature: /gene, else /product, else the first qualifier's value,
/// without parenthesised parts, cut at a word to at most 40 characters.
pub fn feature_name(quals: &[String]) -> String {
    let value = |l: &String| l.split_once('=').map_or(String::new(), |(_, v)| v.trim_matches('"').to_string());
    let find = |k: &str| quals.iter().find(|l| l.split_once('=').is_some_and(|(key, _)| key == k)).map(value);
    let raw = find("gene").or_else(|| find("product")).or_else(|| quals.first().map(value)).unwrap_or_default();
    let mut name = String::new();
    let mut depth = 0;
    for c in raw.chars() {
        match c {
            '(' => depth += 1,
            ')' if depth > 0 => depth -= 1,
            _ if depth == 0 => name.push(c),
            _ => {}
        }
    }
    let mut name = name.split_whitespace().collect::<Vec<_>>().join(" ");
    while name.chars().count() > 40 {
        name.truncate(name.rfind(' ').unwrap_or(40.min(name.len())));
    }
    name
}

/// A new document made from `template` (any SeqBuilder document): its layouts and settings are kept,
/// everything about the sequence and its features is replaced.
pub fn new_document(template: &[u8], name: &str, header: &str, seq: Vec<u8>, circular: bool, feats: Vec<SbdFeat>) -> Result<Sbd, String> {
    let mut s = read(template)?;
    s.set_sequence(seq);
    if let Some(si) = s.container.chunk_mut(b"UVSI") {
        if let Some(at) = si.windows(8).position(|w| w == b"SeqShape") {
            si[at + 12..at + 16].copy_from_slice(&(1 + circular as u32).to_le_bytes());
        }
    }
    s.circular = circular;
    let latin1 = |t: &str| t.chars().map(|c| u8::try_from(c as u32).unwrap_or(b'?')).collect::<Vec<u8>>();
    // A sequence that never was GenBank stores "^^"
    let text = if header.is_empty() { "^^\n".to_string() } else { header.to_string() };
    if let Some(tx) = s.container.chunk_mut(b"UVTX") { *tx = latin1(&text.replace('\n', "\r\n")); }
    s.header = header.to_string();
    if let Some(cb) = s.container.chunk_mut(b"UVCB") {
        let mut v = (name.len() as u32).to_le_bytes().to_vec();
        v.extend(latin1(name));
        edit_all(cb, "SeqName", |_| v.clone());
    }
    let feats: Vec<SbdFeat> = feats.into_iter().enumerate().map(|(i, f)| SbdFeat { id: i as u32, ..f }).collect();
    let fxif = s.container.chunks().find(|c| &c.0 == b"FXIF").map(|c| c.2.to_vec()).ok_or("template has no FXIF")?;
    s.container.replace_all(b"UVFT", feats.iter().map(|f| (f.id, f.encode())).collect()).ok_or("template has no features")?;
    s.container.replace_all(b"FXIF", (1..=feats.len() as u32).map(|i| (i, fxif.clone())).collect());
    // SeqBuilder lists source features hidden
    s.shown = feats.iter().map(|f| (f.id, f.kind != "source")).collect();
    let rows: Vec<(u32, bool)> = feats.iter().map(|f| (f.id, s.shown[&f.id])).collect();
    if let Some(tl) = s.container.chunk_mut(b"UVTL") {
        edit_feature_state(tl, |r| *r = rows.clone());
        // the template's per-feature styles belong to features that are gone
        edit_all(tl, "VSSTYLE", |vs| block(vs, props(vs).into_iter().map(|(k, v)| {
            if k != "OBJINFO" { return (k, v.to_vec()); }
            (k, block(v, props(v).into_iter().filter(|(k, st)| !(*k == "STYLE" && matches!(style(st), Some((_, FEATURE_OBJECT, _))))).map(|(k, st)| (k, st.to_vec()))))
        })));
    }
    s.obj_styles.clear();
    s.feats = feats;
    Ok(s)
}

/// The UVSI chunk of a document (for comparing saves).
#[cfg(test)]
pub fn uvsi(d: &[u8]) -> Vec<u8> {
    Container::parse(d).and_then(|c| c.chunks().find(|c| &c.0 == b"UVSI").map(|c| c.2.to_vec())).unwrap_or_default()
}

pub fn read(d: &[u8]) -> Result<Sbd, String> {
    let container = Container::parse(d).ok_or("not a SeqBuilder document (bad chunk directory)")?;
    let chunks: Vec<_> = container.chunks().collect();
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
        let quals = text(r.get(416..)?).split('\r').filter_map(|l| l.strip_prefix('/').map(str::to_owned)).collect();
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
                    let (Some(ti), Some(raw)) = (child(td, "TINFO"), child(td, "STYLE")) else { continue };
                    let (Some(st), n) = (style(raw), u32le(ti, 0).unwrap_or(0) as usize) else { continue };
                    if let Some(name) = ti.get(4..4 + n) { type_styles.insert(text(name), st.2); }
                },
                _ => {}
            }
        }
    }
    let shown = get(*b"UVTL").next()
        .and_then(|c| ["ALLDSTAT", "CMNSTATE", "VSFEAT", "STATE"].iter().try_fold(c.2, |b, k| child(b, k)))
        .map(|st| feature_state(st).into_iter().collect()).unwrap_or_default();
    drop(chunks);
    Ok(Sbd { container, shown, seq, header, circular, feats, type_styles, obj_styles })
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
        assert_eq!(p1.qual("note").as_deref(), Some("promoter P1 (6)"));
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

    #[test]
    fn feature_records_encode_back_exactly() {
        for f in ["pBR322.sbd", "ALDH2.sbd"] {
            let d = std::fs::read(format!("testdata/{f}")).unwrap();
            let c = Container::parse(&d).unwrap();
            let s = read(&d).unwrap();
            for ((_, id, raw), feat) in c.chunks().filter(|c| &c.0 == b"UVFT").zip(&s.feats) {
                assert!(feat.encode() == raw, "{f} feature {id}");
            }
        }
    }

    fn items_of(s: &Sbd, key: &str) -> Vec<Vec<u8>> {
        let tl = s.container.chunks().find(|c| &c.0 == b"UVTL").unwrap().2.to_vec();
        find_items(&tl, key).into_iter().map(|at| tl[at + 12..at + 12 + u32le(&tl, at + 8).unwrap() as usize].to_vec()).collect()
    }

    /// Adding TEST 100..200 must give the chunks and feature states SeqBuilder Pro wrote in pBR322_add.sbd.
    #[test]
    fn add_feature_matches_seqbuilder() {
        let mut d = read(&std::fs::read("testdata/pBR322.sbd").unwrap()).unwrap();
        let mut want = read(&std::fs::read("testdata/pBR322_add.sbd").unwrap()).unwrap();
        assert_eq!(d.shown.get(&0), Some(&false)); // source is listed but hidden
        let id = d.add_feature(SbdFeat {
            id: 0, name: "TEST".into(), kind: "misc_feature".into(), start: 100, end: 200, complement: false,
            quals: vec!["dnas_title=\"TEST\"".into(), "note=\"TEST_note\"".into()],
        });
        assert_eq!(id, 50);
        let dir = |s: &Sbd| s.container.chunks().map(|c| (c.0, c.1)).collect::<Vec<_>>();
        assert_eq!(dir(&d), dir(&want));
        for typ in [b"UVFT", b"FXIF"] {
            for i in d.container.ids(typ) {
                assert!(d.container.chunk_id_mut(typ, Some(i)).cloned() == want.container.chunk_id_mut(typ, Some(i)).cloned());
            }
        }
        // SeqBuilder also re-saves VSFEAT.TSTYLE with its type entries in another order; only STATE carries the edit.
        let states = |s: &Sbd| items_of(s, "VSFEAT").iter().map(|v| child(v, "STATE").unwrap().to_vec()).collect::<Vec<_>>();
        assert_eq!(states(&d).len(), 9);
        assert!(states(&d) == states(&want));
        let again = read(&d.to_bytes()).unwrap();
        assert_eq!((again.feats.len(), again.feats[50].name.as_str(), again.shown[&50]), (51, "TEST", true));
    }

    /// Deleting echinomycin 24..27 (UVFT 2) must leave the chunks SeqBuilder Pro left in pBR322_del.sbd.
    #[test]
    fn remove_feature_matches_seqbuilder() {
        let mut d = read(&std::fs::read("testdata/pBR322.sbd").unwrap()).unwrap();
        let want = read(&std::fs::read("testdata/pBR322_del.sbd").unwrap()).unwrap();
        d.remove_feature(2);
        let dir = |s: &Sbd| s.container.chunks().map(|c| (c.0, c.1, c.2.to_vec())).filter(|c| &c.0 != b"UVTL" && &c.0 != b"UVCB").collect::<Vec<_>>();
        assert!(dir(&d) == dir(&want));
        assert_eq!(read(&d.to_bytes()).unwrap().feats.len(), 49);
    }

    #[test]
    fn writes_back_unchanged_files_byte_for_byte() {
        for f in ["pBR322.sbd", "pBR322_blaRED.sbd", "ALDH2.sbd"] {
            let d = std::fs::read(format!("testdata/{f}")).unwrap();
            assert!(read(&d).unwrap().to_bytes() == d, "{f}");
        }
    }

    /// Recolouring bla red must give the same style blocks as SeqBuilder Pro saved in pBR322_blaRED.sbd.
    #[test]
    fn custom_style_matches_seqbuilder() {
        let blocks = |s: &Sbd| {
            let tl = s.container.chunks().find(|c| &c.0 == b"UVTL").unwrap().2.to_vec();
            find_items(&tl, "VSSTYLE").into_iter().map(|at| tl[at + 12..at + 12 + u32le(&tl, at + 8).unwrap() as usize].to_vec()).collect::<Vec<_>>()
        };
        let mut d = read(&std::fs::read("testdata/pBR322.sbd").unwrap()).unwrap();
        let red = read(&std::fs::read("testdata/pBR322_blaRED.sbd").unwrap()).unwrap();
        assert_eq!(blocks(&d).len(), 8);
        d.set_feature_style(38, Some(red.obj_styles[&38]));
        assert!(blocks(&d) == blocks(&red));

        // Putting the type's colours back keeps the object style, exactly as SeqBuilder does.
        let reset = read(&std::fs::read("testdata/pBR322_blaRESET.sbd").unwrap()).unwrap();
        let mut back = read(&std::fs::read("testdata/pBR322_blaRED.sbd").unwrap()).unwrap();
        back.set_feature_style(38, Some(reset.type_styles["CDS"]));
        assert!(blocks(&back) == blocks(&reset));

        let again = read(&d.to_bytes()).unwrap();
        assert_eq!(again.obj_styles[&38], red.obj_styles[&38]);
        assert_eq!((again.seq.len(), again.feats.len()), (4361, 50));
        d.set_feature_style(38, None);
        assert!(read(&d.to_bytes()).unwrap().to_bytes() == std::fs::read("testdata/pBR322.sbd").unwrap());
    }
}
