// SeqBuilder-like DNA viewer laid out like SeqBuilder Pro: tool strip, map page with view tabs,
// feature table with tabs, Style/Enzymes and Settings/Details side panels.
mod enzymes;
mod genbank;
mod sbd;

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, Vec2};
use enzymes::{Enzyme, Selector};
use std::collections::HashSet;
use std::f32::consts::{PI, TAU};
use std::path::{Path, PathBuf};

const DNASTAR_DIR: &str = "/Applications/DNASTAR/Lasergene19";
// New .sbd documents start from this SeqBuilder Pro document (New DNA "acgt" with one feature).
const SBD_TEMPLATE: &[u8] = include_bytes!("../testdata/template.sbd");

#[derive(Clone)]
struct Feat {
    kind: String,
    name: String,
    label: String, // name shortened for drawing
    desc: String,
    range: String, // GenBank location text, shown in the table
    start: usize,  // 0-based, inclusive
    end: usize,    // exclusive; end < start means it wraps the origin
    strand: i8,    // 1 forward, -1 reverse
    lane: usize,
    show: bool,
    style: Option<sbd::Style>,         // this feature's own style, from the document or the Style tab
    base: Option<sbd::Style>,          // the document's style for this feature type
    sbd_id: Option<u32>,               // UVFT id when read from a .sbd
    edited: bool,                      // style changed since the last save
    shape: Option<sbd::Shape>,         // this feature's own shape, linkage and label
    base_shape: Option<sbd::Shape>,    // the document's shape for this feature type
    shape_edited: bool,
    quals: Vec<String>,                // "key=value" qualifier lines, quoting included
    dirty: bool,                       // name, type, range or qualifiers changed since the last save
}

/// A feature as stored in .sbd (1-based, start > end on the complement strand).
fn sbd_feat(f: &Feat, id: u32) -> sbd::SbdFeat {
    let (lo, hi, comp) = (f.start + 1, f.end, f.strand < 0);
    let (start, end) = if comp { (hi, lo) } else { (lo, hi) };
    sbd::SbdFeat { id, name: f.name.clone(), kind: f.kind.clone(), start, end, complement: comp, quals: f.quals.clone() }
}

/// GenBank location: as read when it holds what the coordinates cannot (join, order, partial ends)
/// and the feature was not edited, otherwise from the coordinates.
fn gb_location(f: &Feat) -> String {
    let complex = ["join(", "order(", "<", ">"].iter().any(|m| f.range.contains(m));
    if f.start >= f.end || (complex && !f.dirty) { return f.range.clone(); }
    let span = format!("{}..{}", f.start + 1, f.end); // SeqBuilder writes a single base as 426..426
    if f.strand < 0 { format!("complement({span})") } else { span }
}

/// The document as a GenBank flat file, as SeqBuilder exports it: features by start position and the
/// document name as /dnas_title on the first source feature.
fn export_genbank(d: &Doc) -> String {
    let header = if d.header.starts_with("^^") { "" } else { &d.header };
    let title = d.path.file_stem().map_or(d.name.clone(), |s| s.to_string_lossy().into_owned());
    let first_source = d.feats.iter().position(|f| f.kind == "source");
    let quals: Vec<Vec<String>> = d.feats.iter().enumerate().map(|(i, f)| {
        let mut q = f.quals.clone();
        if Some(i) == first_source && !q.iter().any(|l| l.starts_with("dnas_title=")) { q.insert(0, format!("dnas_title=\"{title}\"")); }
        q
    }).collect();
    let mut order: Vec<usize> = (0..d.feats.len()).collect();
    order.sort_by_key(|&i| d.feats[i].start); // stable: equal starts keep document order
    let feats: Vec<genbank::Feature> = order.iter()
        .map(|&i| genbank::Feature { kind: &d.feats[i].kind, location: gb_location(&d.feats[i]), quals: &quals[i] }).collect();
    genbank::write(header, &d.name, &d.seq, d.circular, &feats)
}

/// The Range column for a feature that does not wrap the origin.
fn range_text(f: &Feat) -> String {
    let (lo, hi) = (f.start + 1, f.end);
    if f.strand < 0 { format!("complement({lo}..{hi})") } else { format!("{lo}..{hi}") }
}

/// Set (or with an empty value remove) a qualifier line.
fn set_qual(quals: &mut Vec<String>, key: &str, value: &str) {
    let line = format!("{key}=\"{value}\"");
    match quals.iter().position(|l| l.strip_prefix(key).is_some_and(|r| r.starts_with('='))) {
        Some(i) if value.is_empty() => { quals.remove(i); }
        Some(i) => quals[i] = line,
        None if !value.is_empty() => quals.push(line),
        None => {}
    }
}

struct Doc {
    name: String,
    definition: String,
    accession: String,
    comment: String,
    seq: Vec<u8>,
    circular: bool,
    feats: Vec<Feat>,
    lanes: usize,
    hits: Vec<Vec<usize>>, // per enzyme: 0-based site positions
    path: PathBuf,
    sbd: Option<sbd::Sbd>, // kept to write the document back
    deleted: Vec<u32>,     // .sbd feature ids deleted since the last save
    seq_dirty: bool,       // sequence edited since the last save
    header: String,        // GenBank header (LOCUS..COMMENT) for a new .sbd; empty for FASTA
    topology_dirty: bool,  // switched linear/circular since the last save
}

struct Site {
    enz: usize, // index into App::enz
    pos: usize, // 0-based start of the recognition sequence
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G', b'U' => b'A',
        b'a' => b't', b't' => b'a', b'g' => b'c', b'c' => b'g', b'u' => b'a',
        b'R' => b'Y', b'Y' => b'R', b'K' => b'M', b'M' => b'K',
        b'B' => b'V', b'V' => b'B', b'D' => b'H', b'H' => b'D',
        b'N' | b'n' | b'S' | b'W' => b,
        _ => b'N',
    }
}

fn gc_percent(s: &[u8]) -> f32 {
    if s.is_empty() { return 0.0; }
    let gc = s.iter().filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c')).count();
    100.0 * gc as f32 / s.len() as f32
}

// IUPAC code as a set of bases: A=1, C=2, G=4, T=8.
fn mask(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 1, b'C' => 2, b'G' => 4, b'T' | b'U' => 8,
        b'M' => 3, b'R' => 5, b'W' => 9, b'S' => 6, b'Y' => 10, b'K' => 12,
        b'V' => 7, b'H' => 11, b'D' => 13, b'B' => 14,
        _ => 15,
    }
}

/// Sequence base `b` fits site letter `p`. An ambiguous base only fits when it surely does,
/// unless `potential` ("Show Potential Cuts at Ambiguous Bases") accepts any overlap.
fn iupac(p: u8, b: u8, potential: bool) -> bool {
    let (p, b) = (mask(p), mask(b));
    if potential { p & b != 0 } else { b & !p == 0 }
}

// ponytail: naive scan, O(len × enzymes × site); fine for plasmids, needs an index for whole genomes.
fn find_sites(seq: &[u8], circular: bool, site: &str, potential: bool) -> Vec<usize> {
    let fwd = site.as_bytes();
    let rev: Vec<u8> = fwd.iter().rev().map(|&b| complement(b)).collect();
    let n = seq.len();
    if n < fwd.len() { return vec![]; }
    let last = if circular { n } else { n - fwd.len() + 1 };
    (0..last).filter(|&i| {
        let hit = |p: &[u8]| p.iter().enumerate().all(|(k, &c)| iupac(c, seq[(i + k) % n], potential));
        hit(fwd) || (rev != fwd && hit(&rev))
    }).collect()
}

fn covers(f: &Feat, a: usize, b: usize) -> bool {
    if f.start < f.end { f.start < b && a < f.end } else { f.start < b || a < f.end }
}

// Greedy lane packing so overlapping visible features don't draw on top of each other.
fn assign_lanes(feats: &mut [Feat], len: usize) -> usize {
    let mut lanes: Vec<Vec<(usize, usize)>> = vec![];
    for f in feats.iter_mut().filter(|f| f.show) {
        let spans = if f.start < f.end { vec![(f.start, f.end)] } else { vec![(f.start, len), (0, f.end)] };
        let free = |lane: &Vec<(usize, usize)>| spans.iter().all(|&(a, b)| lane.iter().all(|&(c, d)| b <= c || d <= a));
        let i = lanes.iter().position(free).unwrap_or_else(|| { lanes.push(vec![]); lanes.len() - 1 });
        lanes[i].extend(&spans);
        f.lane = i;
    }
    lanes.len()
}

fn location_strand(l: &gb_io::seq::Location) -> i8 {
    use gb_io::seq::Location::*;
    match l {
        Complement(_) => -1,
        Join(v) | Order(v) if !v.is_empty() && v.iter().all(|x| matches!(x, Complement(_))) => -1,
        _ => 1,
    }
}

fn short(label: &str) -> String {
    if label.chars().count() > 28 { label.chars().take(27).collect::<String>() + "…" } else { label.to_string() }
}

fn load(path: &Path) -> Result<Doc, String> {
    let name = path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_default();
    let text = std::fs::read(path).map_err(|e| e.to_string())?;
    let mut doc = if text.first() == Some(&b'>') {
        let seq: Vec<u8> = text.split(|&b| b == b'\n').filter(|l| !l.starts_with(b">"))
            .flatten().copied().filter(u8::is_ascii_alphabetic).collect();
        Doc { name, definition: String::new(), accession: String::new(), comment: String::new(), seq, circular: false, feats: vec![], lanes: 0, hits: vec![], path: path.into(), sbd: None, deleted: vec![], seq_dirty: false, header: String::new(), topology_dirty: false }
    } else if path.extension().is_some_and(|e| e.eq_ignore_ascii_case("sbd")) {
        let s = sbd::read(&text)?;
        let feats = s.feats.iter().map(|f| {
            let (lo, hi) = (f.start.min(f.end), f.start.max(f.end));
            let q = |k: &str| f.qual(k);
            Feat {
                kind: f.kind.clone(), name: f.name.clone(), label: short(&f.name), desc: q("note").or_else(|| q("product")).unwrap_or_default(),
                range: if f.complement { format!("complement({lo}..{hi})") } else { format!("{lo}..{hi}") },
                start: lo.saturating_sub(1), end: hi, strand: if f.complement { -1 } else { 1 }, lane: 0,
                show: s.shown.get(&f.id).copied().unwrap_or(f.kind != "source"), style: s.obj_styles.get(&f.id).copied(),
                shape: s.obj_shapes.get(&f.id).copied(), base_shape: s.type_shapes.get(&f.kind).copied(), shape_edited: false,
                base: s.type_styles.get(&f.kind).copied(), sbd_id: Some(f.id), edited: false,
                quals: f.quals.clone(), dirty: false,
            }
        }).collect();
        let locus = s.header_field("LOCUS", " ");
        Doc {
            name: locus.split_whitespace().next().map_or(name, str::to_owned),
            definition: s.header_field("DEFINITION", " "), accession: s.header_field("ACCESSION", " "),
            comment: s.header_field("COMMENT", "\n"),
            seq: s.seq.clone(), circular: s.circular, feats, lanes: 0, hits: vec![], path: path.into(), header: s.header.clone(), sbd: Some(s), deleted: vec![], seq_dirty: false, topology_dirty: false,
        }
    } else {
        let rec = gb_io::reader::parse_slice(&text).map_err(|e| e.to_string())?
            .into_iter().next().ok_or("no sequence in file")?;
        let len = rec.seq.len();
        let feats = rec.features.iter().filter_map(|f| {
            let (s, e) = f.location.find_bounds().ok()?;
            let q = |k: &str| f.qualifier_values(k).next().map(str::to_owned);
            // INSDC qualifiers whose values are written without quotes
            const BARE: [&str; 13] = ["anticodon", "citation", "codon_start", "compare", "direction", "estimated_length",
                "mod_base", "number", "rpt_type", "rpt_unit_range", "tag_peptide", "transl_except", "transl_table"];
            let quals: Vec<String> = f.qualifiers.iter().map(|(k, v)| match v {
                Some(v) if BARE.contains(&k.as_ref()) => format!("{k}={v}"),
                // continuation lines join with a space, except inside a protein translation
                Some(v) => format!("{k}=\"{}\"", v.replace('\n', if k == "translation" { "" } else { " " })),
                None => k.to_string(),
            }).collect();
            let name = sbd::feature_name(&quals);
            // find_bounds unwraps origin-spanning joins past len; fold them back.
            let (start, end) = (s.rem_euclid(len as i64) as usize, e as usize);
            let end = if end > len { end - len } else { end };
            Some(Feat {
                kind: f.kind.to_string(), label: short(&name), name, desc: q("note").or_else(|| q("product")).unwrap_or_default(),
                range: f.location.to_gb_format(), start, end, strand: location_strand(&f.location), lane: 0,
                show: f.kind != "source", // SeqBuilder lists source but leaves it hidden
                style: None, base: TYPE_LOOKS.0.get(f.kind.as_ref()).copied(), sbd_id: None, edited: false, quals, dirty: false,
                shape: None, base_shape: TYPE_LOOKS.1.get(f.kind.as_ref()).copied(), shape_edited: false,
            })
        }).collect();
        let circular = rec.is_circular();
        Doc {
            name: rec.name.clone().unwrap_or(name), definition: rec.definition.clone().unwrap_or_default(),
            accession: rec.accession.clone().unwrap_or_default(), comment: rec.comments.join("\n"),
            seq: rec.seq, circular, feats, lanes: 0, hits: vec![], path: path.into(), sbd: None, deleted: vec![], seq_dirty: false,
            header: text.iter().map(|&b| b as char).collect::<String>().split("\nFEATURES").next().unwrap_or_default().to_string() + "\n",
            topology_dirty: false,
        }
    };
    doc.lanes = assign_lanes(&mut doc.feats, doc.seq.len());
    Ok(doc)
}

/// Style for documents that carry none (GenBank, FASTA).
fn kind_style(kind: &str) -> sbd::Style {
    let (fill, line) = match kind {
        "CDS" | "gene" | "mat_peptide" => ([215, 228, 255], [70, 110, 220]),
        "promoter" => ([170, 225, 160], [40, 130, 40]),
        "rep_origin" => ([240, 235, 120], [150, 140, 20]),
        "terminator" => ([240, 170, 170], [170, 40, 40]),
        "primer_bind" => ([170, 210, 245], [30, 110, 190]),
        "misc_binding" | "protein_bind" => ([215, 190, 245], [130, 60, 200]),
        "misc_feature" => ([230, 200, 160], [130, 75, 20]),
        _ => ([215, 215, 220], [110, 110, 125]),
    };
    sbd::Style { line: Some(line), weight: 1.0, fill: Some(fill), shadow: None }
}

/// Type styles and shapes of a SeqBuilder document, used for documents that carry none (GenBank, FASTA),
/// so they are drawn as SeqBuilder draws them after import.
static TYPE_LOOKS: std::sync::LazyLock<(std::collections::HashMap<String, sbd::Style>, std::collections::HashMap<String, sbd::Shape>)> =
    std::sync::LazyLock::new(|| sbd::read(SBD_TEMPLATE).map(|t| (t.type_styles, t.type_shapes)).unwrap_or_default());

const SHAPES: [&str; 6] = ["Box", "Filled Arrow", "Filled Double Arrow", "Triangle Glyph", "Stroked Arrow", "Stroked Goalpost"];
const LINKAGES: [&str; 3] = ["Solid", "Dashed", "None"];
const LABEL_POSITIONS: [&str; 3] = ["Above", "Inside", "Below"];

fn shape(f: &Feat) -> sbd::Shape {
    f.shape.or(f.base_shape).unwrap_or(sbd::Shape { shape: 2, linkage: 1, label: 2, attached: true })
}

fn look(f: &Feat) -> sbd::Style {
    f.style.or(f.base).unwrap_or_else(|| kind_style(&f.kind))
}

fn rgb(c: [u8; 3]) -> Color32 {
    Color32::from_rgb(c[0], c[1], c[2])
}

/// (fill, outline) as drawn; a style without fill or line draws none.
fn colors(f: &Feat) -> (Color32, Color32) {
    let s = look(f);
    (s.fill.map_or(Color32::TRANSPARENT, rgb), s.line.map_or(Color32::TRANSPARENT, rgb))
}

/// Base composition rows of the Details panel: single bases, A+T, G+C, ambiguous, CpG and all
/// dinucleotides, as "percent [count]".
fn base_stats(s: &[u8]) -> Vec<(String, String)> {
    let up: Vec<u8> = s.iter().map(|b| b.to_ascii_uppercase()).collect();
    let n = up.len();
    let count = |c: u8| up.iter().filter(|&&b| b == c).count();
    let (a, c, g, t) = (count(b'A'), count(b'C'), count(b'G'), count(b'T') + count(b'U'));
    let pct = |k: usize, of: usize| format!("{:.2} [{k}]", if of == 0 { 0.0 } else { 100.0 * k as f64 / of as f64 });
    let pairs = n.saturating_sub(1);
    let di = |x: u8, y: u8| up.windows(2).filter(|w| w[0] == x && w[1] == y).count();
    let mut rows = vec![
        ("% A".to_string(), pct(a, n)), ("% A+T".into(), pct(a + t, n)), ("% Ambiguous".into(), pct(n - a - c - g - t, n)),
        ("% C".into(), pct(c, n)), ("% G".into(), pct(g, n)), ("% G+C".into(), pct(g + c, n)), ("% T".into(), pct(t, n)),
        ("% diCG".into(), pct(di(b'C', b'G'), pairs)), ("Total number of bases".into(), n.to_string()),
    ];
    for x in *b"ACGT" {
        for y in *b"ACGT" { rows.push((format!("% {}{}", x as char, y as char), pct(di(x, y), pairs))); }
    }
    rows
}

fn statistics(ui: &mut egui::Ui, s: &[u8]) {
    ui.label(egui::RichText::new("Statistics:").italics());
    egui::Grid::new("statistics").striped(true).num_columns(2).show(ui, |ui| {
        for (k, v) in base_stats(s) { ui.label(k); ui.label(v); ui.end_row(); }
    });
}

/// A primer's own sequence: its /primer qualifier, else the bases it covers read on its strand.
fn primer_sequence(f: &Feat, seq: &[u8]) -> Vec<u8> {
    if let Some(p) = f.quals.iter().find_map(|q| q.strip_prefix("primer=")) { return p.trim_matches('"').as_bytes().to_vec(); }
    if f.start >= f.end { return vec![]; }
    let s = &seq[f.start..f.end];
    if f.strand < 0 { s.iter().rev().map(|&b| complement(b)).collect() } else { s.to_vec() }
}

// ponytail: nearest-neighbour Tm (SantaLucia 1998, 50 nM primer, 0.1 M Na+) reproduces SeqBuilder's 77.7 °C for
// ALDH2's 38-mers, but one value fits several models; check against more SeqBuilder primers.
fn primer_tm(p: &[u8]) -> f64 {
    let up: Vec<u8> = p.iter().map(|b| b.to_ascii_uppercase()).collect();
    if up.len() < 2 || up.iter().any(|b| !b"ACGT".contains(b)) { return f64::NAN; }
    let nn = |d: &[u8]| -> (f64, f64) {
        match d {
            b"AA" | b"TT" => (-7.9, -22.2), b"AT" => (-7.2, -20.4), b"TA" => (-7.2, -21.3),
            b"CA" | b"TG" => (-8.5, -22.7), b"GT" | b"AC" => (-8.4, -22.4), b"CT" | b"AG" => (-7.8, -21.0),
            b"GA" | b"TC" => (-8.2, -22.2), b"CG" => (-10.6, -27.2), b"GC" => (-9.8, -24.4), _ => (-8.0, -19.9), // GG, CC
        }
    };
    let (mut dh, mut ds) = up.windows(2).map(nn).fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
    for end in [up[0], up[up.len() - 1]] {
        let (h, s) = if matches!(end, b'G' | b'C') { (0.1, -2.8) } else { (2.3, 4.1) };
        (dh, ds) = (dh + h, ds + s);
    }
    let (na, ct) = (0.1f64, 50e-9f64);
    ds += 0.368 * up.len() as f64 * na.ln();
    dh * 1000.0 / (ds + 1.987 * (ct / 2.0).ln()) - 273.15
}

/// Standard genetic code, TCAG order.
const CODONS: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

fn translate(codon: &[u8]) -> char {
    let n = |b: u8| match b.to_ascii_uppercase() { b'T' | b'U' => Some(0), b'C' => Some(1), b'A' => Some(2), b'G' => Some(3), _ => None };
    match (n(codon[0]), n(codon[1]), n(codon[2])) {
        (Some(a), Some(b), Some(c)) => CODONS[a * 16 + b * 4 + c] as char,
        _ => 'X',
    }
}

/// The amino acid of reading frame `k` (0-2 forward from base 1, 3-5 reverse from the last base)
/// whose codon has base `i` in the middle, if there is one.
fn frame_residue(seq: &[u8], i: usize, k: usize) -> Option<char> {
    if i == 0 || i + 1 >= seq.len() { return None; }
    let codon = &seq[i - 1..=i + 1];
    if k < 3 {
        ((i - 1) % 3 == k).then(|| translate(codon))
    } else {
        let j = seq.len() - 2 - i; // codon start on the reverse strand
        (j % 3 == k - 3).then(|| translate(&codon.iter().rev().map(|&b| complement(b)).collect::<Vec<u8>>()))
    }
}

/// One feature span on a straight track (Linear view), in its Style tab shape.
fn draw_track_feature(p: &egui::Painter, r: Rect, strand: i8, shape: u32, fill: Color32, line: Color32, st: sbd::Style) {
    let fill = if shape == 5 { Color32::TRANSPARENT } else { fill };
    let outline = Stroke::new(st.weight as f32, line);
    let head = 7.0f32.min(r.width() / 2.0);
    let (lt, rt) = match shape {
        2 | 5 => (strand < 0, strand >= 0),
        3 => (true, true),
        _ => (false, false),
    };
    let body = Rect::from_x_y_ranges(r.left() + if lt { head } else { 0.0 }..=r.right() - if rt { head } else { 0.0 }, r.y_range());
    let tri = |tip: f32, back: f32| vec![Pos2::new(tip, r.center().y), Pos2::new(back, r.top() - 2.0), Pos2::new(back, r.bottom() + 2.0)];
    let mut parts: Vec<Vec<Pos2>> = vec![];
    if lt { parts.push(tri(r.left(), body.left())); }
    if rt { parts.push(tri(r.right(), body.right())); }
    match shape {
        4 => {
            let c = r.center().x;
            let t = if strand < 0 { tri(c - head / 2.0, c + head / 2.0) } else { tri(c + head / 2.0, c - head / 2.0) };
            if let Some(sh) = st.shadow { p.add(Shape::convex_polygon(t.iter().map(|q| *q + Vec2::splat(2.0)).collect(), rgb(sh), Stroke::NONE)); }
            p.add(Shape::convex_polygon(t, fill, outline));
        }
        6 => {
            let stroke = Stroke::new((st.weight as f32).max(1.5), if line == Color32::TRANSPARENT { fill } else { line });
            p.line_segment([r.left_center(), r.right_center()], stroke);
            p.line_segment([r.left_top(), r.left_bottom()], stroke);
            p.line_segment([r.right_top(), r.right_bottom()], stroke);
        }
        _ => {
            if let Some(sh) = st.shadow {
                p.rect_filled(body.translate(Vec2::splat(2.0)), 1.0, rgb(sh));
                for t in &parts { p.add(Shape::convex_polygon(t.iter().map(|q| *q + Vec2::splat(2.0)).collect(), rgb(sh), Stroke::NONE)); }
            }
            p.rect(body, 1.0, fill, outline, egui::StrokeKind::Inside);
            for t in parts { p.add(Shape::convex_polygon(t, fill, Stroke::new(outline.width.min(2.0), line))); }
        }
    }
}

/// SeqBuilder shows short sequences in bp and long ones in kb ("4361 bp", "114.842 kb").
fn size_label(len: usize) -> String {
    if len >= 10_000 { format!("{:.3} kb", len as f64 / 1000.0) } else { format!("{len} bp") }
}

fn tick_step(len: f32) -> usize {
    [100, 250, 500, 1000, 2000, 5000, 10000].into_iter().find(|s| len / *s as f32 <= 12.0).unwrap_or(50000)
}

/// Tabs that SeqBuilder Pro has but this app does not implement yet: shown greyed out to keep the layout.
fn todo_tab(ui: &mut egui::Ui, label: &str) {
    ui.add_enabled(false, egui::Button::new(label).frame(false)).on_disabled_hover_text("Not implemented yet");
}

#[derive(Default, PartialEq, Clone, Copy)]
enum View { #[default] Circular, Linear, Sequence }

#[derive(Default, PartialEq, Clone, Copy)]
enum Tab { #[default] Features, Primers, Comment, Minimap, SiteSummary }

#[derive(Default, PartialEq, Clone, Copy)]
enum SideTop { Style, #[default] Enzymes }

#[derive(Default, PartialEq, Clone, Copy)]
enum SideBottom { Settings, #[default] Details }

#[derive(Default)]
struct App {
    doc: Option<Doc>,
    error: Option<String>,
    enz: Vec<Enzyme>,
    sels: Vec<Selector>,
    shown: HashSet<usize>,        // enzymes checked in the Enzymes tab
    enz_filter: String,
    potential: bool,              // Show Potential Cuts at Ambiguous Bases
    selmgr: Option<u32>,          // Enzyme Selector Manager open, with this selector picked
    sites: Vec<Site>,             // shown sites, sorted by position
    sel: Option<(usize, usize)>,  // half-open
    sel_feat: Option<usize>,
    anchor: usize,
    grouping: usize, // 0 = none, 3, 10
    complement: bool, // show the complementary strand in the Sequence view
    hide_enzymes: bool,
    translations: bool,
    hide_rulers: bool,
    spacing: bool,
    hide_side: bool,
    hide_table: bool,
    zoom: f32,
    view: View,
    tab: Tab,
    side_top: SideTop,
    side_bottom: SideBottom,
    goto_open: bool,
    edit_open: bool, // Features ▸ Edit Feature… dialog
    goto_text: String,
    scroll_to: Option<usize>,
    caret: Option<usize>,                // insertion point in the Sequence view (before this base)
    undo: Vec<(Vec<u8>, Vec<Feat>)>,
    redo: Vec<(Vec<u8>, Vec<Feat>)>,
    frames: u32,
}

impl App {
    fn open(&mut self, path: &Path) {
        match load(path) {
            Ok(d) => {
                self.view = if d.circular { View::Circular } else { View::Sequence };
                self.doc = Some(d);
                self.error = None;
                self.sel = None;
                self.sel_feat = None;
                self.compute_hits();
                // SeqBuilder opens a sequence showing its default selector's enzymes.
                if let Some(id) = self.sels.iter().find(|s| s.name.contains("Unique Six+")).map(|s| s.id) { self.apply_selector(id); }
            }
            Err(e) => self.error = Some(format!("{}: {e}", path.display())),
        }
    }

    /// Details ▸ Feature selection, as SeqBuilder shows it: type, location, orientation, qualifiers, statistics.
    fn feature_details(&self, ui: &mut egui::Ui) {
        let (Some(d), Some(i)) = (&self.doc, self.sel_feat) else { return };
        let Some(f) = d.feats.get(i) else { return };
        egui::ScrollArea::vertical().show(ui, |ui| {
            ui.heading("Feature selection");
            ui.separator();
            ui.label(egui::RichText::new("Feature Information:").italics());
            egui::Grid::new("feature_info").num_columns(2).show(ui, |ui| {
                let orientation = if f.strand < 0 { "Reverse" } else { "Forward" };
                for (k, v) in [("Feature Type:", f.kind.as_str()), ("Location:", f.range.as_str()), ("Orientation:", orientation)] {
                    ui.label(k); ui.label(v); ui.end_row();
                }
                for q in f.quals.iter().filter(|q| !q.starts_with("dnas_title=")) {
                    let (k, v) = q.split_once('=').unwrap_or((q, ""));
                    ui.label(format!("{k}:"));
                    ui.add(egui::Label::new(v.trim_matches('"')).wrap_mode(egui::TextWrapMode::Wrap));
                    ui.end_row();
                }
            });
            if f.start < f.end {
                let mut bases = d.seq[f.start..f.end].to_vec();
                if f.strand < 0 { bases = bases.iter().rev().map(|&b| complement(b)).collect(); }
                statistics(ui, &bases);
            }
        });
    }

    /// Details for the feature picked in the table, editable (name, type, range, strand, note).
    fn feature_editor(&mut self, ui: &mut egui::Ui) {
        let Some(d) = &mut self.doc else { return };
        let len = d.seq.len();
        let Some(f) = self.sel_feat.and_then(|i| d.feats.get_mut(i)) else { return };
        ui.heading("Feature");
        ui.separator();
        let (mut lo, mut hi, mut comp) = (f.start + 1, f.end, f.strand < 0);
        let mut note = f.quals.iter().find_map(|l| l.strip_prefix("note=")).map(|v| v.trim_matches('"').to_string()).unwrap_or_default();
        let wraps = f.start >= f.end; // ponytail: origin-spanning ranges are shown, not edited
        let mut changed = false;
        egui::Grid::new("feature_editor").num_columns(2).show(ui, |ui| {
            ui.label("Name"); changed |= ui.text_edit_singleline(&mut f.name).changed(); ui.end_row();
            ui.label("Type"); changed |= ui.text_edit_singleline(&mut f.kind).changed(); ui.end_row();
            ui.label("Start"); changed |= ui.add_enabled(!wraps, egui::DragValue::new(&mut lo).range(1..=hi)).changed(); ui.end_row();
            ui.label("End"); changed |= ui.add_enabled(!wraps, egui::DragValue::new(&mut hi).range(lo..=len)).changed(); ui.end_row();
            ui.label("Strand"); changed |= ui.checkbox(&mut comp, "Complement").changed(); ui.end_row();
            ui.label("Note"); changed |= ui.text_edit_multiline(&mut note).changed(); ui.end_row();
        });
        if changed {
            if !wraps { (f.start, f.end) = (lo - 1, hi); }
            f.strand = if comp { -1 } else { 1 };
            f.range = range_text(f);
            f.label = short(&f.name);
            if f.quals.iter().any(|q| q.starts_with("dnas_title=")) { set_qual(&mut f.quals, "dnas_title", &f.name.clone()); }
            set_qual(&mut f.quals, "note", &note);
            f.desc = if note.is_empty() { f.desc.clone() } else { note };
            f.dirty = true;
            d.lanes = assign_lanes(&mut d.feats, len);
        }
    }

    /// Replace bases [a, b) with `ins` and move features the way SeqBuilder does: coordinates past the edit
    /// shift, coordinates inside a deleted stretch clamp to its start, features left empty are deleted.
    fn edit_seq(&mut self, a: usize, b: usize, ins: &[u8]) {
        let Some(d) = &mut self.doc else { return };
        self.undo.push((d.seq.clone(), d.feats.clone()));
        self.redo.clear();
        let (cut, add) = (b - a, ins.len());
        d.seq.splice(a..b, ins.iter().copied());
        let start = |s: usize| if s >= b { s - cut + add } else if s >= a { a } else { s };
        let end = |e: usize| if e > b { e - cut + add } else if e > a { a } else { e };
        let mut gone = vec![];
        for (i, f) in d.feats.iter_mut().enumerate() {
            let (s, e) = (start(f.start), end(f.end));
            if (s, e) == (f.start, f.end) { continue; }
            if f.start < f.end && s >= e { gone.push(i); continue; }
            (f.start, f.end, f.dirty) = (s, e, true);
            if s < e { f.range = range_text(f); }
        }
        for i in gone.into_iter().rev() { d.deleted.extend(d.feats.remove(i).sbd_id); }
        d.seq_dirty = true;
        d.lanes = assign_lanes(&mut d.feats, d.seq.len());
        (self.sel, self.caret, self.sel_feat) = (None, Some(a + add), None);
        self.compute_hits();
    }

    fn undo_redo(&mut self, redo: bool) {
        let Some(d) = &mut self.doc else { return };
        let (from, to) = if redo { (&mut self.redo, &mut self.undo) } else { (&mut self.undo, &mut self.redo) };
        let Some((seq, feats)) = from.pop() else { return };
        to.push((std::mem::replace(&mut d.seq, seq), std::mem::replace(&mut d.feats, feats)));
        // Rewrite everything on the next save; ids of features that came back must not be deleted.
        d.feats.iter_mut().for_each(|f| f.dirty = true);
        d.deleted.retain(|id| !d.feats.iter().any(|f| f.sbd_id == Some(*id)));
        d.seq_dirty = true;
        d.lanes = assign_lanes(&mut d.feats, d.seq.len());
        (self.sel, self.caret, self.sel_feat) = (None, None, None);
        self.compute_hits();
    }

    /// Typing, deleting, cut and paste in the Sequence view.
    fn sequence_keys(&mut self, ctx: &egui::Context) {
        if self.view != View::Sequence || self.doc.is_none() || ctx.egui_wants_keyboard_input() { return; }
        let bases = |t: &str| t.bytes().filter(|b| b"ACGTUNRYKMSWBDHVacgtunrykmswbdhv".contains(b)).collect::<Vec<u8>>();
        let len = self.doc.as_ref().map_or(0, |d| d.seq.len());
        let span = |app: &App| app.sel.or(app.caret.map(|c| (c, c)));
        for ev in ctx.input(|i| i.events.clone()) {
            match ev {
                egui::Event::Text(t) | egui::Event::Paste(t) => {
                    let ins = bases(&t);
                    if let (false, Some((a, b))) = (ins.is_empty(), span(self)) { self.edit_seq(a, b, &ins); }
                }
                egui::Event::Cut => if let Some((a, b)) = self.sel {
                    if let Some(s) = self.selected() { ctx.copy_text(String::from_utf8_lossy(s).into_owned()); }
                    self.edit_seq(a, b, b"");
                },
                egui::Event::Key { key, pressed: true, modifiers, .. } if modifiers.is_none() => match (key, self.sel, self.caret) {
                    (egui::Key::Backspace | egui::Key::Delete, Some((a, b)), _) => self.edit_seq(a, b, b""),
                    (egui::Key::Backspace, None, Some(c)) if c > 0 => self.edit_seq(c - 1, c, b""),
                    (egui::Key::Delete, None, Some(c)) if c < len => self.edit_seq(c, c + 1, b""),
                    (egui::Key::ArrowLeft, _, Some(c)) => self.caret = Some(c.saturating_sub(1)),
                    (egui::Key::ArrowRight, _, Some(c)) => self.caret = Some((c + 1).min(len)),
                    _ => {}
                },
                _ => {}
            }
        }
    }

    /// File ▸ New DNA (⌘N): an empty, untitled linear sequence ready for typing.
    fn new_dna(&mut self) {
        self.doc = Some(Doc {
            name: "Untitled".into(), definition: String::new(), accession: String::new(), comment: String::new(),
            seq: vec![], circular: false, feats: vec![], lanes: 0, hits: vec![], path: "Untitled".into(),
            sbd: None, deleted: vec![], seq_dirty: false, header: String::new(), topology_dirty: false,
        });
        (self.view, self.sel, self.sel_feat, self.caret, self.error) = (View::Sequence, None, None, Some(0), None);
        (self.undo, self.redo) = (vec![], vec![]);
        self.compute_hits();
    }

    /// Sequence ▸ Linear / Circular.
    fn set_circular(&mut self, circular: bool) {
        let Some(d) = &mut self.doc else { return };
        if d.circular == circular { return; }
        (d.circular, d.topology_dirty) = (circular, true);
        if !circular && self.view == View::Circular { self.view = View::Linear; }
        self.compute_hits(); // sites across the origin come and go
    }

    /// Features ▸ New Feature (⌘K): a misc_feature over the selection, picked for editing in Details.
    fn new_feature(&mut self) {
        let (Some(d), Some((a, b))) = (&mut self.doc, self.sel) else { return };
        let name = "New Feature".to_string();
        d.feats.push(Feat {
            kind: "misc_feature".into(), label: name.clone(), desc: String::new(), range: format!("{}..{b}", a + 1),
            start: a, end: b, strand: 1, lane: 0, show: true, style: None, sbd_id: None, edited: false,
            base: TYPE_LOOKS.0.get("misc_feature").copied(), shape: None, base_shape: TYPE_LOOKS.1.get("misc_feature").copied(), shape_edited: false,
            quals: vec![format!("dnas_title=\"{name}\"")], dirty: true, name,
        });
        d.lanes = assign_lanes(&mut d.feats, d.seq.len());
        self.sel_feat = Some(d.feats.len() - 1);
        self.side_bottom = SideBottom::Details;
    }

    fn delete_feature(&mut self) {
        let (Some(d), Some(i)) = (&mut self.doc, self.sel_feat) else { return };
        if i >= d.feats.len() { return; }
        let f = d.feats.remove(i);
        d.deleted.extend(f.sbd_id);
        d.lanes = assign_lanes(&mut d.feats, d.seq.len());
        self.sel_feat = None;
    }

    /// Write the open .sbd back (to `to`, or where it came from) with this session's edits:
    /// deleted, new and changed features, then styles and Show flags.
    fn save(&mut self, to: Option<PathBuf>) {
        let Some(d) = &mut self.doc else { return };
        if d.sbd.is_none() {
            // A GenBank or FASTA document becomes an .sbd built from the template.
            let Some(path) = &to else { return self.save_as() };
            let name = path.file_stem().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
            let feats = d.feats.iter().map(|f| sbd_feat(f, 0)).collect();
            match sbd::new_document(SBD_TEMPLATE, &name, &d.header, d.seq.clone(), d.circular, feats) {
                Ok(s) => d.sbd = Some(s),
                Err(e) => { self.error = Some(e); return; }
            }
            for (i, f) in d.feats.iter_mut().enumerate() { (f.sbd_id, f.dirty) = (Some(i as u32), false); }
            (d.seq_dirty, d.deleted) = (false, vec![]);
        }
        let Some(s) = &mut d.sbd else { return };
        for id in d.deleted.drain(..) { s.remove_feature(id); }
        if d.seq_dirty {
            s.set_sequence(d.seq.clone());
            d.seq_dirty = false;
        }
        if d.topology_dirty {
            s.set_circular(d.circular);
            d.topology_dirty = false;
        }
        for f in d.feats.iter_mut().filter(|f| f.dirty) {
            let feat = sbd_feat(f, f.sbd_id.unwrap_or(0));
            match f.sbd_id {
                Some(_) => if let Err(e) = s.update_feature(feat) { self.error = Some(e); return; },
                None => f.sbd_id = Some(s.add_feature(feat)),
            }
            f.dirty = false;
        }
        for f in &d.feats {
            let Some(id) = f.sbd_id else { continue };
            if s.shown.get(&id) != Some(&f.show) { s.set_feature_shown(id, f.show); }
        }
        for f in d.feats.iter_mut().filter(|f| f.edited) {
            let Some(id) = f.sbd_id else { continue };
            s.set_feature_style(id, f.style);
            f.edited = false;
        }
        for f in d.feats.iter_mut().filter(|f| f.shape_edited) {
            let (Some(id), Some(sh)) = (f.sbd_id, f.shape) else { continue };
            s.set_feature_shape(id, sh);
            f.shape_edited = false;
            f.edited = false;
        }
        let path = to.unwrap_or_else(|| d.path.clone());
        match std::fs::write(&path, s.to_bytes()) {
            Ok(()) => { d.path = path; self.error = None; }
            Err(e) => self.error = Some(format!("{}: {e}", path.display())),
        }
    }

    /// File ▸ Export GenBank…
    fn export(&mut self) {
        let Some(d) = &self.doc else { return };
        let name = d.path.file_stem().map(|n| n.to_string_lossy().into_owned() + ".gb").unwrap_or_default();
        if let Some(p) = rfd::FileDialog::new().add_filter("GenBank", &["gb", "gbk"]).set_file_name(name).save_file() {
            if let Err(e) = std::fs::write(&p, export_genbank(d)) { self.error = Some(format!("{}: {e}", p.display())); }
        }
    }

    fn save_as(&mut self) {
        let name = self.doc.as_ref().and_then(|d| d.path.file_stem()).map(|n| n.to_string_lossy().into_owned() + ".sbd").unwrap_or_default();
        if let Some(p) = rfd::FileDialog::new().add_filter("SeqBuilder document", &["sbd"]).set_file_name(name).save_file() {
            self.save(Some(p));
        }
    }

    fn load_enzymes(&mut self, ezd: &Path, sel: &Path) {
        match std::fs::read(ezd) {
            Ok(b) => self.enz = enzymes::parse_ezd(&b),
            Err(e) => self.error = Some(format!("{}: {e}", ezd.display())),
        }
        match std::fs::read_to_string(sel).map_err(|e| e.to_string()).and_then(|t| enzymes::parse_sel(&t)) {
            Ok(s) => self.sels = s,
            Err(e) => self.error = Some(format!("{}: {e}", sel.display())),
        }
        self.shown.clear();
        self.compute_hits();
    }

    fn compute_hits(&mut self) {
        if let Some(d) = &mut self.doc {
            d.hits = self.enz.iter().map(|e| find_sites(&d.seq, d.circular, &e.site, self.potential)).collect();
        }
        self.refresh_sites();
    }

    /// Enzymes passing selector `id` on the current sequence.
    fn selector_result(&self, id: u32) -> Vec<usize> {
        let Some(d) = &self.doc else { return vec![] };
        (0..self.enz.len()).filter(|&i| enzymes::selects(&self.sels, id, &self.enz[i], &d.hits[i])).collect()
    }

    /// Check the selector's enzymes that cut; SeqBuilder leaves absent ones unchecked.
    fn apply_selector(&mut self, id: u32) {
        let Some(d) = &self.doc else { return };
        self.shown = self.selector_result(id).into_iter().filter(|&i| !d.hits[i].is_empty()).collect();
        self.refresh_sites();
    }

    fn refresh_sites(&mut self) {
        self.sites.clear();
        let Some(d) = &self.doc else { return };
        for &i in &self.shown {
            self.sites.extend(d.hits[i].iter().map(|&pos| Site { enz: i, pos }));
        }
        self.sites.sort_by_key(|s| (s.pos, s.enz));
    }

    fn selected(&self) -> Option<&[u8]> {
        let (a, b) = self.sel?;
        Some(&self.doc.as_ref()?.seq[a..b])
    }

    fn menu(&mut self, ui: &mut egui::Ui) {
        egui::MenuBar::new().ui(ui, |ui| {
            ui.menu_button("File", |ui| {
                if ui.button("New DNA  ⌘N").clicked() { self.new_dna(); }
                if ui.button("Open…  ⌘O").clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("Sequence", &["sbd", "gb", "gbk", "genbank", "fasta", "fa", "fas", "seq"]).pick_file() {
                        self.open(&p);
                    }
                }
                if ui.button("Save  ⌘S").clicked() { self.save(None); }
                if ui.button("Save As…").clicked() { self.save_as(); }
                if ui.add_enabled(self.doc.is_some(), egui::Button::new("Export GenBank…")).clicked() { self.export(); }
                if ui.button("Open Enzyme Files…").clicked() {
                    let ezd = rfd::FileDialog::new().add_filter("Enzyme file", &["ezd"]).pick_file();
                    let sel = ezd.as_ref().and_then(|_| rfd::FileDialog::new().add_filter("Selector file", &["sel"]).pick_file());
                    if let (Some(e), Some(s)) = (ezd, sel) { self.load_enzymes(&e, &s); }
                }
            });
            ui.menu_button("Edit", |ui| {
                if ui.add_enabled(!self.undo.is_empty(), egui::Button::new("Undo  ⌘Z")).clicked() { self.undo_redo(false); }
                if ui.add_enabled(!self.redo.is_empty(), egui::Button::new("Redo  ⇧⌘Z")).clicked() { self.undo_redo(true); }
                ui.separator();
                if ui.button("Copy  ⌘C").clicked() {
                    if let Some(s) = self.selected() { ui.ctx().copy_text(String::from_utf8_lossy(s).into_owned()); }
                }
                if ui.button("Select All  ⌘A").clicked() {
                    if let Some(d) = &self.doc { self.sel = Some((0, d.seq.len())); }
                }
                if ui.button("Go to Position…  ⌘J").clicked() { self.goto_open = true; }
            });
            ui.menu_button("Features", |ui| {
                if ui.add_enabled(self.sel.is_some(), egui::Button::new("New Feature  ⌘K")).clicked() { self.new_feature(); }
                if ui.add_enabled(self.sel_feat.is_some(), egui::Button::new("Edit Feature…")).clicked() { self.edit_open = true; }
                if ui.add_enabled(self.sel_feat.is_some(), egui::Button::new("Delete Feature")).clicked() { self.delete_feature(); }
            });
            ui.menu_button("Sequence", |ui| {
                let circular = self.doc.as_ref().map(|d| d.circular);
                if ui.add_enabled(circular.is_some(), egui::Button::selectable(circular == Some(false), "Linear")).clicked() { self.set_circular(false); }
                if ui.add_enabled(circular.is_some(), egui::Button::selectable(circular == Some(true), "Circular")).clicked() { self.set_circular(true); }
            });
            ui.menu_button("Enzymes", |ui| {
                if ui.button("Selector Manager…").clicked() { self.selmgr = self.sels.first().map(|s| s.id); }
                ui.menu_button("Apply Selector", |ui| {
                    for (id, name) in self.sels.iter().map(|s| (s.id, s.name.trim().to_string())).collect::<Vec<_>>() {
                        if ui.button(name).clicked() { self.apply_selector(id); }
                    }
                });
            });
            ui.menu_button("Format", |ui| {
                ui.label("Residue Grouping");
                ui.radio_value(&mut self.grouping, 0, "None");
                ui.radio_value(&mut self.grouping, 3, "Group by 3");
                ui.radio_value(&mut self.grouping, 10, "Group by 10");
            });
        });
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::S)) { self.save(None); }
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::K)) { self.new_feature(); }
        if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::N)) { self.new_dna(); }
        if !ctx.egui_wants_keyboard_input() {
            if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND | egui::Modifiers::SHIFT, egui::Key::Z)) { self.undo_redo(true); }
            if ctx.input_mut(|i| i.consume_key(egui::Modifiers::COMMAND, egui::Key::Z)) { self.undo_redo(false); }
        }
        self.sequence_keys(ctx);
        let (open, all, goto) = ctx.input_mut(|i| (
            i.consume_key(egui::Modifiers::COMMAND, egui::Key::O),
            i.consume_key(egui::Modifiers::COMMAND, egui::Key::A),
            i.consume_key(egui::Modifiers::COMMAND, egui::Key::J),
        ));
        if open {
            if let Some(p) = rfd::FileDialog::new().pick_file() { self.open(&p); }
        }
        if all { if let Some(d) = &self.doc { self.sel = Some((0, d.seq.len())); } }
        if goto { self.goto_open = true; }
        // egui turns ⌘C into an event rather than a key press.
        if ctx.input(|i| i.events.iter().any(|e| matches!(e, egui::Event::Copy))) {
            if let Some(s) = self.selected() { ctx.copy_text(String::from_utf8_lossy(s).into_owned()); }
        }
        if let Some(f) = ctx.input(|i| i.raw.dropped_files.first().map(|f| f.path().to_path_buf())) { self.open(&f); }
    }

    /// SeqBuilder's tool strip, top to bottom; tools not made yet are shown disabled.
    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
            let toggle = |ui: &mut egui::Ui, on: &mut bool, icon: &str, tip: &str| {
                if ui.add(egui::Button::selectable(*on, icon)).on_hover_text(tip).clicked() { *on = !*on; }
            };
            let mut enzymes = !self.hide_enzymes;
            let mut rulers = !self.hide_rulers;
            let mut side = !self.hide_side;
            toggle(ui, &mut self.complement, "⇄", "Strand");
            toggle(ui, &mut enzymes, "✂", "Enzymes");
            toggle(ui, &mut self.translations, "Ⓣ", "Translations and ORFs");
            toggle(ui, &mut rulers, "📏", "Rulers");
            toggle(ui, &mut self.spacing, "↕", "Vertical Spacing");
            toggle(ui, &mut side, "◫", "Show/Hide Side Panels");
            (self.hide_enzymes, self.hide_rulers, self.hide_side) = (!enzymes, !rulers, !side);
            ui.separator();
            for (icon, tip) in [("▤", "Open Agarose Gel"), ("★", "Auto-annotate this sequence"), ("📚", "Open Feature Library"),
                ("📁", "Vector Catalog"), ("ELN", "Send document to selected SciNote sub-folder")] {
                ui.add_enabled(false, egui::Button::new(icon)).on_disabled_hover_text(format!("{tip} (not implemented yet)"));
            }
            ui.separator();
            if ui.button("▦").on_hover_text("Compact EditSeq Layout").clicked() { self.layout(View::Sequence, false, true, true); }
            if ui.button("☰").on_hover_text("Full Sequence Layout").clicked() { self.layout(View::Sequence, true, true, false); }
            ui.add_enabled(false, egui::Button::new("⚟")).on_disabled_hover_text("Primer Design Layout (not implemented yet)");
            if ui.button("↺").on_hover_text("Default Layout").clicked() {
                let circular = self.doc.as_ref().is_some_and(|d| d.circular);
                self.layout(if circular { View::Circular } else { View::Linear }, false, false, false);
            }
            ui.separator();
            ui.add(egui::Button::selectable(true, "↖")).on_hover_text("Range");
            if ui.button("🔍+").on_hover_text("Zoom In").clicked() { self.zoom = (self.zoom * 1.25).min(8.0); }
            if ui.button("🔍−").on_hover_text("Zoom Out").clicked() { self.zoom = (self.zoom / 1.25).max(0.25); }
        });
    }

    /// The layout buttons: which view, which strands, which panels.
    fn layout(&mut self, view: View, complement: bool, hide_table: bool, hide_side: bool) {
        (self.view, self.complement, self.hide_table, self.hide_side) = (view, complement, hide_table, hide_side);
    }

    fn info_bar(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            match self.sel {
                Some((a, b)) => { ui.small(format!("Selection: Forward  {}..{}, length = {}", a + 1, b, b - a)); }
                None if self.view == View::Sequence => { ui.small("Position: 1"); }
                None => { ui.small(format!("Size: {:.1}x    Position: 1", self.zoom)); }
            }
            if let Some(e) = &self.error { ui.colored_label(Color32::RED, e); }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if let Some(d) = &self.doc { ui.small(size_label(d.seq.len())); }
            });
        });
    }

    fn feature_table(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("table_tabs").show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.tab, Tab::Features, "Features");
                ui.selectable_value(&mut self.tab, Tab::Primers, "Primers");
                ui.selectable_value(&mut self.tab, Tab::Comment, "Comment");
                ui.selectable_value(&mut self.tab, Tab::Minimap, "Minimap");
                ui.selectable_value(&mut self.tab, Tab::SiteSummary, "Site Summary");
                todo_tab(ui, "History");
            });
        });
        let Some(doc) = &mut self.doc else { return };
        let mut pick = None;
        let mut jump = None;
        egui::ScrollArea::both().auto_shrink(false).show(ui, |ui| match self.tab {
            Tab::Features => {
                let mut toggled = false;
                egui::Grid::new("feats").striped(true).num_columns(7).spacing([16.0, 4.0]).show(ui, |ui| {
                    for h in ["Show", "Type", "Name", "Range", "Strand", "Length", "Description"] { ui.strong(h); }
                    ui.end_row();
                    for (i, f) in doc.feats.iter_mut().enumerate() {
                        toggled |= ui.checkbox(&mut f.show, "").changed();
                        ui.label(&f.kind);
                        let name = ui.selectable_label(self.sel_feat == Some(i), &f.label);
                        if name.clicked() { pick = Some(i); }
                        if name.double_clicked() { self.edit_open = true; }
                        ui.label(&f.range);
                        ui.label(if f.strand < 0 { "◀" } else { "▶" });
                        ui.label((if f.start < f.end { f.end - f.start } else { doc.seq.len() - f.start + f.end }).to_string());
                        ui.label(&f.desc);
                        ui.end_row();
                    }
                });
                if toggled { doc.lanes = assign_lanes(&mut doc.feats, doc.seq.len()); }
            }
            Tab::Comment => { ui.label(if doc.comment.is_empty() { "No comment." } else { &doc.comment }); }
            Tab::Primers => {
                // SeqBuilder lists PCR_primer features; Score, dTm and Ta OP need its primer design, not made yet.
                egui::Grid::new("primers").striped(true).num_columns(11).spacing([14.0, 4.0]).show(ui, |ui| {
                    for h in ["Show", "Active", "Set Name", "Name", "F/R", "Primer Sequence", "Score", "Length", "Tm", "dTm", "Ta OP"] { ui.strong(h); }
                    ui.end_row();
                    for (i, f) in doc.feats.iter_mut().enumerate().filter(|(_, f)| f.kind == "PCR_primer") {
                        let primer = primer_sequence(f, &doc.seq);
                        ui.checkbox(&mut f.show, "");
                        ui.label("");
                        ui.label(f.quals.iter().find_map(|q| q.strip_prefix("pair=")).unwrap_or("").trim_matches('"'));
                        if ui.selectable_label(self.sel_feat == Some(i), &f.label).clicked() { pick = Some(i); }
                        ui.label(if f.strand < 0 { "R" } else { "F" });
                        ui.monospace(String::from_utf8_lossy(&primer).to_uppercase());
                        ui.label("");
                        ui.label(primer.len().to_string());
                        ui.label(format!("{:.1}", primer_tm(&primer)));
                        ui.label("");
                        ui.label("");
                        ui.end_row();
                    }
                });
            }
            Tab::Minimap => {
                // One row per checked enzyme: name, number of sites, the sites along the whole sequence.
                let mut shown: Vec<usize> = self.shown.iter().copied().collect();
                shown.sort_by(|&a, &b| self.enz[a].name.cmp(&self.enz[b].name));
                let (row_h, name_w, freq_w) = (24.0, 70.0, 40.0);
                let size = Vec2::new(ui.available_width(), 22.0 + row_h * shown.len() as f32);
                let (rect, resp) = ui.allocate_exact_size(size, Sense::click());
                let p = ui.painter_at(rect);
                let v = ui.visuals();
                let (fg, track_bg, site) = (v.text_color(), v.faint_bg_color, Color32::from_rgb(40, 70, 220));
                let (x0, x1) = (rect.left() + name_w + freq_w, rect.right() - 12.0);
                let len = doc.seq.len().max(1);
                let x = |pos: usize| x0 + pos as f32 / len as f32 * (x1 - x0);
                let top = rect.top() + 20.0;
                p.text(Pos2::new(rect.left() + 2.0, top - 2.0), Align2::LEFT_BOTTOM, "Name", FontId::proportional(12.0), fg);
                p.text(Pos2::new(rect.left() + name_w + 2.0, top - 2.0), Align2::LEFT_BOTTOM, "Freq", FontId::proportional(12.0), fg);
                p.line_segment([Pos2::new(x0, top), Pos2::new(x1, top)], Stroke::new(1.5, fg));
                let step = tick_step(len as f32);
                for pos in (step..len).step_by(step) {
                    p.line_segment([Pos2::new(x(pos), top - 4.0), Pos2::new(x(pos), top)], Stroke::new(1.0, fg));
                    p.text(Pos2::new(x(pos), top - 4.0), Align2::CENTER_BOTTOM, pos.to_string(), FontId::proportional(9.0), fg);
                }
                for (k, &e) in shown.iter().enumerate() {
                    let y = top + 2.0 + row_h * k as f32;
                    let mid = y + row_h / 2.0;
                    p.text(Pos2::new(rect.left() + 2.0, mid), Align2::LEFT_CENTER, &self.enz[e].name, FontId::proportional(12.0), fg);
                    p.text(Pos2::new(rect.left() + name_w + freq_w / 2.0, mid), Align2::CENTER_CENTER, doc.hits[e].len().to_string(), FontId::proportional(12.0), fg);
                    p.rect_filled(Rect::from_x_y_ranges(x0..=x1, y + 2.0..=y + row_h - 2.0), 0.0, track_bg);
                    p.line_segment([Pos2::new(x0, mid), Pos2::new(x1, mid)], Stroke::new(1.0, fg));
                    for &pos in &doc.hits[e] {
                        p.line_segment([Pos2::new(x(pos), y + 2.0), Pos2::new(x(pos), y + row_h - 2.0)], Stroke::new(1.0, site));
                    }
                }
                // a click on a track goes there in the sequence
                if let Some(pp) = resp.interact_pointer_pos().filter(|pp| resp.clicked() && pp.x >= x0) {
                    jump = Some((((pp.x - x0) / (x1 - x0)) * len as f32).clamp(0.0, len as f32 - 1.0) as usize);
                }
            }
            Tab::SiteSummary => {
                egui::Grid::new("sites").striped(true).num_columns(3).spacing([16.0, 4.0]).show(ui, |ui| {
                    for h in ["Enzyme", "Sites", "Positions", "Recognition"] { ui.strong(h); }
                    ui.end_row();
                    let mut shown: Vec<usize> = self.shown.iter().copied().collect();
                    shown.sort_by(|&a, &b| self.enz[a].name.cmp(&self.enz[b].name));
                    for i in shown {
                        let hits = &doc.hits[i];
                        ui.label(&self.enz[i].name);
                        ui.label(hits.len().to_string());
                        ui.label(hits.iter().take(8).map(|p| (p + 1).to_string()).collect::<Vec<_>>().join(", "));
                        ui.monospace(&self.enz[i].rs);
                        ui.end_row();
                    }
                });
            }
        });
        if let Some(pos) = jump {
            (self.view, self.caret, self.scroll_to, self.sel) = (View::Sequence, Some(pos), Some(pos), None);
        }
        if let Some(i) = pick {
            let f = &doc.feats[i];
            if f.start < f.end {
                self.sel = Some((f.start, f.end));
                self.scroll_to = Some(f.start);
            }
            self.sel_feat = Some(i);
            self.side_top = SideTop::Style;
        }
    }

    fn side_panel(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("side_bottom").resizable(true).default_size(320.0).show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.side_bottom, SideBottom::Settings, "Settings");
                ui.selectable_value(&mut self.side_bottom, SideBottom::Details, "Details");
            });
            ui.separator();
            match self.side_bottom {
                SideBottom::Settings => {
                    ui.checkbox(&mut self.complement, "Show complementary strand");
                    ui.label("Residue grouping");
                    ui.radio_value(&mut self.grouping, 0, "None");
                    ui.radio_value(&mut self.grouping, 3, "Group by 3");
                    ui.radio_value(&mut self.grouping, 10, "Group by 10");
                }
                SideBottom::Details if self.sel_feat.is_some() && self.doc.is_some() => self.feature_details(ui),
                SideBottom::Details => match &self.doc {
                    None => { ui.label("No document open."); }
                    Some(_) if self.sel.is_some() => {
                        let (a, b) = self.sel.unwrap();
                        ui.heading("Region selection");
                        ui.separator();
                        ui.label(format!("Selection: Forward {}..{}, length = {}", a + 1, b, b - a));
                        if let Some(s) = self.selected() { statistics(ui, s); }
                    }
                    Some(d) => {
                        ui.heading(&d.name);
                        if !d.definition.is_empty() { ui.label(&d.definition); }
                        egui::Grid::new("details").num_columns(2).show(ui, |ui| {
                            for (k, v) in [
                                ("Accession", d.accession.clone()), ("Length", format!("{} bp", d.seq.len())),
                                ("Topology", (if d.circular { "circular" } else { "linear" }).into()),
                                ("GC content", format!("{:.1}%", gc_percent(&d.seq))), ("Features", d.feats.len().to_string()),
                            ] { ui.label(k); ui.label(v); ui.end_row(); }
                        });
                    }
                },
            }
        });
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.side_top, SideTop::Style, "Style");
            ui.selectable_value(&mut self.side_top, SideTop::Enzymes, "Enzymes");
        });
        ui.separator();
        match self.side_top {
            SideTop::Style => {
                let feat = self.sel_feat.and_then(|i| self.doc.as_mut()?.feats.get_mut(i));
                let Some(f) = feat else { ui.label("Select a feature in the table."); return };
                // Laid out as SeqBuilder's styleInfo panel (SeqBuilder.xrc).
                ui.strong(&f.label);
                let mut st = look(f);
                let mut changed = false;
                let colour = |ui: &mut egui::Ui, on: &mut Option<[u8; 3]>, label: &str, default: [u8; 3]| {
                    let mut enabled = on.is_some();
                    let mut c = on.unwrap_or(default);
                    let mut ch = ui.checkbox(&mut enabled, label).changed();
                    ch |= ui.color_edit_button_srgb(&mut c).changed();
                    if ch { *on = enabled.then_some(c); }
                    ch
                };
                egui::Grid::new("style").num_columns(4).show(ui, |ui| {
                    let last_fill = st.fill.unwrap_or([255, 255, 255]);
                    changed |= colour(ui, &mut st.line, "Line:", [0, 0, 0]);
                    ui.label("Weight:");
                    egui::ComboBox::from_id_salt("weight").selected_text(format!("{}", st.weight)).width(50.0).show_ui(ui, |ui| {
                        for w in [0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0] { changed |= ui.selectable_value(&mut st.weight, w, format!("{w}")).changed(); }
                    });
                    ui.end_row();
                    changed |= colour(ui, &mut st.fill, "Fill:", last_fill);
                    ui.end_row();
                    changed |= colour(ui, &mut st.shadow, "Shadow:", [192, 192, 192]);
                    ui.end_row();
                });
                if changed { (f.style, f.edited) = (Some(st), true); }
                ui.add_space(6.0);
                ui.small("Feature");
                let mut sh = shape(f);
                let mut ch = false;
                egui::Grid::new("feature_shape").num_columns(2).show(ui, |ui| {
                    for (label, field, names) in [("Shape:", &mut sh.shape, &SHAPES[..]), ("Linkage:", &mut sh.linkage, &LINKAGES[..]), ("Label Position:", &mut sh.label, &LABEL_POSITIONS[..])] {
                        ui.label(label);
                        let current = names.get((*field as usize).wrapping_sub(1)).copied().unwrap_or("?");
                        egui::ComboBox::from_id_salt(label).selected_text(current).show_ui(ui, |ui| {
                            for (k, n) in names.iter().enumerate() { ch |= ui.selectable_value(field, k as u32 + 1, *n).changed(); }
                        });
                        ui.end_row();
                    }
                });
                ch |= ui.checkbox(&mut sh.attached, "Label Attached").changed();
                if ch { (f.shape, f.shape_edited) = (Some(sh), true); }
                ui.add_space(6.0);
                ui.add_enabled_ui(false, |ui| {
                    ui.small("Font");
                    ui.horizontal(|ui| { ui.label("Selected Item"); ui.label("Helvetica 12"); });
                }).response.on_disabled_hover_text("Not implemented yet");
            }
            SideTop::Enzymes => {
                if self.enz.is_empty() { ui.label("No enzymes loaded (File ▸ Open Enzyme Files…)."); return; }
                ui.horizontal(|ui| {
                    ui.add(egui::TextEdit::singleline(&mut self.enz_filter).hint_text("🔍 Filter").desired_width(ui.available_width() - 60.0));
                    if ui.button("Clear").clicked() { self.shown.clear(); self.refresh_sites(); }
                });
                let filter = self.enz_filter.to_lowercase();
                let list: Vec<usize> = (0..self.enz.len()).filter(|&i| self.enz[i].name.to_lowercase().contains(&filter)).collect();
                let cols = ((ui.available_width() / 92.0) as usize).max(1);
                let mut changed = false;
                egui::ScrollArea::vertical().auto_shrink(false).max_height(ui.available_height() - 52.0).show(ui, |ui| {
                    egui::Grid::new("enzymes").num_columns(cols).min_col_width(84.0).show(ui, |ui| {
                        for (n, &i) in list.iter().enumerate() {
                            let e = &self.enz[i];
                            let mut on = self.shown.contains(&i);
                            let text = if on { egui::RichText::new(&e.name).strong() } else { egui::RichText::new(&e.name) };
                            let hits = self.doc.as_ref().map_or(0, |d| d.hits[i].len());
                            if ui.checkbox(&mut on, text).on_hover_text(format!("{}\n{} sites\nIsoschizomers: {}\nVendors: {}", e.rs, hits, e.iso, e.vendors)).changed() {
                                if on { self.shown.insert(i); } else { self.shown.remove(&i); }
                                changed = true;
                            }
                            if (n + 1) % cols == 0 { ui.end_row(); }
                        }
                    });
                });
                ui.separator();
                ui.small("Enzyme Cuts");
                if ui.checkbox(&mut self.potential, "Show Potential Cuts at Ambiguous Bases").changed() { self.compute_hits(); }
                if changed { self.refresh_sites(); }
            }
        }
    }

    /// Read-only view of the selectors (SeqBuilder: Enzymes ▸ Selector Manager); editing comes later.
    fn selector_manager(&mut self, ctx: &egui::Context) {
        let Some(mut picked) = self.selmgr else { return };
        let mut open = true;
        let mut apply = false;
        egui::Window::new("Enzyme Selector Manager").open(&mut open).default_size([760.0, 420.0]).show(ctx, |ui| {
            ui.columns(3, |c| {
                c[0].strong("Selector");
                egui::ScrollArea::vertical().id_salt("sels").show(&mut c[0], |ui| {
                    for s in &self.sels { ui.selectable_value(&mut picked, s.id, s.name.trim()); }
                });
                c[1].strong("Definition");
                if let Some(s) = self.sels.iter().find(|s| s.id == picked) {
                    c[1].label(format!("Name: {}", s.name.trim()));
                    for line in enzymes::describe(&self.sels, s) { c[1].label(line); }
                }
                apply = c[1].button("Show these enzymes").clicked();
                let result = self.selector_result(picked);
                c[2].strong(format!("Result ({})", result.len()));
                egui::ScrollArea::vertical().id_salt("result").show(&mut c[2], |ui| {
                    if self.doc.is_none() { ui.label("Open a sequence to evaluate."); }
                    for i in result { ui.label(&self.enz[i].name); }
                });
            });
        });
        if apply { self.apply_selector(picked); }
        self.selmgr = open.then_some(picked);
    }

    /// White page centred in the view, scaled by the zoom factor; scrolls when zoomed in.
    fn page(&self, ui: &mut egui::Ui, draw: impl FnOnce(&egui::Painter, Rect, &egui::Ui)) {
        let avail = ui.available_size();
        egui::ScrollArea::both().auto_shrink(false).show(ui, |ui| {
            let (rect, _) = ui.allocate_exact_size(avail * self.zoom.max(1.0), Sense::hover());
            let side = rect.width().min(rect.height()) * self.zoom.min(1.0);
            let page = Rect::from_center_size(rect.center(), Vec2::new(rect.width().min(side * 1.25), side).max(Vec2::splat(200.0)) - Vec2::splat(16.0));
            let p = ui.painter_at(rect);
            p.rect(page, 0.0, ui.visuals().extreme_bg_color, Stroke::new(1.0, ui.visuals().weak_text_color()), egui::StrokeKind::Inside);
            draw(&p, page, ui);
        });
    }

    fn circular_map(&self, p: &egui::Painter, rect: Rect, ui: &egui::Ui, doc: &Doc) {
        let c = rect.center();
        let r = (rect.width().min(rect.height()) / 2.0 - 80.0).max(60.0);
        let len = doc.seq.len().max(1) as f32;
        let ang = |pos: usize| pos as f32 / len * TAU - PI / 2.0;
        let at = |a: f32, rad: f32| c + Vec2::angled(a) * rad;
        let fg = ui.visuals().text_color();
        let arc = |a0: f32, a1: f32, rad: f32| -> Vec<Pos2> {
            let n = (((a1 - a0).abs() * rad / 4.0) as usize).max(2);
            (0..=n).map(|i| at(a0 + (a1 - a0) * i as f32 / n as f32, rad)).collect()
        };
        let span = |s: usize, e: usize| { let a0 = ang(s); let mut a1 = ang(e); if a1 <= a0 { a1 += TAU; } (a0, a1) };
        let small = FontId::proportional(9.0);

        // backbone, selection, ticks
        p.circle_stroke(c, r, Stroke::new(3.0, fg));
        if let Some((s, e)) = self.sel {
            let (a0, a1) = span(s, e);
            p.add(Shape::line(arc(a0, a1, r), Stroke::new(7.0, ui.visuals().selection.bg_fill)));
        }
        let step = tick_step(len);
        for pos in (step..doc.seq.len()).step_by(step) {
            let a = ang(pos);
            p.line_segment([at(a, r - 5.0), at(a, r)], Stroke::new(1.0, fg));
            p.text(at(a, r - 12.0), Align2::CENTER_CENTER, pos.to_string(), small.clone(), fg);
        }

        // features inside the backbone, as outlined arrows
        let mut placed: Vec<Rect> = vec![];
        for f in doc.feats.iter().filter(|f| f.show) {
            let rad = r - 28.0 - f.lane as f32 * 14.0;
            if rad < 20.0 { continue; }
            let (a0, a1) = span(f.start, f.end);
            let (fill, line) = colors(f);
            let st = look(f);
            let outline = Stroke::new(st.weight as f32, line);
            if (a1 - a0) * rad < 6.0 {
                // a few bases: a radial mark, since thick outlines on tiny shapes spike out
                let mid = (a0 + a1) / 2.0;
                let mark = if fill == Color32::TRANSPARENT { line } else { fill };
                p.line_segment([at(mid, rad - 6.0), at(mid, rad + 6.0)], Stroke::new(3.0, mark));
                continue;
            }
            let sh = shape(f);
            let fill = if sh.shape == 5 { Color32::TRANSPARENT } else { fill }; // Stroked Arrow: outline only
            let head = (8.0 / rad).min((a1 - a0) / 2.0);
            // arrow heads as (tip, back) angles
            let heads: Vec<(f32, f32)> = match sh.shape {
                2 | 5 => vec![if f.strand < 0 { (a0, a0 + head) } else { (a1, a1 - head) }],
                3 => vec![(a0, a0 + head), (a1, a1 - head)],
                _ => vec![],
            };
            let b0 = if heads.iter().any(|h| h.0 == a0) { a0 + head } else { a0 };
            let b1 = if heads.iter().any(|h| h.0 == a1) { a1 - head } else { a1 };
            let off = |v: &Vec<Pos2>| v.iter().map(|q| *q + Vec2::splat(2.0)).collect::<Vec<_>>();
            let mid = (a0 + a1) / 2.0;
            match sh.shape {
                4 => {
                    // Triangle Glyph: one triangle at the middle, pointing along the strand
                    let h = (6.0 / rad).min((a1 - a0) / 2.0).max(4.0 / rad);
                    let (tip, back) = if f.strand < 0 { (mid - h, mid + h) } else { (mid + h, mid - h) };
                    let tri = vec![at(tip, rad), at(back, rad + 7.0), at(back, rad - 7.0)];
                    if let Some(c) = st.shadow { p.add(Shape::convex_polygon(off(&tri), rgb(c), Stroke::NONE)); }
                    p.add(Shape::convex_polygon(tri, fill, Stroke::new(outline.width.min(2.0), line)));
                }
                6 => {
                    // Stroked Goalpost: the span as a line with posts at both ends
                    let colour = if line == Color32::TRANSPARENT { fill } else { line };
                    let stroke = Stroke::new((st.weight as f32).max(1.5), colour);
                    p.add(Shape::line(arc(a0, a1, rad), stroke));
                    for a in [a0, a1] { p.line_segment([at(a, rad - 5.0), at(a, rad + 5.0)], stroke); }
                }
                _ => {
                    let arrows: Vec<Vec<Pos2>> = heads.iter().map(|&(tip, back)| vec![at(tip, rad), at(back, rad + 7.0), at(back, rad - 7.0)]).collect();
                    if let Some(c) = st.shadow {
                        if b1 > b0 { p.add(Shape::line(off(&arc(b0, b1, rad)), Stroke::new(8.0, rgb(c)))); }
                        for a in &arrows { p.add(Shape::convex_polygon(off(a), rgb(c), Stroke::NONE)); }
                    }
                    if b1 > b0 {
                        let mut band = arc(b0, b1, rad + 4.0);
                        band.extend(arc(b1, b0, rad - 4.0));
                        band.push(band[0]);
                        p.add(Shape::line(arc(b0, b1, rad), Stroke::new(8.0, fill)));
                        p.add(Shape::line(band, outline));
                    }
                    for a in arrows { p.add(Shape::convex_polygon(a, fill, Stroke::new(outline.width.min(2.0), line))); } // sharp tips
                }
            }
            let galley = p.layout_no_wrap(f.label.clone(), small.clone(), fg);
            // Label Position: Above = outward, Inside = on the band when it fits, Below = toward the centre
            let fits = galley.size().x + 4.0 < (b1 - b0) * rad && sh.shape != 4 && sh.shape != 6;
            let lr = match sh.label { 1 => rad + 13.0, 2 if fits => rad, _ => rad - 13.0 };
            let lrect = Align2::CENTER_CENTER.anchor_size(at(mid, lr), galley.size());
            if placed.iter().any(|q| q.intersects(lrect)) || lrect.distance_to_pos(c) < 30.0 { continue; }
            placed.push(lrect);
            p.galley(lrect.min, galley, fg);
        }

        // restriction sites outside, labels stacked to avoid overlap
        let mut placed: Vec<Rect> = vec![];
        for s in &self.sites {
            let a = ang(s.pos);
            let right = Vec2::angled(a).x >= 0.0;
            let align = if right { Align2::LEFT_CENTER } else { Align2::RIGHT_CENTER };
            let galley = p.layout_no_wrap(self.enz[s.enz].name.clone(), small.clone(), fg);
            let base = at(a, r + 28.0);
            let Some(lrect) = (0..12u32).map(|k| {
                let dy = 11.0 * k.div_ceil(2) as f32 * if k % 2 == 0 { 1.0 } else { -1.0 };
                align.anchor_size(base + Vec2::new(0.0, dy), galley.size())
            }).find(|lr| !placed.iter().any(|q| q.intersects(*lr))) else { continue };
            placed.push(lrect);
            let anchor = if right { lrect.left_center() } else { lrect.right_center() };
            p.line_segment([at(a, r + 2.0), at(a, r + 14.0)], Stroke::new(0.6, fg));
            p.line_segment([at(a, r + 14.0), anchor], Stroke::new(0.6, fg));
            p.galley(lrect.min, galley, fg);
        }

        p.text(c - Vec2::new(0.0, 4.0), Align2::CENTER_BOTTOM, &doc.name, FontId::proportional(15.0), fg);
        p.text(c + Vec2::new(0.0, 4.0), Align2::CENTER_TOP, format!("{} bp", doc.seq.len()), FontId::proportional(13.0), fg);
    }

    fn linear_map(&self, p: &egui::Painter, rect: Rect, ui: &egui::Ui, doc: &Doc) {
        let fg = ui.visuals().text_color();
        let small = FontId::proportional(9.0);
        let (x0, x1) = (rect.left() + 30.0, rect.right() - 30.0);
        let y = rect.top() + rect.height() * 0.45;
        let len = doc.seq.len().max(1) as f32;
        let x = |pos: usize| x0 + pos as f32 / len * (x1 - x0);
        p.line_segment([Pos2::new(x0, y), Pos2::new(x1, y)], Stroke::new(3.0, fg));
        if let Some((s, e)) = self.sel {
            p.line_segment([Pos2::new(x(s), y), Pos2::new(x(e), y)], Stroke::new(7.0, ui.visuals().selection.bg_fill));
        }
        let step = tick_step(len);
        for pos in (step..doc.seq.len()).step_by(step) {
            p.line_segment([Pos2::new(x(pos), y), Pos2::new(x(pos), y + 5.0)], Stroke::new(1.0, fg));
            p.text(Pos2::new(x(pos), y + 7.0), Align2::CENTER_TOP, pos.to_string(), small.clone(), fg);
        }
        for f in doc.feats.iter().filter(|f| f.show) {
            let (fill, line) = colors(f);
            let st = look(f);
            let fy = y + 28.0 + f.lane as f32 * 16.0;
            let spans = if f.start < f.end { vec![(f.start, f.end)] } else { vec![(f.start, doc.seq.len()), (0, f.end)] };
            let sh = shape(f);
            for (s, e) in spans {
                let r = Rect::from_min_max(Pos2::new(x(s), fy - 5.0), Pos2::new(x(e).max(x(s) + 2.0), fy + 5.0));
                draw_track_feature(p, r, f.strand, sh.shape, fill, line, st);
                match sh.label {
                    1 => { p.text(r.center_top() - Vec2::new(0.0, 1.0), Align2::CENTER_BOTTOM, &f.label, small.clone(), fg); }
                    3 => { p.text(r.center_bottom() + Vec2::new(0.0, 1.0), Align2::CENTER_TOP, &f.label, small.clone(), fg); }
                    _ => { p.with_clip_rect(r).text(r.center(), Align2::CENTER_CENTER, &f.label, small.clone(), fg); }
                }
            }
        }
        let mut placed: Vec<Rect> = vec![];
        for s in &self.sites {
            let galley = p.layout_no_wrap(self.enz[s.enz].name.clone(), small.clone(), fg);
            let Some(lrect) = (0..10).map(|k| Align2::CENTER_BOTTOM.anchor_size(Pos2::new(x(s.pos), y - 14.0 - k as f32 * 11.0), galley.size()))
                .find(|lr| !placed.iter().any(|q| q.intersects(*lr))) else { continue };
            placed.push(lrect);
            p.line_segment([Pos2::new(x(s.pos), y), lrect.center_bottom()], Stroke::new(0.6, fg));
            p.galley(lrect.min, galley, fg);
        }
    }

    /// SeqBuilder's Sequence view: fixed column ruler, rows of top strand with a tick axis underneath,
    /// track markers (5', o) on the left, last position on the right, selection filling the whole row.
    fn sequence_view(&mut self, ui: &mut egui::Ui) {
        let Some(doc) = &self.doc else { return };
        let font = FontId::monospace(12.0);
        let cw = ui.painter().layout_no_wrap("a".into(), font.clone(), Color32::WHITE).size().x;
        let lh = 16.0;
        let (gutter, right_w) = (cw * 5.0, cw * 8.0);
        let gap = |i: usize, g: usize| if g == 0 { 0.0 } else { (i / g) as f32 * cw * 0.6 };
        let grouping = self.grouping;
        // bases per row: largest multiple of 10 that fits
        let avail = ui.available_width() - gutter - right_w - 16.0;
        let per_col = cw * if grouping == 0 { 1.0 } else { 1.0 + 0.6 / grouping as f32 };
        let bpr = (((avail / per_col) as usize) / 10 * 10).max(10);
        // Lines of a row, as the tool strip switches them: enzymes, 5' strand, ruler axis, 3' strand,
        // reading frames 1-6, then feature lanes (at least one empty track, as SeqBuilder shows).
        let mut next = 0usize;
        let mut take = |on: bool, n: usize| on.then(|| { next += n; next - n });
        let enz_line = take(!self.hide_enzymes, 1);
        let top_line = take(true, 1).unwrap_or(0);
        let axis_line = take(!self.hide_rulers, 1);
        let bottom_line = take(self.complement, 1);
        let frame_line = take(self.translations, 6);
        let feat_line = take(true, doc.lanes.max(1)).unwrap_or(0);
        let tracks = next;
        let row_h = lh * tracks as f32 + if self.spacing { 18.0 } else { 8.0 };
        let rows = doc.seq.len().div_ceil(bpr).max(1); // an empty sequence still shows a row to type into
        let v = ui.visuals();
        let (fg, dim, sel_bg) = (v.text_color(), v.weak_text_color(), v.selection.bg_fill);
        let frame = v.widgets.noninteractive.bg_stroke.color;

        // column ruler, fixed above the scrolling rows
        let (ruler, _) = ui.allocate_exact_size(Vec2::new(ui.available_width(), 18.0), Sense::hover());
        let x0 = ruler.left() + gutter;
        let col_x = |c: usize| x0 + c as f32 * cw + gap(c, grouping);
        {
            let p = ui.painter_at(ruler);
            p.line_segment([ruler.left_bottom(), ruler.right_bottom()], Stroke::new(1.0, fg));
            for c in 0..bpr {
                let x = col_x(c) + cw / 2.0;
                let h = if (c + 1) % 10 == 0 { 6.0 } else if (c + 1) % 5 == 0 { 4.0 } else { 2.0 };
                p.line_segment([Pos2::new(x, ruler.bottom() - h), Pos2::new(x, ruler.bottom())], Stroke::new(0.6, fg));
                if (c + 1) % 10 == 0 {
                    p.text(Pos2::new(x, ruler.bottom() - 6.0), Align2::CENTER_BOTTOM, (c + 1).to_string(), FontId::proportional(10.0), fg);
                }
            }
        }

        let mut area = egui::ScrollArea::vertical().auto_shrink(false);
        if let Some(pos) = self.scroll_to.take() { area = area.vertical_scroll_offset((pos / bpr) as f32 * row_h); }
        let (mut sel, mut anchor, caret) = (self.sel, self.anchor, self.caret);
        area.show_rows(ui, row_h, rows, |ui, range| {
            let origin = ui.cursor().min - Vec2::new(0.0, range.start as f32 * row_h);
            let base_at = |p: Pos2| -> usize {
                let row = ((p.y - origin.y) / row_h).max(0.0) as usize;
                let mut col = 0;
                while col < bpr && col_x(col) + cw / 2.0 < p.x { col += 1; }
                (row * bpr + col).min(doc.seq.len())
            };
            for row in range {
                let y = origin.y + row as f32 * row_h;
                let line_y = |l: usize| y + lh * l as f32;
                let ys = line_y(top_line);
                let (a, b) = (row * bpr, ((row + 1) * bpr).min(doc.seq.len()));
                let rect = Rect::from_min_size(Pos2::new(origin.x, y), Vec2::new(ui.available_width(), row_h));
                let resp = ui.allocate_rect(rect, Sense::click_and_drag());
                let p = ui.painter();
                let x = |i: usize| col_x(i - a);
                let band = Rect::from_min_max(Pos2::new(x0 - 2.0, ys - 1.0), Pos2::new(col_x(bpr - 1) + cw + 2.0, y + row_h - 6.0));
                if let Some((s, e)) = sel {
                    let (s, e) = (s.max(a), e.min(b));
                    if s < e { p.rect_filled(Rect::from_x_y_ranges(x(s)..=x(e - 1) + cw, band.y_range()), 0.0, sel_bg); }
                }
                p.rect_stroke(band, 2.0, Stroke::new(1.0, frame), egui::StrokeKind::Outside);
                // track markers
                let mx = origin.x + cw;
                let mark = |l: usize, t: &str, f: FontId, c: Color32| { p.text(Pos2::new(mx, line_y(l)), Align2::LEFT_TOP, t, f, c); };
                mark(top_line, "5'", font.clone(), fg);
                if let Some(l) = enz_line { mark(l, "o", FontId::proportional(9.0), dim); }
                if let Some(l) = axis_line { mark(l, "o", FontId::proportional(9.0), dim); }
                if let Some(l) = bottom_line { mark(l, "3'", font.clone(), fg); }
                if let Some(l) = frame_line { for k in 0..6 { mark(l + k, &(k + 1).to_string(), FontId::proportional(10.0), dim); } }
                mark(feat_line, "o", FontId::proportional(9.0), dim);
                // enzymes above the strand
                let mut next_free = 0.0;
                for s in self.sites.iter().filter(|s| enz_line.is_some() && (a..b).contains(&s.pos)) {
                    let px = x(s.pos).max(next_free);
                    let r = p.text(Pos2::new(px, y + 1.0), Align2::LEFT_TOP, &self.enz[s.enz].name, FontId::proportional(10.0), fg);
                    p.line_segment([Pos2::new(x(s.pos), r.bottom()), Pos2::new(x(s.pos), ys)], Stroke::new(0.6, fg));
                    next_free = r.right() + 6.0;
                }
                if let Some(c) = caret.filter(|&c| (a..b).contains(&c) || (c == b && b == doc.seq.len())) {
                    let cx = if c < b || b == a { x(c) } else { x(b - 1) + cw } - 1.0;
                    p.line_segment([Pos2::new(cx, ys - 1.0), Pos2::new(cx, ys + lh)], Stroke::new(2.0, fg));
                }
                for i in a..b {
                    p.text(Pos2::new(x(i), ys), Align2::LEFT_TOP, (doc.seq[i] as char).to_string(), font.clone(), fg);
                    if let Some(l) = bottom_line { p.text(Pos2::new(x(i), line_y(l)), Align2::LEFT_TOP, (complement(doc.seq[i]) as char).to_string(), font.clone(), dim); }
                }
                if let Some(l) = frame_line {
                    for i in a..b {
                        for k in 0..6 {
                            let Some(aa) = frame_residue(&doc.seq, i, k) else { continue };
                            let c = if aa == '*' { Color32::RED } else { fg };
                            p.text(Pos2::new(x(i), line_y(l + k)), Align2::LEFT_TOP, aa.to_string(), font.clone(), c);
                        }
                    }
                }
                // tick axis: every base, longer every 5 and 10
                if let Some(l) = axis_line {
                    let ay = line_y(l) + lh / 2.0;
                    if b > a { p.line_segment([Pos2::new(x(a), ay), Pos2::new(x(b - 1) + cw, ay)], Stroke::new(1.0, fg)); }
                    for i in a..b {
                        let h = if (i + 1) % 10 == 0 { 5.0 } else if (i + 1) % 5 == 0 { 3.5 } else { 2.0 };
                        let cx = x(i) + cw / 2.0;
                        p.line_segment([Pos2::new(cx, ay - h), Pos2::new(cx, ay + h)], Stroke::new(0.8, fg));
                    }
                }
                p.text(Pos2::new(band.right() + cw * 2.0, ys + lh / 2.0), Align2::LEFT_CENTER, b.to_string(), FontId::proportional(11.0), fg);
                for f in doc.feats.iter().filter(|f| f.show && covers(f, a, b)) {
                    let (s, e) = if f.start < f.end { (f.start.max(a), f.end.min(b)) } else if f.start < b { (f.start.max(a), b) } else { (a, f.end.min(b)) };
                    let fy = line_y(feat_line + f.lane) + 2.0;
                    let r = Rect::from_min_max(Pos2::new(x(s), fy), Pos2::new(x(e - 1) + cw, fy + lh - 4.0));
                    let (fill, line) = colors(f);
                    p.rect(r, 2.0, fill, Stroke::new(1.0, line), egui::StrokeKind::Inside);
                    p.with_clip_rect(r.shrink(1.0)).text(r.left_center() + Vec2::new(3.0, 0.0), Align2::LEFT_CENTER,
                        format!("{} {}", if f.strand < 0 { "◀" } else { "▶" }, f.label), FontId::proportional(10.0), Color32::BLACK);
                }
                if let Some(pp) = resp.interact_pointer_pos() {
                    let i = base_at(pp);
                    if resp.drag_started() || resp.clicked() {
                        anchor = if ui.input(|inp| inp.modifiers.shift) { sel.map_or(i, |s| s.0) } else { i };
                    }
                    sel = Some((anchor.min(i), anchor.max(i)));
                }
            }
        });
        if let Some((a, b)) = sel {
            // a click without dragging places the caret
            if a == b { (sel, self.caret) = (None, Some(a)); } else { self.caret = None; }
        }
        self.sel = sel;
        self.anchor = anchor;
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let ctx = &ui.ctx().clone();
        self.shortcuts(ctx);
        // Debug aid for visual checks: SEQBUILDER_SHOT=out.ppm saves the window and quits.
        self.frames += 1;
        if let (3, Ok(out)) = (self.frames, std::env::var("SEQBUILDER_SHOT")) {
            ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::ScreenshotCallback::new(move |img| {
                let mut ppm = format!("P6 {} {} 255\n", img.width(), img.height()).into_bytes();
                ppm.extend(img.pixels.iter().flat_map(|c| [c.r(), c.g(), c.b()]));
                let _ = std::fs::write(&out, ppm);
                std::process::exit(0);
            })));
        }
        if self.frames < 4 { ctx.request_repaint(); }
        if self.frames == 1 {
            if let Ok(v) = std::env::var("SEQBUILDER_VIEW") {
                self.view = match v.as_str() { "linear" => View::Linear, "sequence" => View::Sequence, _ => self.view };
            }
            // Debug aid: SEQBUILDER_FEATURE=n opens with that table row picked.
            if let Some(i) = std::env::var("SEQBUILDER_FEATURE").ok().and_then(|v| v.parse().ok()) { (self.sel_feat, self.side_top) = (Some(i), SideTop::Style); }
            if std::env::var("SEQBUILDER_NEW").is_ok() { self.new_dna(); }
            if let Ok(t) = std::env::var("SEQBUILDER_TAB") {
                self.tab = match t.as_str() { "minimap" => Tab::Minimap, "primers" => Tab::Primers, "comment" => Tab::Comment, "sites" => Tab::SiteSummary, _ => Tab::Features };
            }
            if let Ok(t) = std::env::var("SEQBUILDER_TRACKS") { (self.complement, self.translations) = (t.contains("strand"), t.contains("frames")); }
            if let Some(c) = std::env::var("SEQBUILDER_CARET").ok().and_then(|v| v.parse().ok()) { self.caret = Some(c); }
        }

        egui::Panel::top("menu").show(ui, |ui| self.menu(ui));
        egui::Panel::left("tools").exact_size(48.0).resizable(false).show(ui, |ui| self.toolbar(ui));
        if !self.hide_side { egui::Panel::right("side").resizable(true).default_size(300.0).show(ui, |ui| self.side_panel(ui)); }
        if !self.hide_table { egui::Panel::bottom("table").resizable(true).default_size(240.0).show(ui, |ui| self.feature_table(ui)); }
        egui::CentralPanel::default().show(ui, |ui| {
            self.info_bar(ui);
            ui.separator();
            let circular = self.doc.as_ref().is_some_and(|d| d.circular);
            egui::Panel::bottom("view_tabs").show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.add_enabled_ui(circular, |ui| ui.selectable_value(&mut self.view, View::Circular, "Circular"));
                    ui.selectable_value(&mut self.view, View::Linear, "Linear");
                    ui.selectable_value(&mut self.view, View::Sequence, "Sequence");
                    todo_tab(ui, "Primer Design");
                });
            });
            let Some(d) = &self.doc else {
                ui.centered_and_justified(|ui| ui.label("File ▸ Open… or drop a GenBank / FASTA file here"));
                return;
            };
            match self.view {
                View::Circular => self.page(ui, |p, r, ui| self.circular_map(p, r, ui, d)),
                View::Linear => self.page(ui, |p, r, ui| self.linear_map(p, r, ui, d)),
                View::Sequence => self.sequence_view(ui),
            }
        });
        self.selector_manager(ctx);
        if self.edit_open && self.sel_feat.is_some() {
            let mut open = true;
            egui::Window::new("Edit Feature").open(&mut open).collapsible(false).show(ctx, |ui| self.feature_editor(ui));
            self.edit_open = open;
        }
        if self.goto_open {
            let mut open = true;
            egui::Window::new("Go to Position").open(&mut open).collapsible(false).show(ctx, |ui| {
                let r = ui.text_edit_singleline(&mut self.goto_text);
                r.request_focus();
                if r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    if let (Ok(n), Some(d)) = (self.goto_text.trim().parse::<usize>(), &self.doc) {
                        let pos = n.clamp(1, d.seq.len().max(1)) - 1;
                        self.view = View::Sequence;
                        self.scroll_to = Some(pos);
                        self.sel = Some((pos, pos + 1));
                    }
                    self.goto_open = false;
                }
            });
            self.goto_open &= open;
        }
    }
}

fn main() -> eframe::Result {
    // seqbuilder-rs --export-genbank in.sbd out.gb  (batch conversion, no window)
    let args: Vec<String> = std::env::args().collect();
    if let [_, flag, input, output] = &args[..] {
        if flag == "--export-genbank" {
            let written = load(Path::new(input)).and_then(|d| std::fs::write(output, export_genbank(&d)).map_err(|e| e.to_string()));
            if let Err(e) = written { eprintln!("{input}: {e}"); std::process::exit(1); }
            return Ok(());
        }
    }
    let mut app = App { zoom: 1.0, ..Default::default() };
    let dir = std::env::var("SEQBUILDER_ENZYMES").map(PathBuf::from).unwrap_or_else(|_| DNASTAR_DIR.into());
    // User-edited selectors live where DataManager's DMConfig.cfg points; fall back to the shipped file.
    let home = std::env::var("HOME").unwrap_or_default();
    let user_sel = std::fs::read_to_string(format!("{home}/Library/Preferences/DNASTAR/DataManager/1.0/DMConfig.cfg")).ok()
        .and_then(|c| c.lines().find_map(|l| l.strip_prefix("selectorfile=").map(PathBuf::from)))
        .filter(|p| p.exists());
    app.load_enzymes(&dir.join("Enzymes.ezd"), &user_sel.unwrap_or_else(|| dir.join("EnzymeSelectors.sel")));
    if let Some(p) = std::env::args().nth(1) { app.open(Path::new(&p)); }
    let opts = eframe::NativeOptions { viewport: egui::ViewportBuilder::default().with_inner_size([1440.0, 900.0]), ..Default::default() };
    eframe::run_native("SeqBuilder RS", opts, Box::new(|_| Ok(Box::new(app))))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_and_gc() {
        assert_eq!(b"ACGTN".map(complement), *b"TGCAN");
        assert_eq!(gc_percent(b"GGAT"), 50.0);
    }

    #[test]
    fn loads_genbank_features() {
        let d = load(Path::new("testdata/pBR322.gb")).unwrap();
        assert_eq!(d.seq.len(), 4361);
        assert!(d.circular);
        let bla = d.feats.iter().find(|f| f.label == "bla").expect("bla gene");
        assert_eq!(bla.strand, -1);
        assert!(d.feats.iter().any(|f| f.kind == "source" && !f.show));
        // overlapping visible features never share a lane
        let vis: Vec<_> = d.feats.iter().filter(|f| f.show && f.start < f.end).collect();
        for (i, f) in vis.iter().enumerate() {
            for g in &vis[i + 1..] {
                if f.lane == g.lane { assert!(f.end <= g.start || g.end <= f.start, "{} / {}", f.label, g.label); }
            }
        }
    }

    #[test]
    fn site_search_wraps_origin() {
        // EcoRI GAATTC split across the origin of a circular sequence
        assert_eq!(find_sites(b"TTCAAAAGAA", true, "GAATTC", false), vec![7]);
        assert_eq!(find_sites(b"TTCAAAAGAA", false, "GAATTC", false), Vec::<usize>::new());
        assert_eq!(find_sites(b"AACAGTTTCTGAA", false, "CAGNNNCTG", false), vec![2]);
        // an N in the sequence is only a potential cut
        assert_eq!(find_sites(b"TTGANTTCTT", false, "GAATTC", false), Vec::<usize>::new());
        assert_eq!(find_sites(b"TTGANTTCTT", false, "GAATTC", true), vec![2]);
    }

    #[test]
    fn sbd_matches_genbank() {
        let (s, g) = (load(Path::new("testdata/pBR322.sbd")).unwrap(), load(Path::new("testdata/pBR322.gb")).unwrap());
        assert_eq!(s.seq.to_ascii_uppercase(), g.seq.to_ascii_uppercase());
        assert_eq!((s.name.as_str(), s.circular, s.feats.len()), ("SYNPBR322", true, g.feats.len()));
        for (a, b) in s.feats.iter().zip(&g.feats) {
            assert_eq!((&a.kind, a.start, a.end, a.strand), (&b.kind, b.start, b.end, b.strand), "{}", b.range);
        }
        let red = load(Path::new("testdata/pBR322_blaRED.sbd")).unwrap();
        assert_eq!(red.feats[38].style.unwrap().fill, Some([255, 0, 0]));
    }

    #[test]
    fn saves_style_edits_to_sbd() {
        let mut app = App::default();
        app.open(Path::new("testdata/pBR322.sbd"));
        let f = &mut app.doc.as_mut().unwrap().feats[38];
        f.style = Some(sbd::Style { line: Some([0, 0, 255]), weight: 1.0, fill: Some([0, 200, 0]), shadow: None });
        f.edited = true;
        let out = std::env::temp_dir().join("seqbuilder-rs-test.sbd");
        app.save(Some(out.clone()));
        assert!(app.error.is_none(), "{:?}", app.error);
        let back = sbd::read(&std::fs::read(&out).unwrap()).unwrap();
        assert_eq!(back.obj_styles[&38].fill, Some([0, 200, 0]));
        assert_eq!(back.obj_styles[&38].weight, 1.0);
        assert_eq!(back.feats.len(), 50);
    }

    #[test]
    fn saves_feature_edits_to_sbd() {
        let mut app = App::default();
        app.open(Path::new("testdata/pBR322.sbd"));
        let f = &mut app.doc.as_mut().unwrap().feats[26]; // ROP protein CDS 1915..2106
        assert_eq!(f.name, "ROP protein");
        (f.name, f.start, f.end, f.strand) = ("ROP".into(), 1914, 2100, -1);
        set_qual(&mut f.quals, "note", "edited in SeqBuilder RS");
        f.dirty = true;
        let out = std::env::temp_dir().join("seqbuilder-rs-feature-test.sbd");
        app.save(Some(out.clone()));
        assert!(app.error.is_none(), "{:?}", app.error);
        let back = sbd::read(&std::fs::read(&out).unwrap()).unwrap();
        let rop = &back.feats[26];
        assert_eq!((rop.name.as_str(), rop.start, rop.end, rop.complement), ("ROP", 2100, 1915, true));
        assert_eq!(rop.qual("note").as_deref(), Some("edited in SeqBuilder RS"));
        assert_eq!(rop.qual("product").as_deref(), Some("ROP protein"));
        let orig = sbd::read(&std::fs::read("testdata/pBR322.sbd").unwrap()).unwrap();
        for (a, b) in back.feats.iter().zip(&orig.feats).filter(|(a, _)| a.id != 26) {
            assert_eq!((&a.name, a.start, a.end, &a.quals), (&b.name, b.start, b.end, &b.quals));
        }
    }

    #[test]
    fn adds_deletes_and_hides_features_in_sbd() {
        let mut app = App::default();
        app.open(Path::new("testdata/pBR322.sbd"));
        assert!(!app.doc.as_ref().unwrap().feats[0].show, "source starts hidden, as stored in the file");
        app.sel_feat = Some(2); // echinomycin 24..27
        app.delete_feature();
        app.sel = Some((99, 200));
        app.new_feature();
        app.doc.as_mut().unwrap().feats[5].show = false;
        let out = std::env::temp_dir().join("seqbuilder-rs-add-del-test.sbd");
        app.save(Some(out.clone()));
        assert!(app.error.is_none(), "{:?}", app.error);
        let back = sbd::read(&std::fs::read(&out).unwrap()).unwrap();
        assert_eq!(back.feats.len(), 50);
        assert!(back.feats.iter().all(|f| f.id != 2));
        let new = back.feats.iter().find(|f| f.name == "New Feature").unwrap();
        assert_eq!((new.id, new.kind.as_str(), new.start, new.end), (50, "misc_feature", 100, 200));
        let hidden = app.doc.as_ref().unwrap().feats[5].sbd_id.unwrap();
        assert_eq!((back.shown[&hidden], back.shown[&50]), (false, true));
    }

    /// Inserting GGGGG before base 10 and deleting 2000..2009 must give the sequence and features
    /// SeqBuilder Pro saved in pBR322_ins.sbd and pBR322_cut.sbd.
    #[test]
    fn sequence_edits_match_seqbuilder() {
        for (sample, edit) in [("pBR322_ins", (9, 9, &b"GGGGG"[..])), ("pBR322_cut", (1999, 2009, &b""[..]))] {
            let mut app = App::default();
            app.open(Path::new("testdata/pBR322.sbd"));
            app.edit_seq(edit.0, edit.1, edit.2);
            let out = std::env::temp_dir().join(format!("seqbuilder-rs-{sample}.sbd"));
            app.save(Some(out.clone()));
            assert!(app.error.is_none(), "{:?}", app.error);
            let ours = std::fs::read(&out).unwrap();
            let want = std::fs::read(format!("testdata/{sample}.sbd")).unwrap();
            let (got, exp) = (sbd::read(&ours).unwrap(), sbd::read(&want).unwrap());
            assert!(got.seq == exp.seq, "{sample}: sequence");
            assert_eq!(got.feats.len(), exp.feats.len());
            for (g, e) in got.feats.iter().zip(&exp.feats) {
                assert_eq!((g.id, g.start, g.end, &g.quals), (e.id, e.start, e.end, &e.quals), "{sample}: {}", e.name);
            }
            assert!(sbd::uvsi(&ours) == sbd::uvsi(&want), "{sample}: UVSI");
        }
    }

    #[test]
    fn undo_restores_sequence_and_features() {
        let mut app = App::default();
        app.open(Path::new("testdata/pBR322.sbd"));
        let before = (app.doc.as_ref().unwrap().seq.clone(), app.doc.as_ref().unwrap().feats[2].start);
        app.edit_seq(10, 30, b""); // swallows echinomycin 24..27
        assert_eq!(app.doc.as_ref().unwrap().feats.len(), 49);
        app.undo_redo(false);
        let d = app.doc.as_ref().unwrap();
        assert_eq!((d.seq.clone(), d.feats[2].start, d.feats.len()), (before.0, before.1, 50));
        assert!(d.deleted.is_empty());
    }

    /// Saving pBR322.gb as .sbd (from the template) must give what SeqBuilder Pro made of it: pBR322.sbd.
    #[test]
    fn genbank_saved_as_sbd_matches_seqbuilder() {
        let mut app = App::default();
        app.open(Path::new("testdata/pBR322.gb"));
        let out = std::env::temp_dir().join("pBR322.sbd");
        app.save(Some(out.clone()));
        assert!(app.error.is_none(), "{:?}", app.error);
        let ours = std::fs::read(&out).unwrap();
        let want = std::fs::read("testdata/pBR322.sbd").unwrap();
        let (got, exp) = (sbd::read(&ours).unwrap(), sbd::read(&want).unwrap());
        assert!(got.seq == exp.seq);
        assert_eq!((got.circular, got.header == exp.header), (true, true));
        assert_eq!(sbd::sb_linear(&ours), sbd::sb_linear(&want));
        assert!(sbd::uvsi(&ours) == sbd::uvsi(&want));
        assert_eq!(got.feats.len(), exp.feats.len());
        for (g, e) in got.feats.iter().zip(&exp.feats) {
            assert_eq!((g.id, &g.name, &g.kind, g.start, g.end, g.complement), (e.id, &e.name, &e.kind, e.start, e.end, e.complement));
            assert_eq!(g.quals, e.quals, "{}", e.name);
        }
        assert_eq!(got.shown, exp.shown);
        assert!(got.obj_styles.is_empty());
    }

    /// Exporting pBR322_add.sbd must give what SeqBuilder Pro's File > Export wrote: pBR322_add.gbk.
    #[test]
    fn genbank_export_matches_seqbuilder() {
        let want = std::fs::read_to_string("testdata/pBR322_add.gbk").unwrap();
        let out = export_genbank(&load(Path::new("testdata/pBR322_add.sbd")).unwrap());
        if out != want {
            let line = out.lines().zip(want.lines()).position(|(a, b)| a != b);
            panic!("differs at line {line:?}: {:?}", line.map(|i| (out.lines().nth(i), want.lines().nth(i))));
        }
    }

    /// The same data as NCBI's pBR322.gb whether exported from GenBank or from SeqBuilder's .sbd.
    #[test]
    fn genbank_export_keeps_ncbi_content() {
        let ncbi = gb_io::reader::parse_file("testdata/pBR322.gb").unwrap().remove(0);
        for f in ["testdata/pBR322.gb", "testdata/pBR322.sbd"] {
            let back = gb_io::reader::parse_slice(export_genbank(&load(Path::new(f)).unwrap()).as_bytes()).unwrap().remove(0);
            assert!(back.seq == ncbi.seq && back.is_circular() && back.definition == ncbi.definition, "{f}");
            assert_eq!(back.features.len(), ncbi.features.len());
            let key = |f: &gb_io::seq::Feature| (f.location.find_bounds().unwrap(), f.kind.to_string());
            let (mut a, mut b): (Vec<_>, Vec<_>) = (back.features.iter().map(key).collect(), ncbi.features.iter().map(key).collect());
            (a.sort(), b.sort());
            assert_eq!(a, b, "{f}");
        }
        let out = export_genbank(&load(Path::new("testdata/ALDH2.sbd")).unwrap());
        let back = gb_io::reader::parse_slice(out.as_bytes()).unwrap().remove(0);
        assert_eq!((back.seq.len(), back.features.len(), back.is_circular()), (114842, 2, false));
        assert_eq!(back.features[1].location.to_gb_format(), "complement(5087..5124)");
    }

    /// New DNA, typed bases and a feature saved as .sbd; then a first feature added to a document with none.
    #[test]
    fn new_dna_saves_as_sbd() {
        let mut app = App::default();
        app.new_dna();
        let caret = app.caret.unwrap();
        app.edit_seq(caret, caret, b"atgaaacccgggtttTAA");
        app.sel = Some((3, 12));
        app.new_feature();
        let out = std::env::temp_dir().join("seqbuilder-rs-new.sbd");
        app.save(Some(out.clone()));
        assert!(app.error.is_none(), "{:?}", app.error);
        let back = sbd::read(&std::fs::read(&out).unwrap()).unwrap();
        assert_eq!((back.seq.as_slice(), back.circular, back.feats.len()), (&b"atgaaacccgggtttTAA"[..], false, 1));
        assert_eq!((back.feats[0].start, back.feats[0].end, back.header_field("LOCUS", " ")), (4, 12, String::new()));

        let mut empty = App::default();
        empty.new_dna();
        empty.edit_seq(0, 0, b"acgtacgt");
        let out2 = std::env::temp_dir().join("seqbuilder-rs-new-empty.sbd");
        empty.save(Some(out2.clone()));
        empty.open(&out2);
        empty.sel = Some((0, 4));
        empty.new_feature();
        empty.save(None);
        let back = sbd::read(&std::fs::read(&out2).unwrap()).unwrap();
        assert_eq!((back.feats.len(), back.feats[0].end, back.shown[&back.feats[0].id]), (1, 4, true));
    }

    #[test]
    fn reading_frames() {
        let seq = b"ATGAAATAG"; // M K * in frame 1
        let f1: String = (0..seq.len()).filter_map(|i| frame_residue(seq, i, 0)).collect();
        assert_eq!(f1, "MK*");
        // reverse strand CTATTTCAT read from the last base: L F H
        let r1: String = (0..seq.len()).rev().filter_map(|i| frame_residue(seq, i, 3)).collect();
        assert_eq!(r1, "LFH");
    }

    /// Statistics SeqBuilder Pro's Details panel showed for a 37-base feature of an edited ALDH2
    /// (its 1240..1276 is 1253..1289 of the saved ALDH2.sbd).
    #[test]
    fn base_statistics_match_seqbuilder() {
        let d = load(Path::new("testdata/ALDH2.sbd")).unwrap();
        let rows: std::collections::HashMap<String, String> = base_stats(&d.seq[1252..1289]).into_iter().collect();
        for (k, v) in [("% A", "21.62 [8]"), ("% A+T", "56.76 [21]"), ("% Ambiguous", "0.00 [0]"), ("% C", "32.43 [12]"),
            ("% G", "10.81 [4]"), ("% G+C", "43.24 [16]"), ("% T", "35.14 [13]"), ("% diCG", "0.00 [0]"),
            ("Total number of bases", "37"), ("% AA", "2.78 [1]"), ("% AC", "8.33 [3]"), ("% AG", "5.56 [2]"),
            ("% AT", "5.56 [2]"), ("% CA", "11.11 [4]"), ("% CC", "8.33 [3]")] {
            assert_eq!(rows[k], v, "{k}");
        }
    }

    #[test]
    fn primers_tab_values() {
        let d = load(Path::new("testdata/ALDH2.sbd")).unwrap();
        let p = primer_sequence(&d.feats[0], &d.seq);
        assert_eq!(String::from_utf8_lossy(&p).to_uppercase(), "AGGCTGAGGCAGGAGAATGGCATGAACCCGGGAGGCGG");
        assert_eq!(format!("{:.1}", primer_tm(&p)), "77.7"); // SeqBuilder's Primers tab
        let r = primer_sequence(&d.feats[1], &d.seq);
        assert_eq!((r.len(), format!("{:.1}", primer_tm(&r))), (38, "77.7".to_string()));
    }

    /// Matches SeqBuilder Pro 19 showing pBR322 with its default "Unique Six+ cutters" selector.
    #[test]
    fn dnastar_unique_six_cutters_on_pbr322() {
        let dir = Path::new(DNASTAR_DIR);
        let (Ok(ezd), Ok(sel)) = (std::fs::read(dir.join("Enzymes.ezd")), std::fs::read_to_string(dir.join("EnzymeSelectors.sel"))) else {
            eprintln!("DNASTAR not installed, skipping");
            return;
        };
        let mut app = App::default();
        app.enz = enzymes::parse_ezd(&ezd);
        app.sels = enzymes::parse_sel(&sel).unwrap();
        assert_eq!((app.enz.len(), app.sels.len()), (264, 25));
        app.open(Path::new("testdata/pBR322.gb"));
        assert_eq!(app.shown.len(), app.sites.len(), "only cutting enzymes are checked");
        let mut got: Vec<&str> = app.sites.iter().map(|s| app.enz[s.enz].name.as_str()).collect();
        got.sort();
        let mut want = vec![
            "SspI", "ZraI", "AatII", "EcoRI", "ApoI", "ClaI", "HindIII", "EcoRV", "NheI", "BmtI", "BamHI", "SgrAI",
            "SphI", "EcoNI", "SalI", "PshAI", "EagI", "NruI", "BspMI", "BsmI", "StyI", "AvaI", "MscI", "Bpu10I",
            "BsgI", "BspEI", "BsaBI", "PvuII", "PfoI", "BsmBI", "Tth111I", "BsaAI", "NdeI", "BstZ17I", "SapI",
            "AflIII", "PciI", "AlwNI", "AhdI", "BsaI", "AseI", "PstI", "PvuI", "ScaI",
        ];
        want.sort();
        assert_eq!(got, want);
    }
}
