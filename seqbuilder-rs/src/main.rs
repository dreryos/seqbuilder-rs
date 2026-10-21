// SeqBuilder-like DNA viewer laid out like SeqBuilder Pro: tool strip, map page with view tabs,
// feature table with tabs, Style/Enzymes and Settings/Details side panels.
mod enzymes;
mod sbd;

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, Vec2};
use enzymes::{Enzyme, Selector};
use std::collections::HashSet;
use std::f32::consts::{PI, TAU};
use std::path::{Path, PathBuf};

const DNASTAR_DIR: &str = "/Applications/DNASTAR/Lasergene19";

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
    style: Option<(Color32, Color32)>, // (fill, line) of this feature, from the document or the Style tab
    base: Option<(Color32, Color32)>,  // the document's style for this feature type
    sbd_id: Option<u32>,               // UVFT id when read from a .sbd
    edited: bool,                      // style changed since the last save
    quals: Vec<String>,                // "key=value" qualifier lines, quoting included
    dirty: bool,                       // name, type, range or qualifiers changed since the last save
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
        Doc { name, definition: String::new(), accession: String::new(), comment: String::new(), seq, circular: false, feats: vec![], lanes: 0, hits: vec![], path: path.into(), sbd: None, deleted: vec![] }
    } else if path.extension().is_some_and(|e| e.eq_ignore_ascii_case("sbd")) {
        let s = sbd::read(&text)?;
        let rgb = |c: [u8; 3]| Color32::from_rgb(c[0], c[1], c[2]);
        let feats = s.feats.iter().map(|f| {
            let (lo, hi) = (f.start.min(f.end), f.start.max(f.end));
            let q = |k: &str| f.qual(k);
            let colors = |st: &sbd::Style| {
                let fill = st.fill.map_or(Color32::TRANSPARENT, rgb);
                (fill, st.line.map_or(fill, rgb))
            };
            Feat {
                kind: f.kind.clone(), name: f.name.clone(), label: short(&f.name), desc: q("note").or_else(|| q("product")).unwrap_or_default(),
                range: if f.complement { format!("complement({lo}..{hi})") } else { format!("{lo}..{hi}") },
                start: lo.saturating_sub(1), end: hi, strand: if f.complement { -1 } else { 1 }, lane: 0,
                show: s.shown.get(&f.id).copied().unwrap_or(f.kind != "source"), style: s.obj_styles.get(&f.id).map(colors),
                base: s.type_styles.get(&f.kind).map(colors), sbd_id: Some(f.id), edited: false,
                quals: f.quals.clone(), dirty: false,
            }
        }).collect();
        let locus = s.header_field("LOCUS", " ");
        Doc {
            name: locus.split_whitespace().next().map_or(name, str::to_owned),
            definition: s.header_field("DEFINITION", " "), accession: s.header_field("ACCESSION", " "),
            comment: s.header_field("COMMENT", "\n"),
            seq: s.seq.clone(), circular: s.circular, feats, lanes: 0, hits: vec![], path: path.into(), sbd: Some(s), deleted: vec![],
        }
    } else {
        let rec = gb_io::reader::parse_slice(&text).map_err(|e| e.to_string())?
            .into_iter().next().ok_or("no sequence in file")?;
        let len = rec.seq.len();
        let feats = rec.features.iter().filter_map(|f| {
            let (s, e) = f.location.find_bounds().ok()?;
            let q = |k: &str| f.qualifier_values(k).next().map(str::to_owned);
            let name = ["label", "gene", "product", "note"].iter().find_map(|k| q(k)).unwrap_or_else(|| f.kind.to_string());
            let quals = f.qualifiers.iter().map(|(k, v)| match v {
                // GenBank leaves numbers and [citations] unquoted
                Some(v) if v.parse::<i64>().is_ok() || v.starts_with('[') => format!("{k}={v}"),
                Some(v) => format!("{k}=\"{v}\""),
                None => k.to_string(),
            }).collect();
            // find_bounds unwraps origin-spanning joins past len; fold them back.
            let (start, end) = (s.rem_euclid(len as i64) as usize, e as usize);
            let end = if end > len { end - len } else { end };
            Some(Feat {
                kind: f.kind.to_string(), label: short(&name), name, desc: q("note").or_else(|| q("product")).unwrap_or_default(),
                range: f.location.to_gb_format(), start, end, strand: location_strand(&f.location), lane: 0,
                show: f.kind != "source", // SeqBuilder lists source but leaves it hidden
                style: None, base: None, sbd_id: None, edited: false, quals, dirty: false,
            })
        }).collect();
        let circular = rec.is_circular();
        Doc {
            name: rec.name.clone().unwrap_or(name), definition: rec.definition.clone().unwrap_or_default(),
            accession: rec.accession.clone().unwrap_or_default(), comment: rec.comments.join("\n"),
            seq: rec.seq, circular, feats, lanes: 0, hits: vec![], path: path.into(), sbd: None, deleted: vec![],
        }
    };
    doc.lanes = assign_lanes(&mut doc.feats, doc.seq.len());
    Ok(doc)
}

fn kind_color(kind: &str) -> (Color32, Color32) {
    let rgb = Color32::from_rgb;
    match kind {
        "CDS" | "gene" | "mat_peptide" => (rgb(215, 228, 255), rgb(70, 110, 220)),
        "promoter" => (rgb(170, 225, 160), rgb(40, 130, 40)),
        "rep_origin" => (rgb(240, 235, 120), rgb(150, 140, 20)),
        "terminator" => (rgb(240, 170, 170), rgb(170, 40, 40)),
        "primer_bind" => (rgb(170, 210, 245), rgb(30, 110, 190)),
        "misc_binding" | "protein_bind" => (rgb(215, 190, 245), rgb(130, 60, 200)),
        "misc_feature" => (rgb(230, 200, 160), rgb(130, 75, 20)),
        _ => (rgb(215, 215, 220), rgb(110, 110, 125)),
    }
}

fn colors(f: &Feat) -> (Color32, Color32) {
    f.style.or(f.base).unwrap_or_else(|| kind_color(&f.kind))
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
enum Tab { #[default] Features, Comment, SiteSummary }

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
    zoom: f32,
    view: View,
    tab: Tab,
    side_top: SideTop,
    side_bottom: SideBottom,
    goto_open: bool,
    goto_text: String,
    scroll_to: Option<usize>,
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
            f.range = if comp { format!("complement({lo}..{hi})") } else { format!("{lo}..{hi}") };
            f.label = short(&f.name);
            if f.quals.iter().any(|q| q.starts_with("dnas_title=")) { set_qual(&mut f.quals, "dnas_title", &f.name.clone()); }
            set_qual(&mut f.quals, "note", &note);
            f.desc = if note.is_empty() { f.desc.clone() } else { note };
            f.dirty = true;
            d.lanes = assign_lanes(&mut d.feats, len);
        }
    }

    /// Features ▸ New Feature (⌘K): a misc_feature over the selection, picked for editing in Details.
    fn new_feature(&mut self) {
        let (Some(d), Some((a, b))) = (&mut self.doc, self.sel) else { return };
        let name = "New Feature".to_string();
        d.feats.push(Feat {
            kind: "misc_feature".into(), label: name.clone(), desc: String::new(), range: format!("{}..{b}", a + 1),
            start: a, end: b, strand: 1, lane: 0, show: true, style: None, base: None, sbd_id: None, edited: false,
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
        let Some(s) = &mut d.sbd else {
            self.error = Some("Saving is only supported for documents opened from .sbd so far".into());
            return;
        };
        for id in d.deleted.drain(..) { s.remove_feature(id); }
        for f in d.feats.iter_mut().filter(|f| f.dirty) {
            let (lo, hi) = (f.start + 1, f.end);
            let comp = f.strand < 0;
            let (start, end) = if comp { (hi, lo) } else { (lo, hi) };
            let feat = sbd::SbdFeat { id: f.sbd_id.unwrap_or(0), name: f.name.clone(), kind: f.kind.clone(), start, end, complement: comp, quals: f.quals.clone() };
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
        let rgb = |c: Color32| [c.r(), c.g(), c.b()];
        for f in d.feats.iter_mut().filter(|f| f.edited) {
            let Some(id) = f.sbd_id else { continue };
            let weight = s.obj_styles.get(&id).or(s.type_styles.get(&f.kind)).map_or(1.0, |st| st.weight);
            s.set_feature_style(id, f.style.map(|(fill, line)| sbd::Style { fill: (fill.a() > 0).then(|| rgb(fill)), line: Some(rgb(line)), weight }));
            f.edited = false;
        }
        let path = to.unwrap_or_else(|| d.path.clone());
        match std::fs::write(&path, s.to_bytes()) {
            Ok(()) => { d.path = path; self.error = None; }
            Err(e) => self.error = Some(format!("{}: {e}", path.display())),
        }
    }

    fn save_as(&mut self) {
        let name = self.doc.as_ref().and_then(|d| d.path.file_name()).map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
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
                if ui.button("Open…  ⌘O").clicked() {
                    if let Some(p) = rfd::FileDialog::new().add_filter("Sequence", &["sbd", "gb", "gbk", "genbank", "fasta", "fa", "fas", "seq"]).pick_file() {
                        self.open(&p);
                    }
                }
                if ui.button("Save  ⌘S").clicked() { self.save(None); }
                if ui.button("Save As…").clicked() { self.save_as(); }
                if ui.button("Open Enzyme Files…").clicked() {
                    let ezd = rfd::FileDialog::new().add_filter("Enzyme file", &["ezd"]).pick_file();
                    let sel = ezd.as_ref().and_then(|_| rfd::FileDialog::new().add_filter("Selector file", &["sel"]).pick_file());
                    if let (Some(e), Some(s)) = (ezd, sel) { self.load_enzymes(&e, &s); }
                }
            });
            ui.menu_button("Edit", |ui| {
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
                if ui.add_enabled(self.sel_feat.is_some(), egui::Button::new("Delete Feature  ⌫")).clicked() { self.delete_feature(); }
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
        if !ctx.egui_wants_keyboard_input() && ctx.input(|i| i.key_pressed(egui::Key::Backspace) || i.key_pressed(egui::Key::Delete)) {
            self.delete_feature();
        }
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

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.add(egui::Button::selectable(true, "↖")).on_hover_text("Select");
            ui.separator();
            if ui.button("+").on_hover_text("Zoom in").clicked() { self.zoom = (self.zoom * 1.25).min(8.0); }
            if ui.button("−").on_hover_text("Zoom out").clicked() { self.zoom = (self.zoom / 1.25).max(0.25); }
            if ui.button("⤢").on_hover_text("Fit to view").clicked() { self.zoom = 1.0; }
        });
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
                todo_tab(ui, "Primers");
                ui.selectable_value(&mut self.tab, Tab::Comment, "Comment");
                todo_tab(ui, "Minimap");
                ui.selectable_value(&mut self.tab, Tab::SiteSummary, "Site Summary");
                todo_tab(ui, "History");
            });
        });
        let Some(doc) = &mut self.doc else { return };
        let mut pick = None;
        egui::ScrollArea::both().auto_shrink(false).show(ui, |ui| match self.tab {
            Tab::Features => {
                let mut toggled = false;
                egui::Grid::new("feats").striped(true).num_columns(7).spacing([16.0, 4.0]).show(ui, |ui| {
                    for h in ["Show", "Type", "Name", "Range", "Strand", "Length", "Description"] { ui.strong(h); }
                    ui.end_row();
                    for (i, f) in doc.feats.iter_mut().enumerate() {
                        toggled |= ui.checkbox(&mut f.show, "").changed();
                        ui.label(&f.kind);
                        if ui.selectable_label(self.sel_feat == Some(i), &f.label).clicked() { pick = Some(i); }
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
                SideBottom::Details if self.sel_feat.is_some() && self.doc.is_some() => self.feature_editor(ui),
                SideBottom::Details => match &self.doc {
                    None => { ui.label("No document open."); }
                    Some(_) if self.sel.is_some() => {
                        let (a, b) = self.sel.unwrap();
                        ui.heading("Region selection");
                        ui.separator();
                        ui.label(format!("Selection: Forward {}..{}, length = {}", a + 1, b, b - a));
                        if let Some(s) = self.selected() { ui.label(format!("GC content: {:.1}%", gc_percent(s))); }
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
                ui.strong(&f.label);
                ui.label(&f.kind);
                let (mut fill, mut line) = colors(f);
                let mut changed = false;
                egui::Grid::new("style").num_columns(2).show(ui, |ui| {
                    ui.label("Fill:"); changed |= ui.color_edit_button_srgba(&mut fill).changed(); ui.end_row();
                    ui.label("Line:"); changed |= ui.color_edit_button_srgba(&mut line).changed(); ui.end_row();
                });
                if changed { f.style = Some((fill, line)); f.edited = true; }
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
            let head = (8.0 / rad).min(a1 - a0);
            let (b0, b1, tip, back) = if f.strand < 0 { (a0 + head, a1, a0, a0 + head) } else { (a0, a1 - head, a1, a1 - head) };
            if b1 > b0 {
                let mut band = arc(b0, b1, rad + 4.0);
                band.extend(arc(b1, b0, rad - 4.0));
                band.push(band[0]);
                p.add(Shape::line(arc(b0, b1, rad), Stroke::new(8.0, fill)));
                p.add(Shape::line(band, Stroke::new(1.0, line)));
            }
            p.add(Shape::convex_polygon(vec![at(tip, rad), at(back, rad + 7.0), at(back, rad - 7.0)], fill, Stroke::new(1.0, line)));
            let mid = (a0 + a1) / 2.0;
            let galley = p.layout_no_wrap(f.label.clone(), small.clone(), fg);
            let lrect = Align2::CENTER_CENTER.anchor_size(at(mid, rad - 13.0), galley.size());
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
            let fy = y + 28.0 + f.lane as f32 * 16.0;
            let spans = if f.start < f.end { vec![(f.start, f.end)] } else { vec![(f.start, doc.seq.len()), (0, f.end)] };
            for (s, e) in spans {
                let r = Rect::from_min_max(Pos2::new(x(s), fy - 5.0), Pos2::new(x(e).max(x(s) + 2.0), fy + 5.0));
                p.rect(r, 1.0, fill, Stroke::new(1.0, line), egui::StrokeKind::Inside);
                p.with_clip_rect(r).text(r.center(), Align2::CENTER_CENTER, &f.label, small.clone(), fg);
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
        let strands = if self.complement { 2 } else { 1 };
        // lines per row: enzymes, strand(s), axis, feature lanes (at least one empty track, as SeqBuilder shows)
        let tracks = 1 + strands + 1 + doc.lanes.max(1);
        let row_h = lh * tracks as f32 + 8.0;
        let rows = doc.seq.len().div_ceil(bpr);
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
        let (mut sel, mut anchor, show_comp) = (self.sel, self.anchor, self.complement);
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
                let ys = y + lh; // top strand line
                let ya = ys + lh * strands as f32; // axis line
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
                p.text(Pos2::new(mx, ys), Align2::LEFT_TOP, "5'", font.clone(), fg);
                for t in 1..tracks - 1 { p.text(Pos2::new(mx, ys + lh * t as f32), Align2::LEFT_TOP, "o", FontId::proportional(9.0), dim); }
                // enzymes above the strand
                let mut next_free = 0.0;
                for s in self.sites.iter().filter(|s| (a..b).contains(&s.pos)) {
                    let px = x(s.pos).max(next_free);
                    let r = p.text(Pos2::new(px, y + 1.0), Align2::LEFT_TOP, &self.enz[s.enz].name, FontId::proportional(10.0), fg);
                    p.line_segment([Pos2::new(x(s.pos), r.bottom()), Pos2::new(x(s.pos), ys)], Stroke::new(0.6, fg));
                    next_free = r.right() + 6.0;
                }
                for i in a..b {
                    p.text(Pos2::new(x(i), ys), Align2::LEFT_TOP, (doc.seq[i] as char).to_string(), font.clone(), fg);
                    if show_comp { p.text(Pos2::new(x(i), ys + lh), Align2::LEFT_TOP, (complement(doc.seq[i]) as char).to_string(), font.clone(), dim); }
                }
                // tick axis: every base, longer every 5 and 10
                let ay = ya + lh / 2.0;
                p.line_segment([Pos2::new(x(a), ay), Pos2::new(x(b - 1) + cw, ay)], Stroke::new(1.0, fg));
                for i in a..b {
                    let h = if (i + 1) % 10 == 0 { 5.0 } else if (i + 1) % 5 == 0 { 3.5 } else { 2.0 };
                    let cx = x(i) + cw / 2.0;
                    p.line_segment([Pos2::new(cx, ay - h), Pos2::new(cx, ay + h)], Stroke::new(0.8, fg));
                }
                p.text(Pos2::new(band.right() + cw * 2.0, ys + lh / 2.0), Align2::LEFT_CENTER, b.to_string(), FontId::proportional(11.0), fg);
                for f in doc.feats.iter().filter(|f| f.show && covers(f, a, b)) {
                    let (s, e) = if f.start < f.end { (f.start.max(a), f.end.min(b)) } else if f.start < b { (f.start.max(a), b) } else { (a, f.end.min(b)) };
                    let fy = ya + lh * (1 + f.lane) as f32 + 2.0;
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
        if sel.is_some_and(|(a, b)| a == b) { sel = None; }
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
            if let Some(i) = std::env::var("SEQBUILDER_FEATURE").ok().and_then(|v| v.parse().ok()) { self.sel_feat = Some(i); }
        }

        egui::Panel::top("menu").show(ui, |ui| self.menu(ui));
        egui::Panel::left("tools").exact_size(36.0).resizable(false).show(ui, |ui| self.toolbar(ui));
        egui::Panel::right("side").resizable(true).default_size(300.0).show(ui, |ui| self.side_panel(ui));
        egui::Panel::bottom("table").resizable(true).default_size(240.0).show(ui, |ui| self.feature_table(ui));
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
        assert_eq!(red.feats[38].style.unwrap().0, Color32::from_rgb(255, 0, 0));
    }

    #[test]
    fn saves_style_edits_to_sbd() {
        let mut app = App::default();
        app.open(Path::new("testdata/pBR322.sbd"));
        let f = &mut app.doc.as_mut().unwrap().feats[38];
        f.style = Some((Color32::from_rgb(0, 200, 0), Color32::from_rgb(0, 0, 255)));
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
