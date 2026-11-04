//! Agarose gel simulation (File ▸ New ▸ New Agarose Gel in SeqBuilder Pro).
//!
//! Hardware data is bundled from the original app's resources:
//! `Markers.tdl` (molecular rulers / ladders) and `Gelbox.tdl` (gel box
//! dimensions and well counts), both plain TSV.

/// Ladders from the original `Markers.tdl` (124 molecular rulers).
pub const MARKERS_TDL: &str = include_str!("../assets/Markers.tdl");
/// Gel boxes from the original `Gelbox.tdl`.
pub const GELBOX_TDL: &str = include_str!("../assets/Gelbox.tdl");

/// A molecular ruler from `Markers.tdl`: brand, product line, name and band
/// sizes in bp.
#[derive(Clone, Debug, PartialEq)]
pub struct Ladder {
    pub brand: String,
    pub product: String,
    pub name: String,
    pub sizes: Vec<u32>,
}

/// A gel box from `Gelbox.tdl`: model name, manufacturer, gel width and
/// length in cm, the well counts it supports and the maximum combs.
#[derive(Clone, Debug, PartialEq)]
pub struct GelBox {
    pub model: String,
    pub manufacturer: String,
    pub width_cm: f32,
    pub length_cm: f32,
    pub lanes: Vec<u32>,
    pub max_combs: u32,
}

/// Parse `Markers.tdl`: lines of `Brand \t Product \t Name \t sizes...`,
/// a comment header starting with `;`, trailing empty fields ignored.
pub fn parse_markers_tdl(text: &str) -> Vec<Ladder> {
    let mut out = vec![];
    for line in text.lines() {
        let line = line.trim_start_matches('\u{feff}');
        if line.is_empty() || line.starts_with(';') { continue; }
        let f: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
        if f.len() < 4 { continue; }
        let sizes: Vec<u32> = f[3..].iter().filter_map(|s| s.parse().ok()).collect();
        if sizes.is_empty() { continue; }
        out.push(Ladder { brand: f[0].into(), product: f[1].into(), name: f[2].into(), sizes });
    }
    out
}

/// Parse `Gelbox.tdl`: `Model \t Manufacturer \t Brand \t Width \t Length \t
/// Lanes|Lanes... \t Max Combs`. Lane counts are pipe-separated alternatives.
/// Line endings are normalised (the original file uses bare carriage returns).
pub fn parse_gelbox_tdl(text: &str) -> Vec<GelBox> {
    let text = text.replace("\r\n", "\n").replace('\r', "\n");
    let mut out = vec![];
    for line in text.lines() {
        let line = line.trim_start_matches('\u{feff}');
        if line.is_empty() || line.starts_with(';') { continue; }
        let f: Vec<&str> = line.split('\t').map(|s| s.trim()).collect();
        if f.len() < 7 { continue; }
        let (Ok(w), Ok(l), Ok(combs)) = (f[3].parse::<f32>(), f[4].parse::<f32>(), f[6].parse::<u32>()) else { continue };
        let lanes: Vec<u32> = f[5].split('|').filter_map(|s| s.trim().parse().ok()).collect();
        if lanes.is_empty() { continue; }
        out.push(GelBox { model: f[0].into(), manufacturer: f[1].into(), width_cm: w, length_cm: l, lanes, max_combs: combs });
    }
    out
}

/// Fragment lengths of a digest: top-strand cut positions (already resolved
/// through `cloning::hit_top_cut`, forward and reverse hits included) for a
/// circular or linear molecule. A circular molecule cut at 0 or 1 points
/// stays whole.
pub fn digest_fragments(seq_len: usize, circular: bool, cuts: &[usize]) -> Vec<u32> {
    let mut cuts: Vec<usize> = cuts.iter().copied().filter(|&c| c < seq_len).collect();
    cuts.sort_unstable();
    cuts.dedup();
    if cuts.is_empty() { return vec![seq_len as u32]; }
    let mut out = vec![];
    if circular {
        for w in cuts.windows(2) { out.push((w[1] - w[0]) as u32); }
        out.push((seq_len - cuts[cuts.len() - 1] + cuts[0]) as u32);
    } else {
        out.push(cuts[0] as u32);
        for w in cuts.windows(2) { out.push((w[1] - w[0]) as u32); }
        out.push((seq_len - cuts[cuts.len() - 1]) as u32);
    }
    out
}

/// Fraction of the gel run length a band of `bp` bases has migrated from the
/// well, on a gel resolving from `min_bp` to `max_bp`: log-linear, so bands
/// spread evenly across the ruler's range like on a real agarose gel.
pub fn migration_fraction(bp: u32, min_bp: u32, max_bp: u32) -> f32 {
    let (bp, min_bp, max_bp) = (bp.max(1) as f64, min_bp.max(1) as f64, max_bp.max(2) as f64);
    let t = (max_bp.log10() - bp.log10()) / (max_bp.log10() - min_bp.log10());
    t.clamp(0.0, 1.0) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    const MARKERS: &str = ";Brand\tProduct\tName\tsizes\nBio-Rad\t\t2.5 kb Molecular Ruler\t35000\t32500\t2500\nBio-Rad\t\t1 kb Molecular Ruler\t15000\t1000\nBioline\tHyperLadder\tHyperLadder 1kb\t10037\t8000\t200\n";
    const GELBOX: &str = ";Gel box\tManufacturer\tBrand\tWidth (cm)\tLength (cm)\tLanes\tMax Combs\rHGB-7 Horizontal Gel Box\tAxygen\tCorning Incorporated\t7\t7\t8\t2\rSub-Cell GT (15x10)\tBio-Rad\tBio-Rad Laboratories\t15\t10\t10|15|20|30\t4\r";

    #[test]
    fn markers_parse() {
        let l = parse_markers_tdl(MARKERS);
        assert_eq!(l.len(), 3);
        assert_eq!(l[0].name, "2.5 kb Molecular Ruler");
        assert_eq!(l[0].sizes, vec![35000, 32500, 2500]);
        assert_eq!(l[1].brand, "Bio-Rad");
        assert_eq!(l[1].product, "");
        assert_eq!(l[2].sizes, vec![10037, 8000, 200]);
    }

    #[test]
    fn gelbox_parse() {
        let b = parse_gelbox_tdl(GELBOX);
        assert_eq!(b.len(), 2);
        assert_eq!(b[0].model, "HGB-7 Horizontal Gel Box");
        assert_eq!(b[0].lanes, vec![8]);
        assert_eq!(b[1].width_cm, 15.0);
        assert_eq!(b[1].lanes, vec![10, 15, 20, 30]);
        assert_eq!(b[1].max_combs, 4);
    }

    #[test]
    fn bundled_data_is_readable() {
        let l = parse_markers_tdl(include_str!("../assets/Markers.tdl"));
        assert!(l.len() > 100, "Markers.tdl carries {} ladders", l.len());
        assert!(l.iter().any(|x| x.name.contains("1 kb DNA Ladder")));
        let b = parse_gelbox_tdl(include_str!("../assets/Gelbox.tdl"));
        assert!(b.len() > 10, "Gelbox.tdl carries {} boxes", b.len());
        assert!(b.iter().any(|x| x.model.contains("Sub-Cell GT")));
    }

    #[test]
    fn digest_fragments_linear_and_circular() {
        // linear 100 cut at 30 and 70: 30, 40, 30
        assert_eq!(digest_fragments(100, false, &[30, 70]), vec![30, 40, 30]);
        // circular 100 cut at 30 and 70: 40, 30 + wrap 60
        assert_eq!(digest_fragments(100, true, &[30, 70]), vec![40, 60]);
        // circular cut at 0 only: whole molecule
        assert_eq!(digest_fragments(100, true, &[0]), vec![100]);
        // no cuts: whole molecule
        assert_eq!(digest_fragments(100, false, &[]), vec![100]);
        // cuts outside the sequence are ignored
        assert_eq!(digest_fragments(100, false, &[140, 60]), vec![60, 40]);
    }

    #[test]
    fn migration_is_log_linear_and_ordered() {
        // on a 20 bp .. 48502 bp gel
        let far = migration_fraction(20, 20, 48502);
        let near = migration_fraction(48502, 20, 48502);
        assert_eq!(near, 0.0); // largest band stays at the well
        assert_eq!(far, 1.0); // smallest band runs to the bottom
        // a decade spans the same distance anywhere on the log scale
        let a = migration_fraction(5000, 20, 48502) - migration_fraction(500, 20, 48502);
        let b = migration_fraction(500, 20, 48502) - migration_fraction(50, 20, 48502);
        assert!((a - b).abs() < 1e-6);
        // clamped at the ends
        assert_eq!(migration_fraction(60000, 20, 48502), 0.0);
        assert_eq!(migration_fraction(10, 20, 48502), 1.0);
    }
}
