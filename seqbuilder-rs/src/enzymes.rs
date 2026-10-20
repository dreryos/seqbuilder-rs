// DNASTAR enzyme database (.ezd) and enzyme selector (.sel) files.
//
// .ezd: CR-separated records `ID;NAME:x;ISO:a,b;RS:G>AATT<C;CAV:Y;CUT:Pre;OWN:N;COS:48;UNI:50;VEN:NEB ...`.
//   RS: uppercase IUPAC = recognition site, lowercase n = flanking positions,
//   `>` = top-strand cut, `<` = bottom-strand cut (both in top-strand coordinates).
// .sel: XML <selectordata><selector><type/><id/><name/>...</selector>...; the field
//   meanings below come from the selector dialogs in SeqBuilder.xrc.

pub struct Enzyme {
    pub id: String,
    pub name: String,
    pub iso: String,
    pub rs: String,   // as stored, with cut markers
    pub site: String, // recognition sequence, uppercase IUPAC, flanking n's removed
    pub top_cut: i32, // cut positions relative to the first base of `site`
    pub bot_cut: i32,
    pub commercial: bool,
    pub on_hand: bool,
    pub vendors: String,
    degenerate: Degeneracy,
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Degeneracy { Specific, Partial, Complete }

impl Enzyme {
    /// (overhang length, true if 5' overhang); length 0 = blunt.
    pub fn overhang(&self) -> (i32, bool) {
        let d = self.bot_cut - self.top_cut;
        (d.abs(), d > 0)
    }

    // ponytail: SeqBuilder counts degenerate bases as full positions (AvaI CYCGRG shows as a six-cutter); N counts 0.
    fn complexity(&self) -> i32 {
        self.site.bytes().filter(|&b| b != b'N').count() as i32
    }
}

pub fn parse_ezd(bytes: &[u8]) -> Vec<Enzyme> {
    let text: String = bytes.iter().map(|&b| b as char).collect(); // Latin-1
    text.split(['\r', '\n']).filter(|l| !l.is_empty() && !l.starts_with(';')).filter_map(|line| {
        let mut parts = line.split(';');
        let id = parts.next()?.to_string();
        let field = |k: &str| line.split(';').find_map(|p| p.strip_prefix(k).and_then(|v| v.strip_prefix(':'))).unwrap_or("");
        let rs = field("RS").to_string();
        let (mut bases, mut top, mut bot) = (Vec::new(), 0i32, 0i32);
        for b in rs.bytes() {
            match b {
                b'>' => top = bases.len() as i32,
                b'<' => bot = bases.len() as i32,
                _ => bases.push(b),
            }
        }
        let lead = bases.iter().take_while(|b| b.is_ascii_lowercase()).count();
        let trail = bases[lead..].iter().rev().take_while(|b| b.is_ascii_lowercase()).count();
        let site = String::from_utf8_lossy(&bases[lead..bases.len() - trail]).to_ascii_uppercase();
        if site.is_empty() { return None; }
        let (a, b) = (top.min(bot).max(0) as usize, top.max(bot).min(bases.len() as i32) as usize);
        let over = &bases[a..b];
        let degenerate = if over.iter().all(|c| c.is_ascii_lowercase()) && !over.is_empty() { Degeneracy::Complete }
            else if over.iter().all(|c| b"ACGT".contains(c)) { Degeneracy::Specific }
            else { Degeneracy::Partial };
        Some(Enzyme {
            name: field("NAME").to_string(), iso: field("ISO").to_string(), site,
            top_cut: top - lead as i32, bot_cut: bot - lead as i32,
            commercial: field("CAV") == "Y", on_hand: field("OWN") == "Y",
            vendors: field("VEN").to_string(), rs, id, degenerate,
        })
    }).collect()
}

pub enum Rule {
    Manual(Vec<String>), // enzyme IDs
    Overhang { use_len: bool, min: i32, use_max: bool, max: i32, five_prime: bool, degeneracy: u8 },
    Frequency { min: usize, use_max: bool, max: usize, range: Option<(usize, usize)> },
    Class { use_cx: bool, min_cx: i32, use_max_cx: bool, max_cx: i32, commercial: bool, on_hand: bool },
    Combine(Vec<(u8, u32)>), // (operator, selector id); operators: 0 first, 1 not, 2 or, 3 and, 4 or not, 5 and not
}

pub struct Selector {
    pub id: u32,
    pub name: String,
    pub rule: Rule,
}

pub fn parse_sel(text: &str) -> Result<Vec<Selector>, String> {
    let doc = roxmltree::Document::parse(text).map_err(|e| e.to_string())?;
    Ok(doc.root_element().children().filter(|n| n.has_tag_name("selector")).filter_map(|s| {
        let all = |tag: &str| -> Vec<&str> { s.children().filter(|c| c.has_tag_name(tag)).map(|c| c.text().unwrap_or("").trim()).collect() };
        let txt = |tag: &str| all(tag).first().copied().unwrap_or("");
        let num = |tag: &str| txt(tag).parse::<i64>().unwrap_or(0);
        let flag = |tag: &str| num(tag) != 0;
        let rule = match num("type") {
            1 => Rule::Manual(all("manenzyid").into_iter().map(str::to_owned).collect()),
            2 => Rule::Overhang {
                use_len: flag("useoverlength"), min: num("minoverhang") as i32, use_max: flag("usemaxover"),
                max: num("maxoverhang") as i32, five_prime: flag("overdir"), degeneracy: num("degeneracy") as u8,
            },
            3 => Rule::Frequency {
                min: num("mincuts") as usize, use_max: flag("usemaxcuts"), max: num("maxcuts") as usize,
                range: flag("userange").then(|| (num("lowbase") as usize, num("highbase") as usize)),
            },
            // ponytail: precise/random class and cost filters are ignored; every shipped selector leaves them off.
            4 => Rule::Class {
                use_cx: flag("usecomplexity"), min_cx: num("mincomplexity") as i32, use_max_cx: flag("usemaxcomplexity"),
                max_cx: num("maxcomplexity") as i32, commercial: flag("available"), on_hand: flag("onhand"),
            },
            5 => Rule::Combine(all("comboper").into_iter().zip(all("combid")).map(|(o, i)| (o.parse().unwrap_or(0), i.parse().unwrap_or(0))).collect()),
            _ => return None,
        };
        Some(Selector { id: num("id") as u32, name: txt("name").to_string(), rule })
    }).collect())
}

/// Human-readable definition lines, as SeqBuilder's Selector Manager shows them.
pub fn describe(sels: &[Selector], s: &Selector) -> Vec<String> {
    let name = |id: u32| sels.iter().find(|s| s.id == id).map_or(format!("#{id}"), |s| s.name.trim().to_string());
    let range = |on: bool, min: String, use_max: bool, max: String| match (on, use_max) {
        (false, _) => "any".to_string(),
        (true, true) => format!("{min}..{max}"),
        (true, false) => format!("≥ {min}"),
    };
    match &s.rule {
        Rule::Manual(ids) => vec!["Type: Manual".into(), format!("{} enzymes listed", ids.len())],
        Rule::Overhang { use_len, min, use_max, max, five_prime, degeneracy } => vec![
            "Type: Overhang".into(),
            format!("Overhang length: {} ({})", range(*use_len, min.to_string(), *use_max, max.to_string()), if *five_prime { "5'" } else { "3'" }),
            format!("Degeneracy: {}", ["Specific only", "Degenerate only", "Both", "Golden Gate"].get(*degeneracy as usize).unwrap_or(&"?")),
        ],
        Rule::Frequency { min, use_max, max, range: r } => {
            let mut v = vec!["Type: Frequency".into(), format!("Cut frequency: {}", range(true, min.to_string(), *use_max, max.to_string()))];
            if let Some((lo, hi)) = r { v.push(format!("Sequence range: {lo}..{hi}")); }
            v
        }
        Rule::Class { use_cx, min_cx, use_max_cx, max_cx, commercial, on_hand } => vec![
            "Type: Class".into(),
            format!("Site complexity: {}", range(*use_cx, min_cx.to_string(), *use_max_cx, max_cx.to_string())),
            format!("Commercial: {}   On hand: {}", if *commercial { "yes" } else { "any" }, if *on_hand { "yes" } else { "any" }),
        ],
        Rule::Combine(items) => std::iter::once("Type: Combine".to_string()).chain(items.iter().map(|&(op, id)| {
            format!("{} {}", ["", "not", "or", "and", "or not", "and not"].get(op as usize).unwrap_or(&"?"), name(id)).trim().to_string()
        })).collect(),
    }
}

/// Does enzyme `e` pass selector `id`? `hits` are the enzyme's 0-based site positions in the current sequence.
pub fn selects(sels: &[Selector], id: u32, e: &Enzyme, hits: &[usize]) -> bool {
    eval(sels, id, e, hits, 0)
}

fn eval(sels: &[Selector], id: u32, e: &Enzyme, hits: &[usize], depth: u8) -> bool {
    let Some(s) = sels.iter().find(|s| s.id == id) else { return false };
    if depth > 16 { return false; } // selectors referencing each other in a loop
    match &s.rule {
        Rule::Manual(ids) => ids.iter().any(|i| *i == e.id),
        &Rule::Overhang { use_len, min, use_max, max, five_prime, degeneracy } => {
            let (len, five) = e.overhang();
            let len_ok = !use_len || (len >= min && (!use_max || len <= max));
            let dir_ok = !use_len || len == 0 || five == five_prime;
            let deg_ok = match degeneracy {
                0 => e.degenerate == Degeneracy::Specific,
                1 => e.degenerate != Degeneracy::Specific,
                3 => e.degenerate == Degeneracy::Complete,
                _ => true,
            };
            len_ok && dir_ok && deg_ok
        }
        &Rule::Frequency { min, use_max, max, range } => {
            let n = match range {
                Some((lo, hi)) => hits.iter().filter(|&&p| p + 1 >= lo && p < hi).count(),
                None => hits.len(),
            };
            n >= min && (!use_max || n <= max)
        }
        &Rule::Class { use_cx, min_cx, use_max_cx, max_cx, commercial, on_hand } => {
            let cx = e.complexity();
            (!use_cx || (cx >= min_cx && (!use_max_cx || cx <= max_cx))) && (!commercial || e.commercial) && (!on_hand || e.on_hand)
        }
        Rule::Combine(items) => items.iter().fold(false, |acc, &(op, sid)| {
            let v = eval(sels, sid, e, hits, depth + 1);
            match op { 0 => v, 1 => !v, 2 => acc || v, 3 => acc && v, 4 => acc || !v, 5 => acc && !v, _ => acc }
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_cut_markers() {
        let e = parse_ezd(b";Enzyme file version 18.0.0\rECO1;NAME:EcoRI;ISO:;RS:G>AATT<C;CAV:Y;OWN:N\rBSA1;NAME:BsaI;ISO:;RS:GGTCTCn>nnnn<;CAV:Y;OWN:N\rPST1;NAME:PstI;ISO:;RS:C<TGCA>G;CAV:Y;OWN:N\r");
        assert_eq!(e.len(), 3);
        assert_eq!((e[0].site.as_str(), e[0].top_cut, e[0].bot_cut, e[0].overhang()), ("GAATTC", 1, 5, (4, true)));
        assert_eq!((e[1].site.as_str(), e[1].overhang(), e[1].degenerate), ("GGTCTC", (4, true), Degeneracy::Complete));
        assert_eq!(e[2].overhang(), (4, false));
    }
}
