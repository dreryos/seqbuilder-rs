// GenBank flat-file writer laid out like NCBI's: the header is kept as written (only LOCUS is
// regenerated for the current length and topology), features come from qualifier lines that already
// carry their quoting, the sequence is written 60 bases a line in groups of 10.

pub struct Feature<'a> {
    pub kind: &'a str,
    pub location: String,
    pub quals: &'a [String], // "key=value" lines, quoting included
}

const QUAL_WIDTH: usize = 58; // columns 22..79

/// `header` is the GenBank text before FEATURES (may be empty); `name` is used when it has no LOCUS.
pub fn write(header: &str, name: &str, seq: &[u8], circular: bool, features: &[Feature]) -> String {
    let mut out = String::new();
    let mut lines = header.lines().filter(|l| !l.is_empty()).peekable();
    let old_locus = lines.next_if(|l| l.starts_with("LOCUS"));
    out += &locus(old_locus, name, seq.len(), circular);
    for l in lines { out += l; out.push('\n'); }

    out += "FEATURES             Location/Qualifiers\n";
    for f in features {
        out += &format!("     {:<15} {}\n", f.kind, f.location);
        for q in f.quals {
            for chunk in wrap(&format!("/{q}"), q.starts_with("translation=")) {
                out += &format!("{:21}{chunk}\n", "");
            }
        }
    }

    out += "ORIGIN      \n";
    for (row, bases) in seq.chunks(60).enumerate() {
        let groups: Vec<String> = bases.chunks(10).map(|g| String::from_utf8_lossy(g).to_lowercase()).collect();
        out += &format!("{:>9} {}\n", row * 60 + 1, groups.join(" "));
    }
    out += "//\n";
    out
}

/// LOCUS with the current length and topology; name, molecule, division and date are kept when known.
fn locus(old: Option<&str>, name: &str, len: usize, circular: bool) -> String {
    let words: Vec<&str> = old.map(|l| l.split_whitespace().collect()).unwrap_or_default();
    let bp = words.iter().position(|w| *w == "bp" || *w == "aa");
    let field = |k: usize| bp.and_then(|b| words.get(b + k)).copied();
    let name = words.get(1).copied().filter(|_| bp != Some(1)).unwrap_or(name);
    let mol = field(1).unwrap_or("DNA");
    let (div, date) = match (field(3), field(4)) {
        (Some(d), Some(t)) => (d.to_string(), t.to_string()),
        _ => ("UNK".to_string(), today()),
    };
    let topology = if circular { "circular" } else { "linear" };
    // name and length share columns 13..40, the length right-aligned
    let len = len.to_string();
    let gap = " ".repeat(28usize.saturating_sub(name.len() + len.len()).max(1));
    format!("LOCUS       {name}{gap}{len} bp    {mol:<7} {topology:<8} {div} {date}\n")
}

/// Break a qualifier into lines of at most 58 columns: at spaces for text, anywhere for translations.
fn wrap(text: &str, hard: bool) -> Vec<String> {
    let mut lines = vec![];
    let mut rest: Vec<char> = text.chars().collect();
    while rest.len() > QUAL_WIDTH {
        let cut = if hard { None } else { rest[..=QUAL_WIDTH].iter().rposition(|&c| c == ' ') };
        match cut {
            Some(at) if at > 0 => {
                lines.push(rest[..at].iter().collect());
                rest.drain(..=at);
            }
            _ => {
                lines.push(rest[..QUAL_WIDTH].iter().collect());
                rest.drain(..QUAL_WIDTH);
            }
        }
    }
    lines.push(rest.into_iter().collect());
    lines
}

/// Today as GenBank writes dates (21-OCT-2026), from the system clock in UTC.
fn today() -> String {
    let days = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs() / 86_400) as i64;
    // civil-from-days (Howard Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + (month <= 2) as i64;
    const MONTHS: [&str; 12] = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"];
    format!("{day:02}-{}-{year}", MONTHS[month as usize - 1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wraps_like_ncbi() {
        let long = "note=\"corresponds to one of the 38bp repeats found in Tn3 (bp 1-38 of Tn3)\"";
        assert_eq!(wrap(&format!("/{long}"), false), ["/note=\"corresponds to one of the 38bp repeats found in Tn3", "(bp 1-38 of Tn3)\""]);
        assert!(locus(None, "ALDH2", 114842, false).starts_with(&format!("LOCUS       ALDH2{}114842 bp    DNA     linear   UNK ", " ".repeat(17))));
    }
}
