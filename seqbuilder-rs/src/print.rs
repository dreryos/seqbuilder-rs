//! Printing (PrintCommands / PageForm in SeqBuilder Pro): the printed
//! sequence layout rendered as a PDF — header with name, topology and
//! length, the sequence in numbered 60-base rows grouped in tens, and the
//! feature table.

/// Fonts of the base-14 set the PDF viewer provides.
#[derive(Clone, Copy, PartialEq)]
pub enum Font { Helvetica, HelveticaBold, Courier }

pub struct Pdf {
    pages: Vec<Vec<String>>, // content streams
    ops: Vec<String>,         // current page
}

impl Default for Pdf {
    fn default() -> Self { Self::new() }
}

fn escape(text: &str) -> String {
    // PDF text strings are single-byte for the base-14 fonts: keep the
    // escapes, transliterate the few punctuation marks we emit, and replace
    // anything else non-ASCII with '?'
    let mut out = String::new();
    for ch in text.chars() {
        match ch {
            '\\' => out.push_str(r"\\"),
            '(' => out.push_str(r"\("),
            ')' => out.push_str(r"\)"),
            '\t' => out.push(' '),
            '\r' | '\n' => out.push(' '),
            '—' => out.push('-'),
            '–' => out.push('-'),
            '…' => out.push_str("..."),
            c if (c as u32) < 128 => out.push(c),
            _ => out.push('?'),
        }
    }
    out
}

impl Pdf {
    pub fn new() -> Self { Pdf { pages: vec![], ops: vec![] } }

    /// Start a new page; operations drawn before this land on none.
    pub fn page(&mut self) {
        if !self.ops.is_empty() { self.flush(); }
        self.ops.clear();
    }

    fn flush(&mut self) {
        self.pages.push(std::mem::take(&mut self.ops));
    }

    /// Draw `text` with the top-left corner at (x, y), in points measured
    /// from the bottom-left of the page.
    pub fn text(&mut self, x: f32, y_from_bottom: f32, size: f32, font: Font, text: &str) {
        let f = match font { Font::Helvetica => "/F1", Font::HelveticaBold => "/F2", Font::Courier => "/F3" };
        self.ops.push(format!(
            "BT {f} {size} Tf {x} {y_from_bottom} Td ({}) Tj ET", escape(text)));
    }

    /// Draw a 0.75pt line between two points.
    pub fn line(&mut self, x1: f32, y1: f32, x2: f32, y2: f32) {
        self.ops.push(format!("0.75 w {x1} {y1} m {x2} {y2} l S"));
    }

    /// Fill a rectangle with an RGB colour (0-255 per channel), then restore
    /// black as the fill colour so following text stays readable.
    pub fn rect(&mut self, x: f32, y: f32, w: f32, h: f32, rgb: (u8, u8, u8)) {
        let (r, g, b) = (rgb.0 as f32 / 255.0, rgb.1 as f32 / 255.0, rgb.2 as f32 / 255.0);
        self.ops.push(format!("{r} {g} {b} rg {x} {y} {w} {h} re f 0 0 0 rg"));
    }

    /// Assemble the PDF file (A4 portrait, 595×842pt).
    pub fn finish(mut self) -> Vec<u8> {
        if !self.ops.is_empty() { self.flush(); }
        let n_pages = self.pages.len();
        // ids: 1 catalog, 2 pages tree, then per page (content stream, page),
        // finally the three base-14 fonts — all sequential, so every xref
        // entry has a real offset
        let f1 = 3 + 2 * n_pages;
        let (f2, f3) = (f1 + 1, f1 + 2);
        let max_id = f3 + 1;

        let mut out = b"%PDF-1.4\n".to_vec();
        let mut offsets = vec![0usize; max_id]; // offsets[id], objects are 1-based
        let write = |out: &mut Vec<u8>, offsets: &mut [usize], id: usize, body: String| {
            offsets[id] = out.len();
            out.extend_from_slice(format!("{id} 0 obj\n{body}\nendobj\n").as_bytes());
        };
        write(&mut out, &mut offsets, 1, "<< /Type /Catalog /Pages 2 0 R >>".into());
        let mut kids: Vec<String> = vec![];
        for (i, content) in self.pages.iter().enumerate() {
            let stream_id = 3 + 2 * i;
            let stream: String = content.iter().map(|s| s.to_string()).collect::<Vec<_>>().join("\n");
            write(&mut out, &mut offsets, stream_id, format!("<< /Length {} >>\nstream\n{stream}\nendstream", stream.len() + 1));
            let page_id = stream_id + 1;
            write(&mut out, &mut offsets, page_id, format!(
                "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 {f1} 0 R /F2 {f2} 0 R /F3 {f3} 0 R >> >> /Contents {stream_id} 0 R >>"));
            kids.push(format!("{page_id} 0 R"));
        }
        write(&mut out, &mut offsets, 2, format!("<< /Type /Pages /Kids [{}] /Count {n_pages} >>", kids.join(" ")));
        write(&mut out, &mut offsets, f1, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>".into());
        write(&mut out, &mut offsets, f2, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>".into());
        write(&mut out, &mut offsets, f3, "<< /Type /Font /Subtype /Type1 /BaseFont /Courier >>".into());

        let xref_pos = out.len();
        out.extend_from_slice(format!("xref\n0 {max_id}\n").as_bytes());
        out.extend_from_slice(b"0000000000 65535 f \n");
        for id in 1..max_id {
            let off = offsets[id];
            if off > 0 {
                out.extend_from_slice(format!("{off:010} 00000 n \n").as_bytes());
            } else {
                out.extend_from_slice(b"0000000000 65535 f \n");
            }
        }
        out.extend_from_slice(format!("trailer\n<< /Size {max_id} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF").as_bytes());
        out
    }
}

/// The document header shown on page 1 of the printout.
pub struct PrintDoc<'a> {
    pub name: &'a str,
    pub topology: &'a str,
    pub seq: &'a [u8],
    pub features: &'a [(String, u32, u32, i32)], // (name, start, end, strand)
    pub first_base: usize, // 0-based index of the first printed base
}

/// Lay the document out into `pdf`, returning the number of pages used.
/// A4 portrait: header, sequence in numbered 60-base rows (blocks of 10),
/// then the feature table.
pub fn layout(pdf: &mut Pdf, doc: &PrintDoc) -> usize {
    const W: f32 = 595.0;
    const TOP: f32 = 800.0;
    const BOTTOM: f32 = 60.0;
    const ROW_H: f32 = 12.0;
    pdf.page();
    let mut y = TOP;
    let header = format!("{} — {} — {} bp", doc.name, doc.topology, doc.seq.len());
    pdf.text(40.0, y, 14.0, Font::HelveticaBold, &header);
    y -= 10.0;
    pdf.line(40.0, y, W - 40.0, y);
    y -= 20.0;

    let row_len = 60;
    let seq = std::str::from_utf8(doc.seq).unwrap_or("");
    let seq: String = seq.bytes().filter(|b| !b"\r\n ".contains(b)).map(|b| b as char).collect();
    let total = doc.first_base + seq.len();
    let mut pos = doc.first_base;
    while pos < total {
        let row = &seq[pos - doc.first_base..(pos - doc.first_base + row_len).min(seq.len())];
        let mut line = String::new();
        for (i, chunk) in row.as_bytes().chunks(10).enumerate() {
            if i > 0 { line.push(' '); }
            line.push_str(&String::from_utf8_lossy(chunk));
        }
        pdf.text(40.0, y, 12.0, Font::Courier, &format!("{:>6}  {}", pos + 1, line));
        y -= ROW_H;
        pos += row_len;
        if y < BOTTOM {
            pdf.page();
            y = TOP;
        }
    }

    // feature table, its header on a light band
    y -= ROW_H;
    if y < BOTTOM { pdf.page(); y = TOP; }
    pdf.text(40.0, y, 12.0, Font::HelveticaBold, "Features");
    y -= ROW_H;
    pdf.rect(40.0, y - 3.0, W - 80.0, 12.0, (232, 236, 242));
    pdf.text(40.0, y, 9.0, Font::Helvetica, "Name                                   Range         Strand");
    y -= ROW_H;
    for (name, start, end, strand) in doc.features {
        if y < BOTTOM { pdf.page(); y = TOP; }
        let s = format!("{:<38} {:>6}..{:<6} {}", name, start + 1, end, if *strand < 0 { '-' } else { '+' });
        pdf.text(40.0, y, 9.0, Font::Courier, &s);
        y -= ROW_H;
    }
    // the page under construction is only flushed by finish()
    pdf.pages.len() + usize::from(!pdf.ops.is_empty())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pdf_has_valid_structure() {
        let mut pdf = Pdf::new();
        pdf.page();
        pdf.text(40.0, 800.0, 12.0, Font::Courier, "ACGT");
        pdf.text(40.0, 780.0, 12.0, Font::Courier, "pUC — x…y"); // non-ASCII must not corrupt the string
        pdf.line(40.0, 790.0, 500.0, 790.0);
        let bytes = pdf.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.starts_with("%PDF-1.4"));
        assert!(s.contains("/Type /Catalog"));
        assert!(s.contains("(ACGT) Tj"));
        assert!(s.contains("(pUC - x...y) Tj"));
        assert!(s.ends_with("%%EOF"));
        assert!(s.contains("startxref"));
        // every in-use xref entry must resolve to that exact object header
        let sx = s.rfind("startxref").unwrap();
        let xref_pos: usize = s[sx + 10..].trim().lines().next().unwrap().parse().unwrap();
        assert!(s[xref_pos..].starts_with("xref"));
        let xref = &s[xref_pos..xref_pos + s[xref_pos..].find("trailer").unwrap()];
        let mut entries = xref.lines();
        assert_eq!(entries.next().unwrap(), "xref");
        let count: usize = entries.next().unwrap().split_whitespace().nth(1).unwrap().parse().unwrap();
        let lines: Vec<&str> = entries.collect();
        assert_eq!(lines.len(), count);
        assert!(lines[0].trim().ends_with('f'));
        for (id, line) in lines.iter().enumerate().skip(1) {
            if line.trim().ends_with('n') {
                let off: usize = line[..10].parse().unwrap();
                assert!(s[off..].starts_with(&format!("{id} 0 obj")), "xref entry for object {id} points at offset {off}: {}", &s[off..off + 24]);
            }
        }
    }

    #[test]
    fn layout_pages_a_long_sequence() {
        let seq = vec![b'A'; 60 * 70]; // 70 rows -> multiple pages
        let feats = vec![("bla".to_string(), 100u32, 900u32, -1i32)];
        let doc = PrintDoc { name: "pBR322", topology: "circular", seq: &seq, features: &feats, first_base: 0 };
        let mut pdf = Pdf::new();
        let pages = layout(&mut pdf, &doc);
        assert!(pages >= 2, "70 rows do not fit one page: {pages}");
        let bytes = pdf.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("pBR322 - circular - 4200 bp"));
        assert!(s.contains("Features"));
        assert!(s.contains("101..900"));
    }

    #[test]
    fn print_selection_offsets_rows() {
        let seq = b"ACGTACGTAC";
        let doc = PrintDoc { name: "sel", topology: "linear", seq, features: &[], first_base: 50 };
        let mut pdf = Pdf::new();
        layout(&mut pdf, &doc);
        let bytes = pdf.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.contains("    51  ACGTACGTAC"), "row starts at 51: {}", s.lines().take(8).collect::<Vec<_>>().join(" | "));
        // the feature table's header band: the light fill, its full width
        // (515pt tall 12pt) and the black restoration that keeps the text
        // over it readable
        let (colour, band) = (s.find("0.9098039 0.9254902 0.9490196 rg"), s.find("515 12 re f"));
        let (colour, band) = (colour.expect("header fill colour"), band.expect("header band dimensions"));
        assert!(colour < band, "the fill colour must be set before the rectangle is filled");
        assert!(s.contains("515 12 re f 0 0 0 rg"), "fill colour restored to black after the band");
    }
}
