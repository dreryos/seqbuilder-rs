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
    text.replace('\\', r"\\").replace('(', r"\(").replace(')', r"\)")
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

    /// Fill a rectangle with an RGB colour (0-255 per channel).
    pub fn rect(&mut self, x: f32, y: f32, w: f32, h: f32, rgb: (u8, u8, u8)) {
        let (r, g, b) = (rgb.0 as f32 / 255.0, rgb.1 as f32 / 255.0, rgb.2 as f32 / 255.0);
        self.ops.push(format!("{r} {g} {b} rg {x} {y} {w} {h} re f"));
    }

    /// Assemble the PDF file (A4 portrait, 595×842pt).
    pub fn finish(mut self) -> Vec<u8> {
        if !self.ops.is_empty() { self.flush(); }
        let n_pages = self.pages.len();
        let mut objs: Vec<String> = vec![]; // object bodies, 1-indexed

        // 1: catalog, 2: pages tree; then per page: content stream + page obj
        let mut page_ids = vec![];
        let mut next_id = 3;
        for content in &self.pages {
            let stream_id = next_id;
            objs.push(format!("<< /Length {} >>\nstream\n{}\nendstream", content.iter().map(|s| s.len() + 1).sum::<usize>(), content.join("\n")));
            let page_id = next_id + 1;
            objs.push(format!("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 1000 0 R /F2 1001 0 R /F3 1002 0 R >> >> /Contents {stream_id} 0 R >>"));
            page_ids.push(page_id);
            next_id += 2;
        }
        let kids: Vec<String> = page_ids.iter().map(|id| format!("{id} 0 R")).collect();
        let pages_obj = format!("<< /Type /Pages /Kids [{}] /Count {n_pages} >>", kids.join(" "));

        let mut out = format!("%PDF-1.4\n").into_bytes();
        let mut offsets = vec![0usize]; // byte offset of object i (1-indexed)
        let mut write_obj = |out: &mut Vec<u8>, offsets: &mut Vec<usize>, id: usize, body: &str| {
            offsets.push(out.len());
            out.extend_from_slice(format!("{id} 0 obj\n{body}\nendobj\n").as_bytes());
        };
        write_obj(&mut out, &mut offsets, 1, "<< /Type /Catalog /Pages 2 0 R >>");
        write_obj(&mut out, &mut offsets, 2, &pages_obj);
        for (i, body) in objs.iter().enumerate() {
            write_obj(&mut out, &mut offsets, 3 + i, body);
        }
        // base-14 fonts at fixed ids 1000-1002
        write_obj(&mut out, &mut offsets, 1000, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>");
        write_obj(&mut out, &mut offsets, 1001, "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>");
        write_obj(&mut out, &mut offsets, 1002, "<< /Type /Font /Subtype /Type1 /BaseFont /Courier >>");

        let xref_pos = out.len();
        let max_id = 1003;
        out.extend_from_slice(format!("xref\n0 {max_id}\n").as_bytes());
        out.extend_from_slice(b"0000000000 65535 f \n");
        for id in 1..max_id {
            let off = offsets.get(id).copied().unwrap_or(0);
            if off > 0 && (id <= 2 + objs.len() || id >= 1000) {
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
    const COL_W: f32 = 7.2; // Courier advance at 12pt is 7.2
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

    // feature table
    y -= ROW_H;
    if y < BOTTOM { pdf.page(); y = TOP; }
    pdf.text(40.0, y, 12.0, Font::HelveticaBold, "Features");
    y -= ROW_H;
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
        pdf.line(40.0, 790.0, 500.0, 790.0);
        let bytes = pdf.finish();
        let s = String::from_utf8_lossy(&bytes);
        assert!(s.starts_with("%PDF-1.4"));
        assert!(s.contains("/Type /Catalog"));
        assert!(s.contains("(ACGT) Tj"));
        assert!(s.ends_with("%%EOF"));
        assert!(s.contains("startxref"));
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
        assert!(s.contains("pBR322 — circular — 4200 bp"));
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
    }
}
