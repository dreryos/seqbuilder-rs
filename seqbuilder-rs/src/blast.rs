//! Online search (SeqBuilder Pro's "Search for Sequences…"): submit the
//! sequence to NCBI's BLAST URL API, poll for completion and parse the XML
//! result. The original scraped the returned HTML through the rule tables
//! in BlastParse.xml; the URL API's XML output carries the same hits
//! without scraping.

const BLAST_URL: &str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";

/// One database hit with its best HSP.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct BlastHit {
    pub accession: String,
    pub definition: String,
    pub identity: f32,   // percent identical positions
    pub length: u32,     // alignment length
    pub bits: f32,       // bit score
    pub evalue: String,
    pub query_from: u32,
    pub query_to: u32,
    pub hit_from: u32,
    pub hit_to: u32,
    pub qseq: String,
    pub hseq: String,
    pub midline: String,
}

/// A finished BLAST run.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct BlastResult {
    pub rid: String,
    pub db: String,
    pub program: String,
    pub query_length: u32,
    pub hits: Vec<BlastHit>,
}

fn urlencode(s: &str) -> String {
    let mut out = String::new();
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' => out.push(b as char),
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

/// Submit a BLAST job; returns the RID used to poll for results. The query
/// goes in a POST body — long sequences would exceed any GET URL limit.
pub fn submit(seq: &[u8], db: &str, program: &str, hitlist: usize) -> Result<String, String> {
    let seq = std::str::from_utf8(seq).map_err(|_| "the sequence is not ASCII text".to_string())?
        .bytes().filter(|b| !b" \r\n\t".contains(b)).map(|b| b as char).collect::<String>();
    if seq.is_empty() { return Err("the sequence is empty".into()); }
    let body = format!("CMD=Put&PROGRAM={}&DATABASE={}&HITLIST_SIZE={}&QUERY={}",
        urlencode(program), urlencode(db), hitlist, urlencode(&seq));
    let resp = ureq::post(BLAST_URL)
        // submissions carry at most a few kb and NCBI answers within
        // seconds; a dropped dialog can only cancel this worker between
        // requests, so the request must not stay in flight for two minutes
        .timeout(std::time::Duration::from_secs(30))
        .set("Content-Type", "application/x-www-form-urlencoded")
        .send_string(&body).map_err(|e| e.to_string())?;
    let body = resp.into_string().map_err(|e| e.to_string())?;
    for line in body.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("RID = ") {
            let rid = rest.trim().to_string();
            if !rid.is_empty() { return Ok(rid); }
        }
    }
    Err("NCBI did not return a request ID".into())
}

/// Poll a submitted job: Ok(None) while it is still running, Ok(Some) with
/// the parsed result once it finishes. UNKNOWN (expired RID) and FAILED are
/// terminal errors.
pub fn poll(rid: &str) -> Result<Option<BlastResult>, String> {
    let url = format!("{BLAST_URL}?CMD=Get&FORMAT_TYPE=XML&RID={}", urlencode(rid));
    let body = http_get(&url)?;
    if body.contains("<BlastOutput") {
        return Ok(Some(parse_xml(&body, rid)?));
    }
    if body.contains("Status=WAITING") {
        return Ok(None);
    }
    if body.contains("Status=UNKNOWN") {
        return Err(format!("NCBI does not know RID {rid} (it may have expired)"));
    }
    if body.contains("Status=FAILED") {
        return Err("NCBI reported the BLAST job as failed".into());
    }
    Err(format!("unexpected BLAST response for RID {rid}"))
}

fn http_get(url: &str) -> Result<String, String> {
    // a short timeout: an NCBI poll answers within seconds, and a dropped
    // dialog can only cancel the worker between requests, so requests must
    // not stay in flight for the old two minutes
    // transient transport failures (a slow answer, a reset connection) are
    // retried — the NCBI job keeps running through them, and the worker's
    // overall deadline bounds the wait; HTTP status errors are terminal
    let mut last = String::new();
    for attempt in 0..3 {
        match ureq::get(url).timeout(std::time::Duration::from_secs(30)).call() {
            Ok(resp) => return resp.into_string().map_err(|e| e.to_string()),
            // transport errors (timeout, reset, DNS) are transient: retry
            Err(e) if attempt < 2 && matches!(e, ureq::Error::Transport(_)) => {
                last = e.to_string();
                std::thread::sleep(std::time::Duration::from_secs(5));
            }
            // HTTP status errors are terminal
            Err(e) => return Err(e.to_string()),
        }
    }
    Err(last)
}

/// Parse a BLAST XML report (CMD=Get&FORMAT_TYPE=XML) into hits.
pub fn parse_xml(xml: &str, rid: &str) -> Result<BlastResult, String> {
    // roxmltree refuses documents carrying a DTD; NCBI's report has one
    let xml: String = xml.lines().filter(|l| !l.trim_start().starts_with("<!DOCTYPE")).collect::<Vec<_>>().join("\n");
    let doc = roxmltree::Document::parse(&xml).map_err(|e| e.to_string())?;
    let root = doc.root_element();
    let text = |n: &roxmltree::Node, tag: &str| n.children().find(|c| c.has_tag_name(tag))
        .and_then(|c| c.text()).unwrap_or("").trim().to_string();
    let num = |n: &roxmltree::Node, tag: &str| text(n, tag).parse::<f32>().unwrap_or(0.0);
    let mut res = BlastResult {
        rid: rid.to_string(),
        db: text(&root, "BlastOutput_db"),
        program: text(&root, "BlastOutput_program"),
        query_length: text(&root, "BlastOutput_query-length").parse().unwrap_or(0),
        hits: vec![],
    };
    for it in root.descendants().filter(|n| n.has_tag_name("Iteration")) {
        for hit in it.children().filter(|n| n.has_tag_name("Iteration_hits"))
            .flat_map(|h| h.children().filter(|c| c.has_tag_name("Hit")))
        {
            let accession = text(&hit, "Hit_accession");
            let definition = text(&hit, "Hit_def");
            let Some(hsp) = hit.descendants().find(|n| n.has_tag_name("Hsp")) else { continue };
            res.hits.push(BlastHit {
                accession, definition,
                identity: num(&hsp, "Hsp_identity"),
                length: num(&hsp, "Hsp_align-len") as u32,
                bits: num(&hsp, "Hsp_bit-score"),
                evalue: text(&hsp, "Hsp_evalue"),
                query_from: num(&hsp, "Hsp_query-from") as u32,
                query_to: num(&hsp, "Hsp_query-to") as u32,
                hit_from: num(&hsp, "Hsp_hit-from") as u32,
                hit_to: num(&hsp, "Hsp_hit-to") as u32,
                qseq: text(&hsp, "Hsp_qseq"),
                hseq: text(&hsp, "Hsp_hseq"),
                midline: text(&hsp, "Hsp_midline"),
            });
        }
    }
    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = r#"<?xml version="1.0"?>
<!DOCTYPE BlastOutput>
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>pBR322</BlastOutput_query-def>
  <BlastOutput_query-length>4361</BlastOutput_query-length>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gnl|BL_ORD_ID|0</Hit_id>
          <Hit_def>Escherichia coli plasmid pBR322, complete sequence</Hit_def>
          <Hit_accession>U00004</Hit_accession>
          <Hit_len>4361</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>8010</Hsp_bit-score>
              <Hsp_score>4185</Hsp_score>
              <Hsp_evalue>0.0</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>4361</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>4361</Hsp_hit-to>
              <Hsp_identity>4361</Hsp_identity>
              <Hsp_positive>4361</Hsp_positive>
              <Hsp_align-len>4361</Hsp_align-len>
              <Hsp_qseq>ACGT</Hsp_qseq>
              <Hsp_hseq>ACGT</Hsp_hseq>
              <Hsp_midline>||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_def>synthetic construct cloning vector</Hit_def>
          <Hit_accession>JX128753</Hit_accession>
          <Hit_len>2700</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_bit-score>2400</Hsp_bit-score>
              <Hsp_evalue>1e-180</Hsp_evalue>
              <Hsp_query-from>3293</Hsp_query-from>
              <Hsp_query-to>4155</Hsp_query-to>
              <Hsp_hit-from>150</Hsp_hit-from>
              <Hsp_hit-to>1012</Hsp_hit-to>
              <Hsp_identity>843</Hsp_identity>
              <Hsp_align-len>861</Hsp_align-len>
              <Hsp_qseq>AC--GT</Hsp_qseq>
              <Hsp_hseq>ACGGGT</Hsp_hseq>
              <Hsp_midline>||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>"#;

    #[test]
    fn parses_hits_and_hsps() {
        let r = parse_xml(SAMPLE, "R1TEST").unwrap();
        assert_eq!(r.program, "blastn");
        assert_eq!(r.db, "nr");
        assert_eq!(r.query_length, 4361);
        assert_eq!(r.rid, "R1TEST");
        assert_eq!(r.hits.len(), 2);
        let h1 = &r.hits[0];
        assert_eq!(h1.accession, "U00004");
        assert_eq!(h1.identity, 4361.0);
        assert_eq!(h1.length, 4361);
        assert_eq!(h1.query_from, 1);
        assert_eq!(h1.hit_to, 4361);
        let h2 = &r.hits[1];
        assert_eq!(h2.accession, "JX128753");
        assert_eq!(h2.evalue, "1e-180");
        assert_eq!(h2.length, 861);
        assert_eq!(h2.hseq, "ACGGGT");
    }

    #[test]
    fn urlencode_keeps_nucleotides() {
        assert_eq!(urlencode("ACGTacgt012"), "ACGTacgt012");
        assert_eq!(urlencode("a b&c"), "a%20b%26c");
    }
}
