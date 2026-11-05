// Genetic codes, matching SeqBuilder Pro's Sequence > Select Genetic Code and
// Edit Genetic Code (GeneticCode/CodonRegistry/GeneticCodeEditor in the
// original; the binary reads GenBank /transl_table qualifiers, i.e. the NCBI
// tables). Codon order in `aas`/`starts` is TCAG.

/// One NCBI translation table as shipped.
pub struct NcbiCode {
    pub id: u8,           // NCBI transl_table id
    pub name: &'static str,
    pub aas: [u8; 64],    // amino acid per codon, TCAG order
    pub starts: [u8; 64], // b'1' marks an initiation codon
}

/// The working copy the app edits (Genetic Code editor mutates it).
#[derive(Clone)]
pub struct GeneticCode {
    pub id: u8,           // NCBI id, 0 once user-edited
    pub name: String,
    aas: Vec<u8>,
    starts: Vec<bool>,
}

/// Codon index in TCAG order: first base T/C/A/G is the high two bits.
pub fn codon_index(codon: &[u8]) -> Option<usize> {
    let n = |b: u8| match b.to_ascii_uppercase() { b'T' | b'U' => Some(0), b'C' => Some(1), b'A' => Some(2), b'G' => Some(3), _ => None };
    let (a, b, c) = (n(codon[0])?, n(codon[1])?, n(codon[2])?);
    Some(a * 16 + b * 4 + c)
}

/// The codon (as three letters) of a TCAG-order index.
pub fn index_codon(i: usize) -> [u8; 3] {
    const B: [u8; 4] = *b"TCAG";
    [B[i / 16], B[i / 4 % 4], B[i % 4]]
}

impl GeneticCode {
    pub fn standard() -> Self { Self::from_ncbi(1).unwrap() }

    pub fn ncbi(id: u8) -> Option<&'static NcbiCode> {
        NCBI_CODES.iter().find(|c| c.id == id)
    }

    pub fn from_ncbi(id: u8) -> Option<Self> {
        let c = Self::ncbi(id)?;
        Some(Self { id: c.id, name: c.name.to_string(), aas: c.aas.to_vec(), starts: c.starts.iter().map(|&b| b == b'1').collect() })
    }

    pub fn translate(&self, codon: &[u8]) -> char {
        match codon_index(codon) {
            Some(i) => self.aas[i] as char,
            None => 'X',
        }
    }

    pub fn is_stop(&self, codon: &[u8]) -> bool {
        codon_index(codon).is_some_and(|i| self.aas[i] == b'*')
    }

    pub fn is_start(&self, codon: &[u8]) -> bool {
        codon_index(codon).is_some_and(|i| self.starts[i])
    }

    /// The one-letter amino acid of a codon index.
    pub fn aa_at(&self, i: usize) -> char { self.aas[i] as char }

    /// Set the amino acid of a codon index (Genetic Code editor).
    pub fn set_aa(&mut self, i: usize, aa: char) { if i < 64 { self.aas[i] = aa as u8; } }

    /// Toggle a codon index as an initiation codon (Edit Start Codon(s)).
    pub fn toggle_start(&mut self, i: usize) { if i < 64 { self.starts[i] = !self.starts[i]; } }

    pub fn is_start_index(&self, i: usize) -> bool { self.starts[i] }
}

impl Default for GeneticCode {
    fn default() -> Self { Self::standard() }
}

// NCBI genetic codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi,
// last update Sep 23, 2024). `aas` holds the amino acid of each codon in TCAG order;
// `starts` marks initiation codons ('1'), which translate as Met when used to start.
pub const NCBI_CODES: &[NcbiCode] = &[
    NcbiCode { id: 1, name: "The Standard Code", aas: *b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0001000000000000000100000000000000010000000000000000000000000000" },
    NcbiCode { id: 2, name: "The Vertebrate Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000011110000000000000001000000000000" },
    NcbiCode { id: 3, name: "The Yeast Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000110000000000000001000000000000" },
    NcbiCode { id: 4, name: "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0011000000000000000100000000000011110000000000000001000000000000" },
    NcbiCode { id: 5, name: "The Invertebrate Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", starts: *b"0001000000000000000000000000000011110000000000000001000000000000" },
    NcbiCode { id: 6, name: "The Ciliate, Dasycladacean and Hexamita Nuclear Code", aas: *b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 9, name: "The Echinoderm and Flatworm Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000001000000000000" },
    NcbiCode { id: 10, name: "The Euplotid Nuclear Code", aas: *b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 11, name: "The Bacterial, Archaeal and Plant Plastid Code", aas: *b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0001000000000000000100000000000011110000000000000001000000000000" },
    NcbiCode { id: 12, name: "The Alternative Yeast Nuclear Code", aas: *b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000100000000000000010000000000000000000000000000" },
    NcbiCode { id: 13, name: "The Ascidian Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", starts: *b"0001000000000000000000000000000000110000000000000001000000000000" },
    NcbiCode { id: 14, name: "The Alternative Flatworm Mitochondrial Code", aas: *b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 15, name: "Blepharisma Nuclear Code", aas: *b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 16, name: "Chlorophycean Mitochondrial Code", aas: *b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 21, name: "Trematode Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000001000000000000" },
    NcbiCode { id: 22, name: "Scenedesmus obliquus Mitochondrial Code", aas: *b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 23, name: "Thraustochytrium Mitochondrial Code", aas: *b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000010010000000000000001000000000000" },
    NcbiCode { id: 24, name: "Rhabdopleuridae Mitochondrial Code", aas: *b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", starts: *b"0001000000000000000100000000000000010000000000000001000000000000" },
    NcbiCode { id: 25, name: "Candidate Division SR1 and Gracilibacteria Code", aas: *b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0001000000000000000000000000000000010000000000000001000000000000" },
    NcbiCode { id: 26, name: "Pachysolen tannophilus Nuclear Code", aas: *b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000100000000000000010000000000000000000000000000" },
    NcbiCode { id: 27, name: "Karyorelict Nuclear Code", aas: *b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 28, name: "Condylostoma Nuclear Code", aas: *b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 29, name: "Mesodinium Nuclear Code", aas: *b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 30, name: "Peritrich Nuclear Code", aas: *b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 31, name: "Blastocrithidia Nuclear Code", aas: *b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0000000000000000000000000000000000010000000000000000000000000000" },
    NcbiCode { id: 32, name: "Balanophoraceae Plastid Code", aas: *b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", starts: *b"0001000000000000000100000000000011110000000000000001000000000000" },
    NcbiCode { id: 33, name: "Cephalodiscidae Mitochondrial UAA-Tyr Code", aas: *b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", starts: *b"0001000000000000000100000000000000010000000000000001000000000000" },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_code_translates() {
        let c = GeneticCode::standard();
        assert_eq!(c.translate(b"ATG"), 'M');
        assert_eq!(c.translate(b"TGG"), 'W');
        assert!(c.is_stop(b"TAA"));
        assert!(c.is_start(b"ATG"));
        assert!(!c.is_start(b"ATT"));
    }

    #[test]
    fn all_tables_have_64_and_valid_aas() {
        for c in NCBI_CODES {
            assert!(c.aas.iter().all(|&a| (a as char).is_ascii_uppercase() || a == b'*'));
            assert!(c.starts.iter().all(|&s| s == b'0' || s == b'1'));
            assert!(c.starts.iter().any(|&s| s == b'1'), "table {} has no starts", c.id);
        }
    }

    #[test]
    fn bacterial_code_11_differs() {
        // NCBI table 1 initiates from ATG, TTG and CTG; table 11 adds more
        let c = GeneticCode::from_ncbi(11).unwrap();
        let s = GeneticCode::standard();
        assert!(s.is_start(b"TTG") && s.is_start(b"CTG") && s.is_start(b"ATG"));
        assert!(c.is_start(b"ATT") && !s.is_start(b"ATT"));
        assert_eq!(c.translate(b"TTG"), 'L');
    }

    #[test]
    fn codon_index_roundtrip() {
        for i in 0..64 {
            let codon = index_codon(i);
            assert_eq!(codon_index(&codon), Some(i), "{codon:?}");
        }
    }

    #[test]
    fn editor_edits_apply() {
        let mut c = GeneticCode::standard();
        c.set_aa(codon_index(b"TGA").unwrap(), 'W'); // UGA -> Trp like code 2
        assert_eq!(c.translate(b"TGA"), 'W');
        assert!(!c.is_stop(b"TGA"));
        c.toggle_start(codon_index(b"GTG").unwrap());
        assert!(c.is_start(b"GTG"));
    }
}
