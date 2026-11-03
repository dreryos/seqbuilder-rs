# SeqBuilder Pro 19 — Gap Analysis vs. SeqBuilder RS

Static reverse-engineering of the shipped SeqBuilder Pro 19 macOS app
(`/Applications/DNASTAR/Lasergene19/SeqBuilder Pro.app`, arm64 slice, 40 MB,
symbols not stripped — 1,381 C++ classes) compared against the current
SeqBuilder RS clone. Evidence: `Contents/Resources/SeqBuilder.xrc` (complete
wxWidgets UI definition: all menus and 46 dialogs), symbol table (`nm`),
binary strings, and bundled resources. No execution (license server offline).

## Key architectural finding

SeqBuilder Pro is a **wxWidgets C++ application**, not Cocoa. The entire UI
(menus, dialogs, panels) is defined in `SeqBuilder.xrc` (513 KB XML), which
serves as a complete interface specification for the clone. Internal C++
classes (see `seqbuilder-pro-classes.txt`) reveal the module structure.

## What the clone already matches

- Views: Circular, Linear, Sequence, PrimerDesign
- Tabs: Features, Primers, Comment, Minimap, Site Summary, History
- Side panels: Style/Enzymes, Settings/Details
- Enzymes: database (.ezd), selectors (.sel), Selector Manager (read-only)
- File I/O: .sbd read/write (round-trip), GenBank import/export, FASTA read
- Primer pair design with formulas taken from SeqBuilder's own code
- Undo/redo, translations (6 frames), grouping, complement strand, styles

## Missing modules (in recommended implementation order)

### 1. ORF detection (`ORFSearcher`, `ORFObject`, `ORFCriteriaDialog`, `ORFStyles`)
Sequence ▸ ORF Detection Options: minimum ORF size (AAs) and a three-way
radio (show only ORFs with start codon / without start codon / both).
ORFs have frame (1–6), strand, "ORF %d" names, a Show/Hide ORFs toggle and a
separate ORF Style editor. Source files: `ORFCriteriaDialog.cpp`,
`DisplayObjORF.cpp`, `DisplayObjLinearORF.cpp`, `LinearORFShape.cpp`,
`RowCreatorORF.cpp`, `ORFDrawUtils.cpp`.

### 2. Genetic codes (`GeneticCode`, `GeneticCodeEditor`, `CodonRegistry`, `CodonShuffler`)
Edit/select genetic code (.rtc files), translation exceptions, codon usage
tables. Opens the door to protein support.

### 3. Feature depth (`Feature*` — 100+ classes)
Auto-Annotate sequence, Feature Library (send to + manager), Join to Feature,
New Feature With Translation, Update All Translations, feature search,
Export features as text.

### 4. Tools menu (`ToolRunner`, `Tools*Panel`, `ScriptBuilder`)
Batch Clone / Reverse Complement / Translate, Concatenate, Extract Features
as Sequences, Export into single file, Sample by (5'/3' end tags, contains
sequence, probability, length, name), Detect Duplicate Features, Create BED,
Change Origin, Convert Ambiguity Characters. Scripting: `.script` templates
shipped in the app bundle (`cloneVerification.template.script`, …).

### 5. Cloning suite (`CloningCommands` 358 methods, `CloningStrategy`, `CloningTree`, `GatewayCloningData`, `GoldenGate*`, `MakeTAVector`, `MakeTOPOVector`)
- Restriction cloning: Copy/Paste Restriction Insert, PCR-directed
  restriction cloning, Change Overhangs
- Modern assembly: Gibson (LIC/SLIC/CPEC/SLiCE), In-Fusion, GeneArt,
  Golden Gate (domestication sites, overhang sets)
- Gateway: BP/LR reactions, entry/expression clones, MultiSite Gateway Pro
- TA / Directional TOPO cloning (incl. "Make a TA/TOPO Vector")
- Verify Clone (markup + map thumbnails), Batch Clone, clone preview/tree

### 6. Gel simulation (`GelView` 116 methods, `AuiGelWindow`, `GelMarker*`)
Agarose gel simulation of restriction digests with real hardware data:
`Markers.tdl` (molecular rulers) and `Gelbox.tdl` (gel box dimensions) are
plain TSV in the app bundle.

### 7. Lower priority
- Projects/collections (`CollectionWindow` 358 methods), vector/primer
  catalogs (.pri, .library), GenVision export (.gnv)
- Protein support (New Protein, AA views)
- Print/page layout (`PrintCommands`, `PageForm`, batch print, layout templates)
- Curtain View, Fragment View, Copy as Picture / formatted / spliced sequence
- RNA toggle, New Origin, Special Translation Ruler
- BLAST/Entrez online search (needs a decision on network features)
- IPC "Send Sequence To" integration with other DNASTAR apps

## File formats in the original not yet in the clone

| Format | Purpose |
|---|---|
| `.pri` | primer catalog |
| `.rtc` | genetic code |
| `.gnv` | GenVision project export |
| `.library` | fragment library |
| `.abi` / `.scf` | chromatogram traces |
| `.tab` / `.text` | tabular export |

## Evidence files

- `docs/seqbuilder-pro-classes.txt` — all 1,381 class names from the arm64
  symbol table (wxWidgets/std internals filtered out)
- Menus extracted from `SeqBuilder.xrc` are reproduced in the chat history;
  the original file remains at
  `/Applications/DNASTAR/Lasergene19/SeqBuilder Pro.app/Contents/Resources/SeqBuilder.xrc`
