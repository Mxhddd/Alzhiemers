# Module 1 — positive controls per target (decided 2026-10-13)

The positive control is a hard gate. The target's instrument must move a
trait that is mechanistically downstream of the target, in the expected
direction, before the AD outcome is read. A target without a defensible
control is **excluded and reported as untestable with current
instruments, not as a null result**.

**Evidence.** GWAS Catalog REST API v2 (EBI), queried 2026-10-13 with
`mapped_gene=<GENE>`, keeping associations with p ≤ 5×10⁻⁸. The raw
summary is in `mr_gwascatalog_targets_2026-10-13.json`. Catalog mapping
is coarse, so these are reasons to expect a control to work, not proof.
The gate itself decides at run time.

| Target | Control | Higher expression should... | Why | Genome-wide evidence at the locus |
|---|---|---|---|---|
| ADRB1 | resting heart rate | raise HR | β1 blockade lowers HR (brief) | ADRB1: SBP/DBP min p 2e-63 and ECG traits; heart-rate GWAS exist |
| SLC6A2 | resting heart rate | lower HR | NET inhibition (atomoxetine) raises HR | **none for HR**: the control may fail and exclude the target |
| ADRA1A | diastolic BP | raise DBP | α1 blockade lowers BP | DBP p 7e-9 (Keaton 2024, GCST90310295); 1 association only |
| ADRA2A | fasting glucose | raise glucose | α2A over-expression impairs insulin secretion (Rosengren 2010, PMID 19965390) | glucose ×12, min p 7e-39; T2D p 2e-15; HOMA-B, insulin |
| DBH | systolic BP | raise SBP | DBH makes NE; less NE lowers BP | SBP min p 8e-21, DBP 2e-35; DBH protein p 3e-295 |
| CACNA1D | systolic BP | raise SBP | Ca-channel blockade lowers BP | SBP min p 2e-38; DBP 1e-32; heart rate 1e-23 |
| HCRTR2 | sleep duration | shorten sleep | orexin antagonists lengthen sleep | sleep duration p 6e-11 (Dashti 2019); chronotype 3e-48 (Jones 2019) |
| LGMN | — | — | no downstream trait of legumain activity | only legumain protein itself (p 2e-298) and unrelated proteins: **untestable** |
| SLC18A2 | — | — | nothing indexes VMAT2 function | 6 associations (blood proteins, ghrelin, glucose): **untestable** |
| MAOA | — | — | X-linked | 1 association, none genome-wide: **untestable** |

Using a protein's own level (e.g. legumain from a cis-pQTL) is not a
positive control: it tests that the instrument is an instrument, not that
the target's biology is captured.
