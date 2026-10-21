# Module 1 extracted values — sources (2026-10-21)

Every value in `data/extracted/mr_*.csv` was extracted by Claude from the
source named in its row, with a locator and the verbatim text, and is marked
`verified: false` until the researcher checks it. Local copies of the sources
are in `data/raw/papers/mr/` (git-ignored).

| source_id | Citation | Access used |
|---|---|---|
| huang2021 | Huang W, Xiao J, Ji J, Chen L. Association of lipid-lowering drugs with COVID-19 outcomes from a Mendelian randomization study. *eLife* 2021;10:e73873. doi:10.7554/eLife.73873. PMID 34866576, PMC8709572. CC BY 4.0 | Europe PMC full-text XML and supplementary files (`elife-73873-supp1.xlsx`, sha256 b75429f0…77cc3f) |
| ovcharenko2005 | Ovcharenko I, Loots GG, Nobrega MA, Hardison RC, Miller W, Stubbs L. Evolution and functional classification of vertebrate gene deserts. *Genome Res* 2005;15:137–145. doi:10.1101/gr.3015505. PMID 15590943, PMC540279 | PMC article page (pmc.ncbi.nlm.nih.gov/articles/PMC540279/, HTTP 200, no login). The publisher's PDF answered 403 and its XML is not distributed by PMC; neither was pursued |
| sun2023 | Sun BB, Chiou J, Traylor M, et al. Plasma proteomic associations with genetics and health in the UK Biobank. *Nature* 2023;622:329–338. doi:10.1038/s41586-023-06592-6. CC BY 4.0 | nature.com article page and Supplementary Tables (`41586_2023_6592_MOESM3_ESM.xlsx`) |
| deklein2023 | de Klein N, Tsai EA, Vochteloo M, et al. Brain expression quantitative trait locus and network analyses reveal downstream effects and putative drivers for brain-related diseases. *Nat Genet* 2023;55:377–388. doi:10.1038/s41588-023-01300-6. CC BY 4.0 | nature.com article page |
| ferkingstad2021 | Ferkingstad E, Sulem P, Atlason BA, et al. Large-scale integration of the plasma proteome with genetics and disease. *Nat Genet* 2021;53:1712–1721. doi:10.1038/s41588-021-00978-w. CC BY 4.0 | nature.com article page and Supplementary Tables (`41588_2021_978_MOESM4_ESM.xlsx`) |

## HMGCR validation: candidates considered

| Candidate | Readable legally? | Reproducible from open data? | Used? |
|---|---|---|---|
| Huang et al. 2021, eLife (PMC8709572), Supplementary file 1 Table 6: IVW OR 1.44 (1.24–1.68), 7 SNPs listed in Table 3 | yes (CC BY) | yes: GLGC 2013 and CARDIoGRAMplusC4D 1000G are open, and the 7 SNPs' GLGC statistics equal the paper's Table 3 | **yes** |
| Ference et al. 2016, NEJM 375:2144 (PCSK9/HMGCR), per 10 mg/dL lower LDL-C: OR 0.81 for cardiovascular events, 1.13 for diabetes (point estimates as given by the researcher) | the article page answered 403 to automated access (Cloudflare check); not read, CIs not taken. The researcher can read it in a browser | no: individual-level data from many studies, and the paper's own composite cardiovascular outcome | no |
| PMC11223406 (J Transl Med 2024, asthma): HMGCR on CHD, IVW OR 1.444 (1.240–1.682), 7 SNPs | yes | same numbers as Huang 2021 to 3 decimals; instruments not listed, gene window not stated | no (the eLife paper lists its instruments) |
| PMC11260000 (Open Med 2024): IVW OR 1.62 (1.36–1.93) | yes | no: instruments from UK Biobank LDL-C (OpenGWAS ieu-b-110), not GLGC | no |
| Williams et al. 2020, Ann Neurol (PMC6944510): CAD positive control | yes | no: principal-components method on correlated variants; CAD estimates only in a figure | no |
| PMC10691537 | — | excluded by the researcher (its printed CI does not contain its estimate) | no |
