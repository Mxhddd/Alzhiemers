# Sources for data/extracted/meta_extraction.csv and meta_validation_moraros2017.csv

Extracted by Claude, 2026-10-13. Every number was copied verbatim from the open full text
(Europe PMC fullTextXML or NCBI efetch db=pmc; supplementary files via Europe PMC).
Each row of the CSV carries the locator and a verbatim snippet.
**verified: false** for every row until a human checks it against the paper.

How these values were obtained:

- The candidate estimates (all of them, including those not used) are in `candidates/`.
- The selection and exclusion rules are in `src/metaanalysis/build_extraction.py` and `prereg/metaanalysis.md` §3.
- The SHA-256 of every source file used is in `candidates/source_files.sha256`.
- The source files themselves are not redistributed: several licences are NC/ND, and author manuscripts are not CC-licensed.

**Licence / access.**

- **OA** means Europe PMC open-access subset, under the stated Creative Commons licence.
- **NIH author manuscript** means the full text is free to read and to text-mine under the funder policy, but it is not CC-licensed. Europe PMC marks these isOpenAccess = N.

| source_id (PMCID) | rows | citation | licence / access |
|---|---|---|---|
| PMC9976965 | 1 | Beaman EE, Bonde AN, Larsen SMU, Ozenne B, Lohela TJ, Nedergaard M, Gíslason GH, Knudsen GM, Holst SC. Blood-brain barrier permeable β-blockers linked to lower risk of Alzheimer's disease in hypertension. Brain : a journal of neurology 2023;146(3):1141-1151. doi:10.1093/brain/awac076 | OA, CC BY-NC |
| PMC13127832 | 1 | Desai RJ, Varma VR, Mahesri M, Lee SB, Freedman A, Gerhard T, Segal J, Vine S, Ritchey MBE, Horton DB, Thambisetty M. Population-Based Validation Results From the Drug Repurposing for Effective Alzheimer's Medicines (DREAM) Study. Clinical pharmacology and therapeutics 2025;117(4):1039-1050. doi:10.1002/cpt.3583 | NIH author manuscript (free full text, text-mining permitted; not CC-licensed) |
| PMC12693173 | 1 | Chen SS, Chang YC, Yu CY, Hsieh TY, Sung WW. Differential Risks of Dementia, Depression, and Injury Among Common α-Blockers, with Tamsulosin as the Reference Drug: A Real-World Cohort Study in Men with Lower Urinary Tract Symptoms. Journal of clinical medicine 2025;14(23):8302. doi:10.3390/jcm14238302 | OA, CC BY |
| PMC11340942 | 1 | Fung KW, Baye F, Baik SH, McDonald CJ. Tamsulosin use in benign prostatic hyperplasia and risks of Parkinson's disease, Alzheimer's disease and mortality: An observational cohort study of elderly Medicare enrollees. PloS one 2024;19(8):e0309222. doi:10.1371/journal.pone.0309222 | OA, CC0 |
| PMC13644072 | 1 | Zhang Q, Schultz JL, Simmering JE, Kirkpatrick BQ, Weber MA, Skuodas S, Hicks T, Pierce G, Laughlin M, George B, Bertolli AX, Larson T, Thangavel R, Oya M, Meyerholz DK, Aldridge G, Fassler J, Narayanan NS. Glycolysis-enhancing α&lt;sub&gt;1&lt;/sub&gt;-adrenergic antagonists have therapeutic potential in Alzheimer's disease. Alzheimer's research & therapy 2026;():. doi:10.1186/s13195-026-02135-z | OA, CC BY |
| PMC9981768 | 1 | Weber MA, Sivakumar K, Tabakovic EE, Oya M, Aldridge GM, Zhang Q, Simmering JE, Narayanan NS. Glycolysis-enhancing α<sub>1</sub>-adrenergic antagonists modify cognitive symptoms related to Parkinson's disease. NPJ Parkinson's disease 2023;9(1):32. doi:10.1038/s41531-023-00477-1 | OA, CC BY |
| PMC9542191 | 1 | Latvala L, Tiihonen M, Murtola TJ, Hartikainen S, Tolppanen AM. Use of α1-adrenoceptor antagonists tamsulosin and alfuzosin and the risk of Alzheimer's disease. Pharmacoepidemiology and drug safety 2022;31(10):1110-1120. doi:10.1002/pds.5503 | OA, CC BY-NC |
| PMC9857015 | 1 | Garcia-Argibay M, Hiyoshi A, Fall K, Montgomery S. Association of 5α-Reductase Inhibitors With Dementia, Depression, and Suicide. JAMA network open 2022;5(12):e2248135. doi:10.1001/jamanetworkopen.2022.48135 | OA, CC BY |
| PMC9708561 | 1 | Colbourne L, Harrison PJ. Brain-penetrant calcium channel blockers are associated with a reduced incidence of neuropsychiatric disorders. Molecular psychiatry 2022;27(9):3904-3912. doi:10.1038/s41380-022-01615-6 | OA, CC BY |
| PMC13572180 | 1 | Wu Q, Li L, Lei Y, Zhou T, Tang H, Zhang B, Lu Y, Salvatore M, Zhang D, Tong J, Wang T, Chen S, Li H, Xu Z, Huang Y, Hu Y, Hong N, Zhou Y, Lin F, O'Brien KS, Chen Y, Li R, Xu H, Wang F, Bian J, Wolk DA, Chen Y. Discovering repurposable drugs for Alzheimer's disease and related dementias: target trial emulation using decentralised real-world data. EBioMedicine 2026;132():106466. doi:10.1016/j.ebiom.2026.106466 | OA, CC BY-NC-ND |
| PMC8751379 | 1 | Fang J, Zhang P, Wang Q, Chiang CW, Zhou Y, Hou Y, Xu J, Chen R, Zhang B, Lewis SJ, Leverenz JB, Pieper AA, Li B, Li L, Cummings J, Cheng F. Artificial intelligence framework identifies candidate targets for drug repurposing in Alzheimer's disease. Alzheimer's research & therapy 2022;14(1):7. doi:10.1186/s13195-021-00951-z | OA, CC BY |
| PMC10191296 | 1 | Song Y, Racette BA, Camacho-Soto A, Searles Nielsen S. Biologic targets of prescription medications and risk of neurodegenerative disease in United States Medicare beneficiaries. PloS one 2023;18(5):e0285011. doi:10.1371/journal.pone.0285011 | OA, CC BY |
| PMC6909196 | 1 | Kern DM, Cepeda MS, Lovestone S, Seabrook GR. Aiding the discovery of new treatments for dementia by uncovering unknown benefits of existing medications. Alzheimer's & dementia (New York, N. Y.) 2019;5():862-870. doi:10.1016/j.trci.2019.07.012 | OA, CC BY-NC-ND |
| PMC5383251 | 1 | Then CK, Chi NF, Chung KH, Kuo L, Liu KH, Hu CJ, Shen SC, Lin YK. Risk analysis of use of different classes of antidepressants on subsequent dementia: A nationwide cohort study in Taiwan. PloS one 2017;12(4):e0175187. doi:10.1371/journal.pone.0175187 | OA, CC BY |
| PMC11772729 | 1 | Ding M, Wennberg AM, Engström G, Modig K. Use of common cardiovascular disease drugs and risk of dementia: A case-control study in Swedish national register data. Alzheimer's & dementia : the journal of the Alzheimer's Association 2025;21(1):e14389. doi:10.1002/alz.14389 | OA, CC BY-NC |
| PMC6690390 | 2 | Wändell P, Carlsson AC, Sundquist J, Sundquist K. Antihypertensive drugs and relevant cardiovascular pharmacotherapies and the risk of incident dementia in patients with atrial fibrillation. International journal of cardiology 2018;272():149-154. doi:10.1016/j.ijcard.2018.07.106 | NIH author manuscript (free full text, text-mining permitted; not CC-licensed) |
| PMC12746425 | 1 | Tonegawa-Kuji R, Karavani E, Danziger M, Zhang P, Hou Y, Zhou Y, Bykova M, Pieper AA, Rosen-Zvi M, Cummings J, Cheng F. Critical evaluation of real-world evidence of repurposable medicines in the Alzheimer's disease drug development pipeline using a target trial emulation. Alzheimer's & dementia (New York, N. Y.) 2026;12(1):e70193. doi:10.1002/trc2.70193 | OA, CC BY-NC |
| PMC5880772 | 3 | Murray MD, Hendrie HC, Lane KA, Zheng M, Ambuehl R, Li S, Unverzagt FW, Callahan CM, Gao S. Antihypertensive Medication and Dementia Risk in Older Adult African Americans with Hypertension: A Prospective Cohort Study. Journal of general internal medicine 2018;33(4):455-462. doi:10.1007/s11606-017-4281-x | OA, CC BY |
| PMC4358759 | 1 | Gray SL, Anderson ML, Dublin S, Hanlon JT, Hubbard R, Walker R, Yu O, Crane PK, Larson EB. Cumulative use of strong anticholinergics and incident dementia: a prospective cohort study. JAMA internal medicine 2015;175(3):401-407. doi:10.1001/jamainternmed.2014.7663 | NIH author manuscript (free full text, text-mining permitted; not CC-licensed) |
| PMC5915701 | 1 | Richardson K, Fox C, Maidment I, Steel N, Loke YK, Arthur A, Myint PK, Grossi CM, Mattishent K, Bennett K, Campbell NL, Boustani M, Robinson L, Brayne C, Matthews FE, Savva GM. Anticholinergic drugs and risk of dementia: case-control study. BMJ (Clinical research ed.) 2018;361():k1315. doi:10.1136/bmj.k1315 | OA, CC BY-NC |
| PMC6593623 | 1 | Coupland CAC, Hill T, Dening T, Morriss R, Moore M, Hippisley-Cox J. Anticholinergic Drug Exposure and the Risk of Dementia: A Nested Case-Control Study. JAMA internal medicine 2019;179(8):1084-1093. doi:10.1001/jamainternmed.2019.0677 | OA, CC BY |
| PMC6802337 | 1 | Grossi CM, Richardson K, Fox C, Maidment I, Steel N, Loke YK, Arthur A, Myint PK, Campbell N, Boustani M, Robinson L, Brayne C, Matthews FE, Savva GM. Anticholinergic and benzodiazepine medication use and risk of incident dementia: a UK cohort study. BMC geriatrics 2019;19(1):276. doi:10.1186/s12877-019-1280-2 | OA, CC BY |
| PMC7391507 | 1 | Suh Y, Ah YM, Han E, Jun K, Hwang S, Choi KH, Lee JY. Dose response relationship of cumulative anticholinergic exposure with incident dementia: validation study of Korean anticholinergic burden scale. BMC geriatrics 2020;20(1):265. doi:10.1186/s12877-020-01671-z | OA, CC BY |
| PMC6494800 | 1 | Joung KI, Kim S, Cho YH, Cho SI. Association of Anticholinergic Use with Incidence of Alzheimer's Disease: Population-based Cohort Study. Scientific reports 2019;9(1):6802. doi:10.1038/s41598-019-43066-0 | OA, CC BY |
| PMC5383287 | 1 | Yang YW, Liu HH, Lin TH, Chuang HY, Hsieh T. Association between different anticholinergic drugs and subsequent dementia risk in patients with diabetes mellitus. PloS one 2017;12(4):e0175335. doi:10.1371/journal.pone.0175335 | OA, CC BY |
| PMC10311860 | 1 | Poonawalla IB, Xu Y, Gaddy R, James A, Ruble M, Burns S, Dixon SW, Suehs BT. Anticholinergic exposure and its association with dementia/Alzheimer's disease and mortality in older adults. BMC geriatrics 2023;23(1):401. doi:10.1186/s12877-023-04095-7 | OA, CC BY |
| PMC12538848 | 1 | Zhu N, Eriksdotter M, Liabeuf S, Tan ECK, Johnell K, Garcia-Ptacek S, Xu H. Anticholinergic burden and incident dementia: a Swedish nationwide case-control study. Alzheimer's research & therapy 2025;17(1):227. doi:10.1186/s13195-025-01883-8 | OA, CC BY |
| PMC11580265 | 1 | Iyen B, Coupland C, Bell BG, Ashcroft DM, Orrell MW, Bishara D, Dening T, Avery AJ. Risk of dementia associated with anticholinergic drugs for overactive bladder in adults aged ≥55 years: nested case-control study. BMJ medicine 2024;3(1):e000799. doi:10.1136/bmjmed-2023-000799 | OA, CC BY |
| PMC12695687 | 1 | Chi D, Liang X, Xiao Z, Zhou X, Zhao Q, Wang B, Ding D. Joint Effects of Anticholinergic Burden and Neurofilament Light on Dementia Risk: The Shanghai Aging Study. CNS neuroscience & therapeutics 2025;31(12):e70691. doi:10.1002/cns.70691 | OA, CC BY |
| PMC9974875 | 1 | Gildengers A, Gildengers A, Stoehr GP, Ran X, Jacobsen E, Teverovsky E, Chang CH, Ganguli M. Anticholinergic Drug Burden and Risk of Incident MCI and Dementia: A Population-based Study. Alzheimer disease and associated disorders 2023;37(1):20-27. doi:10.1097/wad.0000000000000538 | NIH author manuscript (free full text, text-mining permitted; not CC-licensed) |
| PMC7033580 | 1 | Liu YP, Chien WC, Chung CH, Chang HA, Kao YC, Tzeng NS. Are Anticholinergic Medications Associated With Increased Risk of Dementia and Behavioral and Psychological Symptoms of Dementia? A Nationwide 15-Year Follow-Up Cohort Study in Taiwan. Frontiers in pharmacology 2020;11():30. doi:10.3389/fphar.2020.00030 | OA, CC BY |
| PMC12022375 | 1 | Pham Nguyen TP, Thibault D, Gray SL, Weintraub D, Willis AW. Impact of Anticholinergic Burden and Clinical-Demographic Characteristics on Incident Dementia in Parkinson Disease. Journal of geriatric psychiatry and neurology 2025;38(4):254-263. doi:10.1177/08919887241313376 | OA, CC BY-NC |
| PMC13570850 | 1 | Amrein MA, Hengartner MP, Näpflin M, Farcher R, Huber CA. Long-Term Antidepressant Use and Dementia Risk: A Propensity Score-Matched Study of Swiss Older Adults, with an Anticholinergic-Specific Analysis. Drugs - real world outcomes 2026;13(3):329-338. doi:10.1007/s40801-026-00560-y | OA, CC BY-NC |
| PMC12274944 | 1 | Li B, Reid J, McClure A, Clemens KK, Welk B. The Impact of Real-world Use of Overactive Bladder Medications on Dementia Risk in Younger Adults. European urology open science 2025;78():32-40. doi:10.1016/j.euros.2025.06.008 | OA, CC BY-NC-ND |
| PMC5347943 | 6 | Moraros J, Nwankwo C, Patten SB, Mousseau DD. The association of antidepressant drug usage with cognitive impairment or dementia, including Alzheimer disease: A systematic review and meta-analysis. Depression and anxiety 2017;34(3):217-226. doi:10.1002/da.22584 | OA, CC BY |

## Searched, nothing usable found (open access)

| Contrast or class | What was searched | Result |
|---|---|---|
| atomoxetine vs stimulants | Europe PMC, OA filter; full query in `candidates/queries_2026-10-13.txt` | No primary study with an incident-dementia outcome |
| guanfacine (alone) | Europe PMC, OA filter | None; the only central-alpha2 estimate is Murray 2018's class row |
| disulfiram | Europe PMC, OA filter | None |
| isradipine or other brain-penetrant DHP | Europe PMC, OA filter | Only Colbourne 2022 (the active-comparator row) |
| Z-drugs vs DORAs | Europe PMC, OA filter | Only Tonegawa-Kuji 2026 (supplementary Table S14). Gao 2026 (PMC13008954) is a pharmacovigilance disproportionality study, not incidence, and was not used |

Anticholinergic searches, rejected candidates and blocked records are logged in
`candidates/meta_anticholinergic_log.md`. Six records had no retrievable full text through
the allowed routes: Lockery 2021, Carrière 2009, Taylor-Rowan 2021 (Cochrane), Ruxton 2015,
Zheng 2021 and Kotochinsky 2025. Several are paywalled (e.g. Heath 2018, Hong 2019).

## Oddities found in sources (left as printed)

- **Zhang 2026 (PMC13644072)** prints HR 0.84, 95% CI [0.84–0.84], with t = −867.9 on an IPTW pseudo-population. Excluded (E2).
- **Liu 2020 (PMC7033580)** prints CIs that are markedly asymmetric on the log scale. Excluded (E2).
- **Chi 2025 (PMC12695687)**: the abstract calls HR 1.55 (1.08–2.24) adjusted, but in Table 3 that value is Model 1 (univariate). Model 2 is used.
- **Moraros 2017 (PMC5347943)** has three internal inconsistencies:
  - Cherbuin's upper CI is 7.06 in Table 3 but 7.07 in the text. Table 3 is used.
  - Goveas is listed in both depression subgroups.
  - The text calls Goveas 2012 an RCT, while Table 3 says prospective cohort.
- **Desai 2025 (PMC13127832)** gives four analyses with equal stated priority. The first (as-treated) is used (rule S2).
- **Tonegawa-Kuji 2026** prints the HR as zolpidem vs suvorexant. The pipeline inverts it for the DORA class (`invert_ratio`). The number is stored as printed.
