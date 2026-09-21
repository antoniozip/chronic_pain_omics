# Metabolomics screening (n = 39 retrieved, 9 included)

**Screened:** 2026-08-15 · **Record:** `metabolomics_manual_review.csv`
**Retrieval:** MetaboLights via EBI Search — see `plan/2026-08-15-metabolomics-retrieval-repair.md`

Every retrieved record carries a disposition, as PRISMA requires. Verdicts on
the 16 prioritised studies are the reviewer's; the 23 non-prioritised were
excluded as a block on the reviewer's instruction that none carries a pain
phenotype, and each is recorded with the criterion it fails rather than being
deleted.

| | n |
|---|---|
| retrieved | 39 |
| excluded, automatic screen (species or no pain phenotype) | 23 |
| excluded, reviewer | 7 |
| **included** | **9** |

Inclusion here means "carries a pain-versus-control design worth opening". It
does **not** mean usable: no study's files have been examined yet, so per-sample
values, group labels and quantification tables are all still unverified. The
proteomic arm lost 84% at exactly this step.

## Included (n = 9)

| # | Accession | Species | Platform | Year | Title |
|---|---|---|---|---|---|
| 1 | [MTBLS12704](https://www.ebi.ac.uk/metabolights/MTBLS12704) | Homo sapiens | mass spectrometry | 2025 | Exploring Salivary Metabolites as Biomarkers in Chronic Craniofacial and Orofacial Pain: A Metabolomic Analysis. |
| 2 | [MTBLS13513](https://www.ebi.ac.uk/metabolights/MTBLS13513) | Homo sapiens | mass spectrometry | 2025 | Multi-omics analysis identifies a microbiota–bile acid–TLR signaling axis driving bladder injury in interstitial cystitis |
| 3 | [MTBLS13869](https://www.ebi.ac.uk/metabolights/MTBLS13869) | Rattus norvegicus | mass spectrometry | 2026 | Serum Metabolomics in the Rat Model of Recurrent Pelvic Pain |
| 4 | [MTBLS14059](https://www.ebi.ac.uk/metabolights/MTBLS14059) | Mus musculus | mass spectrometry | 2026 | Semaglutide targets muscle mitochondria to regulate osteoarthritis |
| 5 | [MTBLS14295](https://www.ebi.ac.uk/metabolights/MTBLS14295) | Mus musculus | mass spectrometry | 2026 | Semaglutide targets muscle mitochondria to regulate glutamine metabolism and treat osteoarthritis |
| 6 | [MTBLS14319](https://www.ebi.ac.uk/metabolights/MTBLS14319) | Mus musculus | mass spectrometry | 2026 | In-depth characterization of the CaV2.2-knockout mouse line |
| 7 | [MTBLS2774](https://www.ebi.ac.uk/metabolights/MTBLS2774) | Homo sapiens | mass spectrometry | 2021 | Altered metabolome and microbiome features provide clues in understanding irritable bowel syndrome and depression comorbidity |
| 8 | [MTBLS5667](https://www.ebi.ac.uk/metabolights/MTBLS5667) | Rattus norvegicus | mass spectrometry | 2024 | Clinical efficacy of Yiqi Yangxue formula on knee osteoarthritis and unraveling therapeutic mechanism through plasma metabolites in rats. |
| 9 | [MTBLS9662](https://www.ebi.ac.uk/metabolights/MTBLS9662) | Homo sapiens | mass spectrometry | 2025 | Metabolic Profiling of Synovial Fluid in Human Temporomandibular Joint Osteoarthritis |

## Excluded by the reviewer (n = 7)

| # | Accession | Species | Title | Reason code | Proposed reason |
|---|---|---|---|---|---|
| 1 | [MTBLS12737](https://www.ebi.ac.uk/metabolights/MTBLS12737) | Homo sapiens | Prediction of irritable bowel syndrome by integrating urine metabolites and gut microbiota | not_recorded | no_pain_versus_control_contrast |
| 2 | [MTBLS13512](https://www.ebi.ac.uk/metabolights/MTBLS13512) | Homo sapiens | Multi-omics analysis identifies a microbiota–bile acid–TLR signaling axis driving bladder injury in interstitial cystitis | duplicate_of_MTBLS13513 | — |
| 3 | [MTBLS14200](https://www.ebi.ac.uk/metabolights/MTBLS14200) | Homo sapiens | Integrated microbiome and metabolome analysis of milk of lactational mastitis women in remission period | not_recorded | no_pain_phenotype |
| 4 | [MTBLS14816](https://www.ebi.ac.uk/metabolights/MTBLS14816) | Mus musculus | Therapeutic Effect of Modified Meridian-Guided Acupoint Pressing on Lumbar Facet Joint Osteoarthritis: An Integrated Microbiomics and Metabolomics Analysis | not_recorded | no_pain_versus_control_contrast |
| 5 | [MTBLS1894](https://www.ebi.ac.uk/metabolights/MTBLS1894) | Homo sapiens | 1H-NMR-Based Analysis for Exploring Knee Synovial Fluid Metabolite Changes after Local Cryotherapy in Knee Arthritis Patients | not_recorded | intervention_pre_post_not_case_control |
| 6 | [MTBLS6219](https://www.ebi.ac.uk/metabolights/MTBLS6219) | Homo sapiens | Diagnostic utility of clinicodemographic, biochemical and metabolite variables to identify viable pregnancies in a symptomatic cohort during early gestation. | not_recorded | no_pain_phenotype |
| 7 | [MTBLS9665](https://www.ebi.ac.uk/metabolights/MTBLS9665) | Homo sapiens | Metabolic Profiling of Synovial Fluid in Human Temporomandibular Joint Osteoarthritis | duplicate_of_MTBLS9662 | — |

Two exclusions are settled by the record itself — MTBLS13512 and MTBLS9665 are
the twins of included duplicates. The other five carry `not_recorded`: the
reviewer gave a verdict without a reason, and inventing one would put a
fabricated justification into the PRISMA record. The `Proposed reason` column
holds an inference from the title, for the reviewer to confirm or replace.

## Excluded by automatic screen (n = 23)

Species outside {Homo sapiens, Mus musculus, Rattus norvegicus}, or no pain
term in title or description.

| # | Accession | Species | Title | Reason code |
|---|---|---|---|---|
| 1 | [MTBLS11752](https://www.ebi.ac.uk/metabolights/MTBLS11752) | Mus musculus | Deciphering the Role of N-Acetyltransferase 10 in Thalamic Hemorrhage Through Integrative Multi-Omics and Experimental Validation | no_pain_phenotype |
| 2 | [MTBLS11941](https://www.ebi.ac.uk/metabolights/MTBLS11941) | Mus musculus | The effect of an NAT10 inhibitor (remodelin) on inflammation-related genes in the treatment of thalamic hemorrhage based on metabolic data | no_pain_phenotype |
| 3 | [MTBLS12716](https://www.ebi.ac.uk/metabolights/MTBLS12716) | Rattus norvegicus | SADI-S and SG Surgeries Induce Distinct Bile Acid Profiles Linked to Improved Glucose Metabolism via Microbiota Interactions_serum | no_pain_phenotype |
| 4 | [MTBLS12748](https://www.ebi.ac.uk/metabolights/MTBLS12748) | Ovis aries | Saline pasture improve meat quality in Qinghai Tibetan sheep through changes in the rumen microbiota | wrong_species |
| 5 | [MTBLS12829](https://www.ebi.ac.uk/metabolights/MTBLS12829) | Mus musculus | Zingerone treats postmenopausal osteoporosis via increased ferroptosis sensitivity by p53-mediated regulation of Sat1 and Gpx4 expression | no_pain_phenotype |
| 6 | [MTBLS13219](https://www.ebi.ac.uk/metabolights/MTBLS13219) | Mus musculus | The integration of network pharmacology and multiomics reveals the mechanism by which Malus hupehensis leaves inhibit acute liver injury | no_pain_phenotype |
| 7 | [MTBLS1991](https://www.ebi.ac.uk/metabolights/MTBLS1991) | Homo sapiens | A molecular index for biological age identified from the metabolome and senescence-associated secretome in humans. | no_pain_phenotype |
| 8 | [MTBLS2040](https://www.ebi.ac.uk/metabolights/MTBLS2040) | Homo sapiens | Symptomatology and Serum Nuclear Magnetic Resonance Metabolomics; Do They Predict Endometriosis in Fertile Women Undergoing Laparoscopic Sterilisation? A Prospective Cross-sectional Study | no_pain_phenotype |
| 9 | [MTBLS2060](https://www.ebi.ac.uk/metabolights/MTBLS2060) | Drosophila melanogaster | Prediction of complex phenotypes using the Drosophila melanogaster metabolome | wrong_species |
| 10 | [MTBLS2689](https://www.ebi.ac.uk/metabolights/MTBLS2689) | Canis lupus familiaris | Multi-Omics Approach to Elucidate Cerebrospinal Fluid Changes in Dogs with Intervertebral Disc Herniation | wrong_species |
| 11 | [MTBLS3342](https://www.ebi.ac.uk/metabolights/MTBLS3342) | Homo sapiens | Oxidative phosphorylation is a metabolic vulnerability of endocrine therapy and palbociclib resistant metastatic breast cancers | no_pain_phenotype |
| 12 | [MTBLS354](https://www.ebi.ac.uk/metabolights/MTBLS354) | Homo sapiens | Lipid metabolites as potential diagnostic and prognostic biomarkers for acute community acquired pneumonia | no_pain_phenotype |
| 13 | [MTBLS3581](https://www.ebi.ac.uk/metabolights/MTBLS3581) | Mus musculus; Standard | Gut microbiota carcinogen metabolism causes distal tissue tumours | no_pain_phenotype |
| 14 | [MTBLS543](https://www.ebi.ac.uk/metabolights/MTBLS543) | Equus caballus | Synovial Fluid Metabolites Differentiate between Septic and Nonseptic Joint Pathologies | wrong_species |
| 15 | [MTBLS579](https://www.ebi.ac.uk/metabolights/MTBLS579) | Homo sapiens | Diagnostic metabolite biomarkers of chronic typhoid carriage | no_pain_phenotype |
| 16 | [MTBLS6038](https://www.ebi.ac.uk/metabolights/MTBLS6038) | Homo sapiens | Serum organic acid metabolites can be used as potential biomarkers to identify prostatitis, benign prostatic hyperplasia, and prostate cancer (Targeted organic acid assay) | no_pain_phenotype |
| 17 | [MTBLS6039](https://www.ebi.ac.uk/metabolights/MTBLS6039) | Homo sapiens | Serum organic acid metabolites can be used as potential biomarkers to identify prostatitis, benign prostatic hyperplasia, and prostate cancer (Untargeted assay) | no_pain_phenotype |
| 18 | [MTBLS6053](https://www.ebi.ac.uk/metabolights/MTBLS6053) | Canis lupus familiaris | Metabolomic profiling of cerebrospinal fluid from dogs with meningoencephalitis of unknown origin by nuclear magnetic resonance | wrong_species |
| 19 | [MTBLS723](https://www.ebi.ac.uk/metabolights/MTBLS723) | Mus musculus | CXCL12 and MYC control energy metabolism to support adaptive responses after kidney injury | no_pain_phenotype |
| 20 | [MTBLS7759](https://www.ebi.ac.uk/metabolights/MTBLS7759) | Mus musculus | Human umbilical cord-derived mesenchymal stem cells ameliorate perioperative neurocognitive disorder by inhibiting inflammatory responses and activating BDNF/TrkB/CREB signaling pathway in aged mice | no_pain_phenotype |
| 21 | [MTBLS779](https://www.ebi.ac.uk/metabolights/MTBLS779) | Paris polyphylla; Paris fargesii | Comparative analysis of proteomic and metabolomic profiles of different species of Paris | wrong_species |
| 22 | [MTBLS8621](https://www.ebi.ac.uk/metabolights/MTBLS8621) | Homo sapiens | UHPLC-OE-MS metabolic profile of follicular fluid in patients with ovarian endometriosis: a pilot study of an endometriosis cohort | no_pain_phenotype |
| 23 | [MTBLS9119](https://www.ebi.ac.uk/metabolights/MTBLS9119) | Homo sapiens | Serum metabolomics analysis of patients with chronic obstructive pulmonary disease and 'frequent exacerbator' phenotype | no_pain_phenotype |

## Resolved: MTBLS14059 and MTBLS14295 are distinct experiments

Both were included, and the titles, dates and shared abstract opening made them
look like one experiment deposited twice — the pseudo-replication
`conf/analysis/superseded_studies.csv` exists to prevent. Their ISA-Tab settles
it: **they are two experiments from one publication, and both may be kept.**

| | MTBLS14059 | MTBLS14295 |
|---|---|---|
| tissue | muscle | plasma |
| biological samples | 48 (+7 QC) | 16 (+5 QC) |
| design | 2×2×2 factorial | four arms |
| factors | HFD × OA × Semaglutide | OA, S-mito, C-mito |
| arms | CAS CAV CSS CSV HAS HAV HSS HSV, n=6 | CTRL OA SEMA MITO, n=4 |
| submitted | 2026-03-16 | 2026-04-16 |

The designs are structurally incompatible, which is what makes this decidable
rather than a judgement call. MTBLS14295 has mitochondria-transplantation arms
(`S-mito`, `C-mito`) with no counterpart anywhere in MTBLS14059's factorial,
and MTBLS14059 has a diet factor (`HFD`) absent from MTBLS14295 — so each
contains animals that cannot appear in the other. Group sizes differ (6 against
4), and sample names are disjoint apart from the QC injections, which are
pooled apparatus rather than biology.

Neither is registered in `superseded_studies.csv`: nothing is superseded.

**Two caveats to carry forward.** They share a laboratory, a protocol and a
publication (submitter Tian, released 2026-05-06), so they are correlated in a
way a random-effects model treating them as independent will understate — the
same concern as any multi-cohort single-lab contribution, and a reason to check
their influence in a leave-one-out. And they measure different tissues, which
in this project's transcriptomic arm was the factor that dominated variance, so
tissue belongs in the record as a moderator rather than being averaged away.

## Independence of the included set

Setting that pair aside, the 9 included studies are 8 independent experiments at
most. Species split: 4 human,
3 mouse,
2 rat.

With `meta.min_studies = 3`, a pooled metabolomic estimate needs three studies
measuring the same metabolite. Whether that holds depends on identifier overlap
(ChEBI and HMDB spaces both appear in this corpus) and is unknown until the
tables are opened — the same obstacle that made the proteomic arm a
peptide-versus-gene mapping problem before it was a statistics problem.
