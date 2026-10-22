# Multi-cohort validation: analysis plan

Version 1.0, written 2026-10-22, committed before the consolidated
analysis was run. Changes made after this commit are listed under
[Deviations](#deviations) in `docs/VALIDATION.md`, each with its reason.

## 0. What was known when this plan was written

This plan is not a blind pre-registration, and it says so. Exploratory runs
on six cohorts (SRP075100, SRP051201, SRP226726, PRJNA863881, PRJNA786436,
PRJNA736618) were done before it, with an earlier version of the scripts.
Their results were seen, and some are summarised in the task that motivated
this plan:

* the Red complex fraction alone reached an AUC close to the SMDI;
* no index separated non-responding from responding sites at baseline
  (PRJNA786436);
* the naive Bayes and VSEARCH classifiers gave different index values.

To limit the influence of those results:

* the index set, the cut-offs, the unit of analysis and the tests below are
  the ones the published method and the previous `docs/VALIDATION.md`
  already used, or are fixed here before any new cohort was downloaded;
* the cohorts are chosen by the criteria in section 2, applied to a
  documented database search (section 2.3), not by their results;
* the comparison of the SMDI with the Red complex is a pre-specified test,
  whatever its outcome;
* the new cohorts (section 2.4) had not been downloaded when this plan was
  written.

## 1. Questions

1. **Discrimination (primary).** In cohorts independent of the SMDI's
   derivation, how well does the SMDI separate periodontally healthy
   subjects from subjects with periodontitis? Primary estimand: the pooled
   subject-level AUC of the SMDI across independent case-control cohorts.
2. **Added value over the Red complex.** Does the SMDI discriminate better
   than the Red complex fraction, the simplest published alternative?
3. **Published cut-offs.** What sensitivity and specificity do the
   published SMDI cut-offs give in new populations?
4. **Response to therapy.** Do the indices decrease after non-surgical
   periodontal therapy, and more than in untreated controls?
5. **Prognosis.** Do baseline index values predict which sites or subjects
   respond to therapy?
6. **Clinical correlation.** Do the indices correlate with probing depth,
   clinical attachment loss, bleeding on probing and stage?
7. **Robustness.** How much do the answers depend on the 16S region, the
   classifier, sequencing depth and the SMDI pseudocount?

## 2. Cohorts

### 2.1 Inclusion criteria

1. Human subgingival plaque (paper points or curette), sequenced as a 16S
   rRNA gene amplicon on a short-read platform (Illumina, Ion Torrent or
   454).
2. Raw reads in SRA/ENA/DDBJ.
3. One of these designs:
   * **case-control:** periodontally healthy subjects and subjects with
     periodontitis, at least 10 subjects per group;
   * **therapy:** subjects sampled before and after non-surgical
     (mechanical) periodontal therapy, at least 10 paired subjects or
     sites.
4. Group, subject and (for therapy) time point recoverable from public
   metadata (SRA/BioSample records, sample or library names), without
   contacting the authors.

### 2.2 Exclusion criteria

1. Population recruited for a systemic condition (diabetes, rheumatoid
   arthritis, cirrhosis, cancer, …). The disease is a confounder that the
   public metadata do not let us adjust for.
2. Interventions other than mechanical therapy (systemic antibiotics,
   probiotics, drugs, surgery), unless the arm receiving mechanical therapy
   alone is labelled.
3. Children, implants, orthodontics, or special phenotypes (molar-incisor
   pattern, Papillon-Lefèvre syndrome).
4. Long-read (PacBio, Nanopore) amplicons. The plugin targets short
   amplicons, and long reads would need another error model.
5. More than about 30 GB of raw data or more than about 12 h of
   processing on the validation machine (4 cores, 15 GB RAM).
6. Runs already covered by another included accession (duplicates).

### 2.3 Independence from the SMDI

Chen, Marsh & Al-Hebshi (2022) trained the SMDI on Abusleme et al. 2013
(SRP012422), Griffen et al. 2012 (SRP009299) and Kirst et al. 2015
(SRP051201). They tested it on Pei et al. 2020 (SRP226726), the baseline
samples of Bizzarro et al. 2016 (SRP062274) and HMP. These cohorts are
labelled **not independent**. They are analysed and reported separately,
and are never pooled with the independent cohorts.

### 2.4 Search (2026-10-22)

* **ENA portal API**, `result=study`:
  `(study_title="*periodont*" OR study_description="*periodont*") AND
  (study_title="*subgingival*" OR study_description="*subgingival*")`.
  201 studies.
* **NCBI BioProject** (E-utilities):
  `(periodontitis OR periodontal) AND (subgingival OR "dental plaque" OR
  "periodontal pocket")`. 241 projects, 47 of them not in the ENA list.
* For each project, the ENA/SRA run report gave the library strategy,
  platform and size. 97 projects had at least 20 amplicon runs. Their
  BioSample attributes were read to check the design and the labels.
* The SMDI paper's dataset table (Appendix Table 1) was added.

The records are in `validation/search/`. `screening_2026-10-22.tsv`
gives the decision and reason for each project.

| Accession | Design | Region¹ | Independent | Status before this plan |
|---|---|---|---|---|
| SRP075100 | case-control + therapy (sites) | V4 | yes | processed |
| PRJNA863881 | case-control + clinical + stage | V3-V4 | yes | processed |
| PRJNA786436 | therapy, responding vs non-responding sites | V4 | yes | processed |
| PRJNA736618 | therapy, day 0/1/7/90 | V3-V4 | yes | processed |
| PRJEB83770 | case-control (+ gingivitis) | ? | yes | new |
| PRJNA1470369 | case-control, stage I-IV | ? | yes | new |
| PRJNA1294584 | case-control, stage 2-4 | ? | yes | new |
| PRJNA773202 | case-control | V3-V4 | yes | new |
| PRJNA773195 | same samples as PRJNA773202 | V1-V3 | yes | new (region robustness only) |
| PRJNA726254 | case-control + therapy | ? | yes | new |
| PRJEB106816 | therapy vs untreated control | ? | yes | new |
| SRP051201 | case-control | V1-V3 | **no** (training) | processed |
| SRP226726 | case-control | V3-V4 | **no** (test) | processed |
| SRP009299 | case-control, V1-V2 and V4 | V1-V2, V4 | **no** (training) | failed (DADA2 error model) |
| SRP062274 | periodontitis, therapy ± antibiotics | V5-V7 | **no** (test, baseline) | new |
| SRP012422 | case-control | V1-V2 | **no** (training) | excluded: multiplexed runs without public barcodes |

¹ "?" means the region will be determined from the primer sequences found
at the start of the reads, before denoising.

## 3. Processing

The same pipeline is used for every cohort (`validation/process_dataset.sh`,
one configuration file per cohort):

1. Import, then DADA2 (`denoise-paired`, `denoise-single` or
   `denoise-pyro` for 454 and Ion Torrent). Trimming and truncation are
   chosen from the read-quality profiles and the primer positions,
   **before** any index is computed, and recorded in the configuration
   file.
2. Two classifiers, both on eHOMD 16S RefSeq V15.23:
   * **nb (primary).** Naive Bayes trained on the reference trimmed to the
     cohort's primers. This is the standard QIIME 2 workflow that plugin
     users follow.
   * **vsearch (robustness).** Consensus of the top hits at ≥ 98% identity
     against full-length eHOMD. This is the closest to the BLASTN ≥ 98%
     rule used to derive the SMDI.
3. `qiime periodontitis collapse-complexes` with default parameters.

## 4. Outcomes and measures

* **Primary index:** SMDI (49 species).
* **Secondary measures:** SMDI_simplified, PC2014_moderate and the six
  complex fractions (Red, Orange, Yellow, Green, Purple, Blue).
* **Red complex fraction:** the explicit comparator for the SMDI.
* **Unit of analysis:** the subject. Several samples or sites of a subject
  at one time point are averaged. Exceptions:
  * PRJNA786436, whose design compares sites, is analysed at site level,
    with subjects as clusters reported;
  * SRP075100 therapy, where the paired unit is the site (tooth), as
    published.
* **Case-control groups:** health vs periodontitis (any stage or
  severity). Gingivitis is reported descriptively and is not part of the
  primary contrast. In cohorts that also contain post-therapy samples,
  only baseline samples enter the case-control analysis.

## 5. Pre-specified cut-offs

Only the published cut-offs are used:

* SMDI: 0, 0.4 (best in training) and 0.5 (best in test);
* SMDI_simplified: 0.

A sample is classified as periodontitis if the index is greater than the
cut-off. Data-driven cut-offs (Youden) may be reported as descriptive only.

## 6. Statistics

All analyses are done by scripts in `validation/`. Random seeds are fixed.

### 6.1 Discrimination (Q1)

* AUC = P(periodontitis > health), with ties counted as 0.5.
* 95% CI by stratified bootstrap over subjects: 2000 resamples, cases and
  controls resampled separately, percentile interval.
* Per-cohort SE for meta-analysis: Hanley & McNeil (1982), primary;
  bootstrap SE, sensitivity.
* **Random-effects meta-analysis** of the per-cohort AUCs:
  * primary: REML estimate of τ², with a Hartung-Knapp-Sidik-Jonkman CI
    (recommended with few studies);
  * sensitivity: DerSimonian-Laird with a Wald CI;
  * reported: I², τ² and the 95% prediction interval.
* Each cohort enters the primary meta-analysis once. PRJNA773202 (V3-V4)
  represents the Santiago/Oporto cohort, because V3-V4 is the region of
  most other cohorts; PRJNA773195 (V1-V3) is used only for region
  robustness. The non-independent cohorts are pooled in a separate
  meta-analysis, which is reported but not interpreted as validation.

### 6.2 SMDI vs Red complex (Q2)

* Paired DeLong test (DeLong et al. 1988; fast algorithm of Sun & Xu
  2014) of AUC(SMDI) − AUC(Red) in each cohort.
* Random-effects meta-analysis (REML, HKSJ) of the paired differences,
  with the DeLong SE of the difference.
* The same comparison is reported for SMDI_simplified vs Red, as a
  secondary test.

### 6.3 Cut-offs (Q3)

* Sensitivity, specificity and accuracy for each cohort and cut-off,
  with Wilson 95% CIs.
* Pooled sensitivity and specificity: random-effects meta-analysis
  (DerSimonian-Laird) of logit proportions, with a 0.5 continuity
  correction for 0 or 100%.

### 6.4 Therapy (Q4)

* Paired change, follow-up minus baseline, per subject (or site, see
  section 4).
* Wilcoxon signed-rank test, median change, and the fraction that
  decreased, for each cohort and each time point.
* Effect size: matched-pairs rank-biserial correlation.
* PRJEB106816: changes in the treated arm vs the untreated control arm,
  by Mann-Whitney U, with AUC as effect size.
* SRP062274 (not independent): reported separately, by antibiotic arm.

### 6.5 Prognosis (Q5)

* AUC of the baseline index for non-response, with a bootstrap CI.
  * PRJNA786436: non-responding vs responding sites.
  * SRP062274: the responder label (not independent).

### 6.6 Clinical correlation (Q6)

* Spearman ρ of each measure with:
  * probing depth, clinical attachment loss, bleeding on probing and plaque
    index, wherever public;
  * stage, as an ordinal 0 (health) to 4 (stage IV).
* Gingivitis is excluded from the stage correlation, because it is not a
  stage.

### 6.7 Multiplicity

* The primary analysis (pooled SMDI AUC) is a single test and is not
  adjusted.
* Benjamini-Hochberg FDR is applied within each family. A family is all
  the tests of one question (Q2–Q6) across cohorts, measures and time
  points, for one classifier. Both raw p and q are reported.

## 7. Robustness

1. **16S region.**
   * Same samples in two regions: PRJNA773195 (V1-V3) vs PRJNA773202
     (V3-V4), and SRP009299 (V1-V2 vs V4) if it can be processed.
   * Measures: per-sample shift in SMDI, Spearman ρ between regions, and
     the AUC in each region.
   * Meta-regression of the per-cohort AUC on region (V4 vs other), as
     exploratory only.
2. **Classifier.**
   * Every analysis is repeated with vsearch.
   * Agreement between nb and vsearch: Spearman ρ, and the per-sample
     difference (median and IQR).
3. **Rarefaction.**
   * Tables rarefied to 1,000, 5,000 and 10,000 reads, 10 times each
     (fixed seeds). Samples below the depth are dropped.
   * The mean AUC over the 10 draws is reported for the case-control
     cohorts, against the unrarefied AUC. The plugin's documentation says
     not to rarefy. This analysis tests how sensitive the indices are to
     depth, and is not a recommendation.
4. **Pseudocount ε.**
   * The original SMDI (ε = half the smallest non-zero relative abundance
     of the table) is compared with fixed ε = 10⁻⁶, 10⁻⁵ and 10⁻⁴, set
     with the new `--p-smdi-pseudocount` parameter.
   * Measures: per-cohort AUC, Spearman ρ with the original SMDI, and the
     sensitivity and specificity at the published cut-offs.
5. **Leave-one-cohort-out** for the primary meta-analysis.

## 8. Failure handling

* **Processing errors.** A cohort that fails is retried with at most two
  documented changes, chosen from read QC and not from index results (for
  example, a quality-score encoding or error-model setting). If it still
  fails, it is excluded and reported.
* **Budget.** A cohort over about 30 GB or about 12 h is skipped and
  reported.
* **Labels.** If the metadata turn out to be ambiguous (for example,
  groups that cannot be mapped to health and periodontitis), the cohort is
  excluded before any index is computed for it.
* **QC flags.** These are reported, and do not lead to exclusion:
  * median depth after DADA2 below 5,000 reads;
  * median `qc_assigned_fraction` below 0.1;
  * more than 20 of the 49 SMDI species absent from the table.

  A sensitivity analysis drops the flagged cohorts from the primary
  meta-analysis.
* **Too few cohorts.** If fewer than three independent case-control
  cohorts remain, no pooled estimate is reported, only per-cohort
  results.
* **Negative results** are reported with the same detail as positive
  ones.
