# `q2-periodontitis`

**QIIME 2 Plugin for Subgingival Biofilm Ecological Characterization, Socransky & Extended Complexes Collapse, and Periodontal Dysbiosis Indices**

[![CI](https://github.com/LaBiOmicS/q2-periodontitis/actions/workflows/ci.yml/badge.svg)](https://github.com/LaBiOmicS/q2-periodontitis/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![QIIME 2](https://img.shields.io/badge/QIIME%202-2024.10%2B-purple.svg)](https://qiime2.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org)
[![Code Style: PEP8](https://img.shields.io/badge/code%20style-PEP%208-orange.svg)](https://peps.python.org/pep-0008/)

---

## 📋 Table of Contents
1. [Overview & Biological Rationale](#1-overview--biological-rationale)
2. [Key Features](#2-key-features)
3. [Ecological Complexes Table (Socransky & Extended)](#3-ecological-complexes-table-socransky--extended)
4. [Mathematical Formulation of Indices](#4-mathematical-formulation-of-indices)
   - [Periodontal Microbial Dysbiosis Index (pMDI & pMDI_ilr)](#41-periodontal-microbial-dysbiosis-index-pmdi--pmdi_ilr)
   - [Socransky Weighted Load (SWL)](#42-socransky-weighted-load-swl)
   - [Keystone Risk Score (KRS)](#43-keystone-risk-score-krs)
   - [Longitudinal Treatment Deltas (Delta-Indices)](#44-longitudinal-treatment-deltas-delta-indices)
5. [Installation](#5-installation)
6. [Usage Guide & CLI Reference](#6-usage-guide--cli-reference)
   - [1. collapse-complexes](#1-collapse-complexes)
   - [2. delta-indices](#2-delta-indices)
   - [3. summarize](#3-summarize)
7. [QIIME 2 Semantic Types & Directory Formats](#7-qiime-2-semantic-types--directory-formats)
8. [Integration with Downstream Pipelines](#8-integration-with-downstream-pipelines)
9. [Citations & License](#9-citations--license)

See also: [docs/THEORY.md](docs/THEORY.md), [docs/TUTORIAL.md](docs/TUTORIAL.md), [docs/SPEC.md](docs/SPEC.md) and the [CHANGELOG](CHANGELOG.md).

---

## 1. Overview & Biological Rationale

Periodontitis is a multifactorial chronic inflammatory disease characterized by progressive destruction of periodontal supporting tissues (periodontal ligament and alveolar bone), driven by dysbiosis in the subgingival biofilm.

In standard amplicon (16S rRNA) or shotgun metagenomic workflows within QIIME 2, taxonomic profiling often fragments interpretation into hundreds or thousands of sparse ASVs/OTUs. While standard alpha and beta diversity metrics provide broad community insights, they fail to capture the **ecological succession** and **functional microbial consortia** well established in oral microbiology.

### Theoretical Foundation
- **Socransky et al. (1998) Model:** Groups subgingival biofilm species into color-coded microbial complexes (*Red*, *Orange*, *Yellow*, *Green*, *Purple*), reflecting early, intermediate, and late colonization stages. The *Blue* complex (*Actinomyces* species) was added in later work by the same group (Socransky & Haffajee, 2002).
- **Keystone-Pathogen Hypothesis (Hajishengallis et al., 2012):** Establishes that low-abundance keystone pathogens (such as *Porphyromonas gingivalis*) orchestrate tissue homeostasis breakdown and elevate destructive inflammatory responses.
- **Extended Complex (plugin-defined):** Groups periodontitis-associated taxa identified after 1998 by culture-independent methods (e.g., *Filifactor alocis*, *Fretibacterium* spp., *Treponema lecithinolyticum*, *Selenomonas sputigena*, *Desulfobulbus* sp.), several of which have moderate evidence of association with periodontitis in the systematic review by Pérez-Chaparro et al. (2014).

`q2-periodontitis` bridges this analytical gap by mapping feature count tables and taxonomic strings into collapsed ecological complex matrices, quantitative dysbiosis metrics, and interactive `.qzv` visual dashboards.

---

## 2. Key Features

- 🧬 **Database-Agnostic Mapping:** Parses taxonomy strings from **eHOMD**, **SILVA**, **Greengenes**, **Greengenes2**, **GTDB**, **NCBI** and consensus taxonomies (from `q2-consensus-taxonomy`). Database name suffixes (SILVA `Prevotella_7`, `Treponema_2`; GTDB/GG2 `Fusobacterium_A`, `nucleatum_D`) are removed before lookup, placeholder names (`uncultured_bacterium`, `metagenome`, `unidentified`, …) are ignored, and historical names are resolved through a synonym table.
- 🎯 **Configurable Resolution Modes (`--p-assignment-mode`):**
  - `strict`: only species-level matches are assigned to a complex.
  - `genus-majority` *(default)*: uses the species when available; otherwise, the genus-level mapping.
- 🔎 **Confidence filter (`--p-min-confidence`):** features whose classifier confidence is below the threshold are reported as `Other` (`low-confidence`), keeping every sample's complex profile summing to 1.
- 📐 **Multidimensional Dysbiosis Metrics:**
  - **pMDI:** log10 ratio of dysbiotic versus eubiotic complexes.
  - **pMDI_ilr:** compositional (isometric log-ratio) balance between the same two groups.
  - **SWL:** load weighted by virulence and proteolytic potential.
  - **KRS:** sigmoidal risk score centred on a keystone abundance threshold.
  - **Delta-Indices:** paired longitudinal treatment monitoring (e.g., pre- vs. post-treatment).
- 📊 **Interactive `.qzv` Visualizations:** Vega-Lite charts (stacked complex bar chart and pMDI × SWL scatter plot, bundled for offline use) and automated non-parametric statistics (Mann-Whitney U / Kruskal-Wallis) for every categorical metadata column, with Benjamini-Hochberg FDR-adjusted *q*-values.
- 🔄 **QIIME 2 Native:** Outputs registered semantic types that can be viewed as QIIME 2 Metadata (e.g., in `q2-longitudinal` and `qiime metadata tabulate`).

---

## 3. Ecological Complexes Table (Socransky & Extended)

The built-in reference is `q2_periodontitis/data/ehomd_map.tsv`. Every entry has a `source` column: `Socransky1998` for the original members, a `citations.bib` key for literature-backed extensions, or `unverified` for extensions that still need a reference.

| Complex | Hex Color | Reference Target Taxa | Clinical Significance & Ecological Role |
| :--- | :--- | :--- | :--- |
| **Red** | `#E41A1C` | *Porphyromonas gingivalis*, *Tannerella forsythia*, *Treponema denticola* | Late-colonizing anaerobes; associated with deep periodontal pockets, bleeding on probing, and active tissue destruction. |
| **Orange** | `#FF7F00` | *Fusobacterium nucleatum* subspecies, *Prevotella intermedia*, *Prevotella nigrescens*, *Parvimonas micra*, *Campylobacter rectus*, *C. showae*, *C. gracilis*, *Eubacterium nodatum*, *Streptococcus constellatus* | Bridge organisms; facilitate co-aggregation and metabolic succession between early colonizers and strict anaerobes. |
| **Extended** | `#984EA3` | *Filifactor alocis*, *Fretibacterium* spp., *Treponema lecithinolyticum*, *Treponema medium*, *Selenomonas sputigena*, *Desulfobulbus* sp., *Aggregatibacter actinomycetemcomitans*†, *Treponema maltophilum*† | Periodontitis-associated taxa described after 1998 (Pérez-Chaparro et al., 2014). |
| **Yellow** | `#FFFF33` | *Streptococcus gordonii*, *S. mitis*, *S. oralis*, *S. sanguinis*, *S. intermedius* | Early colonizers of the acquired pellicle under healthy conditions. |
| **Blue** | `#377EB8` | *Actinomyces naeslundii*, *A. viscosus*, *A. israelii*, *A. gerencseriae* | Early colonizers associated with periodontal health. |
| **Green** | `#4DAF4A` | *Eikenella corrodens*, *Capnocytophaga gingivalis*, *C. ochracea*, *C. sputigena*, *Campylobacter concisus* | Early-to-intermediate colonizers. |
| **Purple** | `#88419D` | *Veillonella parvula*, *Actinomyces odontolyticus* (*Schaalia odontolytica*) | Early colonizers; lactate fermenters. |
| **Other** | `#999999` | Unmapped, unclassified or low-confidence features | Background microflora. |

† Extension without a verified reference yet (`source = unverified`). Genus-level entries (used by `genus-majority`) and other `unverified` extensions are listed in `ehomd_map.tsv`. The genera *Porphyromonas* and *Prevotella* are intentionally **not** mapped at genus level because they contain health-associated species (e.g., *P. catoniae*, *P. melaninogenica*).

You can add or override entries with `--i-custom-mapping` (see [§6](#1-collapse-complexes)).

---

## 4. Mathematical Formulation of Indices

Given a relative abundance matrix **P** of *N* samples and *M* features, and a binary projection matrix **W** of *M* features and *K* complexes (each feature belongs to exactly one complex, with unmapped and low-confidence features in `Other`), the collapsed complex relative abundance matrix **Y** is:

$$\mathbf{Y} = \mathbf{P} \mathbf{W}$$

Each row of **Y** sums to 1 (or 0 for a sample without reads). Let $\mathcal{D} = \{\text{Red}, \text{Orange}, \text{Extended}\}$ (dysbiotic) and $\mathcal{H} = \{\text{Yellow}, \text{Blue}\}$ (eubiotic), and let $\epsilon$ be the pseudocount (`--p-pseudocount`, default $10^{-5}$), a single value applied to every sample on the relative-abundance scale.

### 4.1. Periodontal Microbial Dysbiosis Index (pMDI & pMDI_ilr)

$$pMDI_i = \log_{10} \left( \frac{\sum_{c \in \mathcal{D}} y_{i, c} + \epsilon}{\sum_{c \in \mathcal{H}} y_{i, c} + \epsilon} \right)$$

#### Compositional balance (pMDI_ilr)
The isometric log-ratio balance between the geometric means of the dysbiotic ($r = 3$) and eubiotic ($s = 2$) parts:

$$pMDI_{\text{ilr}, i} = \sqrt{\frac{r s}{r + s}} \, \ln \frac{\left(\prod_{c \in \mathcal{D}} (y_{i,c} + \epsilon)\right)^{1/r}}{\left(\prod_{c \in \mathcal{H}} (y_{i,c} + \epsilon)\right)^{1/s}}, \qquad \sqrt{\tfrac{rs}{r+s}} = \sqrt{\tfrac{6}{5}}$$

(The former `pMDI_clr` was exactly $pMDI_{\text{ilr}} / \sqrt{6/5}$ and has been removed.)

- **Interpretation:**
  - `pMDI > 0`: dysbiotic complexes dominate.
  - `pMDI ~ 0`: transitional community.
  - `pMDI < 0`: eubiotic (commensal) complexes dominate.

Note that zero-abundance complexes contribute $\ln \epsilon$ to `pMDI_ilr`, so its value depends on the pseudocount; compare samples processed with the same `--p-pseudocount`.

---

### 4.2. Socransky Weighted Load (SWL)

$$SWL_i = \sum_{c \in \mathcal{K}} \omega_c \cdot y_{i, c}$$

Empirical Virulence Weights ($\vec{\omega}$):
- $\omega_{\text{Red}} = 1.00$
- $\omega_{\text{Extended}} = 0.90$
- $\omega_{\text{Orange}} = 0.60$
- $\omega_{\text{Green}} = 0.20$
- $\omega_{\text{Purple}} = 0.10$
- $\omega_{\text{Yellow}} = \omega_{\text{Blue}} = \omega_{\text{Other}} = 0.00$

---

### 4.3. Keystone Risk Score (KRS)

Based on the *Keystone-Pathogen Hypothesis*, **KRS** multiplies the pathogenic load by a sigmoid centred on the critical Red complex abundance $\theta$ (`--p-keystone-threshold`, default 0.005 = 0.5%):

$$KRS_i = \left( 1 + \tanh\left( \alpha \cdot \frac{y_{i, \text{Red}} - \theta}{\theta} \right) \right) \cdot \left( y_{i, \text{Red}} + y_{i, \text{Extended}} + 0.5 \cdot y_{i, \text{Orange}} \right)$$

where $\alpha$ (`--p-keystone-steepness`, default 1.0) controls the steepness. The sigmoid factor equals 1 when $y_{\text{Red}} = \theta$, tends to 2 above the threshold and to $1 - \tanh(\alpha)$ when Red is absent, so $\theta$ acts as a real threshold. (Previous versions used $\tanh(\alpha \, y_{\text{Red}}/\theta)$, in which $\theta$ only rescaled the curve.)

---

### 4.4. Longitudinal Treatment Deltas (Delta-Indices)

For each subject $j$ with exactly one baseline and one target sample, and for every index $X \in \{pMDI, pMDI_{\text{ilr}}, SWL, KRS\}$:

$$\Delta X_j = X_{j, \text{Target}} - X_{j, \text{Baseline}}$$

- **Clinical Interpretation:**
  - `Delta < 0`: reduction of dysbiosis (favorable response).
  - `Delta >= 0`: no improvement, recalcitrance, or recurrence.

---

## 5. Installation

### Prerequisites
The plugin is developed and tested against the official **QIIME 2 amplicon 2024.10** distribution:

```bash
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
conda env create -n q2-amplicon-2024.10 --file qiime2-amplicon-2024.10-py310-linux-conda.yml
conda activate q2-amplicon-2024.10
```

### Direct Installation from GitHub

```bash
pip install git+https://github.com/LaBiOmicS/q2-periodontitis.git
```

### Installation for Local Development

```bash
git clone https://github.com/LaBiOmicS/q2-periodontitis.git
cd q2-periodontitis
conda env update -n q2-amplicon-2024.10 --file environment.yml   # adds pytest; QIIME packages are already present
pip install -e .
pytest
```

### Refresh QIIME 2 Plugin Cache

```bash
qiime dev refresh-cache
qiime periodontitis --help
```

---

## 6. Usage Guide & CLI Reference

### 1. `collapse-complexes`
Maps taxonomic strings to Socransky/Extended complexes, collapses the feature count matrix, and computes the dysbiosis indices table.

```bash
qiime periodontitis collapse-complexes \
  --i-table feature-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-assignment-mode genus-majority \
  --p-min-confidence 0.7 \
  --p-pseudocount 1e-5 \
  --p-keystone-threshold 0.005 \
  --p-keystone-steepness 1.0 \
  --o-complex-table socransky-complexes.qza \
  --o-indices dysbiosis-indices.qza \
  --o-mapping feature-socransky-map.qza
```

#### Inputs
- `--i-table` (`FeatureTable[Frequency]`): ASV/OTU feature count table.
- `--i-taxonomy` (`FeatureData[Taxonomy]`): Taxonomic classifications (`Taxon` and, optionally, `Confidence`).
- `--i-custom-mapping` (`SocranskyReference`, *optional*): Taxon-to-complex table that is added to, and overrides, the built-in reference.

#### Parameters
- `--p-assignment-mode` (`Str % Choices('strict', 'genus-majority')`): `strict` assigns species-level matches only; `genus-majority` uses the species when available, otherwise the genus-level mapping. Default: `genus-majority`.
- `--p-min-confidence` (`Float % Range(0, 1)`): Features with `Confidence` below this value go to `Other` with resolution `low-confidence`. Features without a numeric confidence are never filtered. Default: `0.0`.
- `--p-pseudocount` (`Float > 0`): Pseudocount added to every complex relative abundance before the log-ratios. Default: `1e-5`.
- `--p-keystone-threshold` (`Float > 0`): Critical Red complex abundance $\theta$ at which the KRS sigmoid is centred. Default: `0.005`.
- `--p-keystone-steepness` (`Float`): Steepness $\alpha$ of the KRS sigmoid. Default: `1.0`.

#### Outputs
- `--o-complex-table` (`SocranskyComplexes`): Relative abundance of the 8 complexes per sample (rows sum to 1).
- `--o-indices` (`PeriodontalIndices`): pMDI, pMDI_ilr, SWL and KRS per sample.
- `--o-mapping` (`SocranskyMapping`): Audit table with each feature's original taxonomy string, complex, resolution (`species`, `genus`, `unassigned`, `low-confidence`) and confidence.

#### Custom reference (`SocranskyReference`)
A tab-separated file with a mandatory header containing the columns `taxon` and `complex` (extra columns are allowed and ignored). `complex` must be one of `Red`, `Orange`, `Extended`, `Yellow`, `Blue`, `Green`, `Purple`, `Other`; invalid names are rejected on import.

```text
taxon	complex
Selenomonas noxia	Extended
Eubacterium saphenum	Extended
```

```bash
qiime tools import \
  --type SocranskyReference \
  --input-path my-reference.tsv \
  --output-path my-reference.qza
```

---

### 2. `delta-indices`
Calculates paired longitudinal changes ($\Delta$, target minus baseline) of every index for each subject.

```bash
qiime periodontitis delta-indices \
  --i-indices dysbiosis-indices.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-subject-col SubjectID \
  --p-group-col Timepoint \
  --p-baseline-group Pre-Treatment \
  --p-target-group Post-Treatment \
  --o-delta-table delta-indices.qza
```

#### Parameters
- `--i-indices`: Dysbiosis indices artifact generated by `collapse-complexes`.
- `--m-metadata-file`: Sample metadata containing subject identifiers and time points.
- `--p-subject-col`: Column with subject/patient identifiers.
- `--p-group-col`: Column with the treatment group / time point. Categorical and numeric columns are supported (for a numeric column, `--p-baseline-group 0` matches `0` or `0.0`).
- `--p-baseline-group`: Baseline group value (e.g., `Pre-Treatment`).
- `--p-target-group`: Target follow-up group value (e.g., `Post-Treatment`).

#### Output
- `--o-delta-table` (`SampleData[PeriodontalIndicesDelta]`): One row per paired subject. **The index is the subject ID** (column `subject-id`), not a sample ID.

Subjects with only one of the two groups are skipped. An error is raised if any subject has more than one sample in the baseline or target group (the offending subjects are listed), or if no subject can be paired.

---

### 3. `summarize`
Generates an interactive visualization report (`.qzv`).

```bash
qiime periodontitis summarize \
  --i-complexes socransky-complexes.qza \
  --i-indices dysbiosis-indices.qza \
  --i-mapping feature-socransky-map.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization periodontitis-report.qzv
```

When metadata are provided, each index is compared across the groups of every **categorical** metadata column (Mann-Whitney U for two groups, Kruskal-Wallis for more; groups with fewer than 2 samples are ignored). Numeric columns are not tested. *p*-values are adjusted with the Benjamini-Hochberg procedure across all tests in the report; tests with *q* < 0.05 are highlighted. Tests that cannot be computed (e.g., all values identical) are skipped with a warning.

The JavaScript libraries (Vega 5.30.0, Vega-Lite 5.21.0, Vega-Embed 6.26.0) are bundled inside the visualization, so it works offline.

```bash
qiime tools view periodontitis-report.qzv
```

---

## 7. QIIME 2 Semantic Types & Directory Formats

| Semantic type | Description | Directory format (file) |
| :--- | :--- | :--- |
| `SocranskyComplexes` | Relative abundance (*N* × 8) of the ecological complexes per sample. | `SocranskyComplexesDirectoryFormat` (`complexes.tsv`) |
| `PeriodontalIndices` | pMDI, pMDI_ilr, SWL and KRS per sample. | `PeriodontalIndicesDirectoryFormat` (`indices.tsv`) |
| `SampleData[PeriodontalIndicesDelta]` | Per-subject deltas; indexed by subject ID (`subject-id`). | `PeriodontalIndicesDeltaDirectoryFormat` (`deltas.tsv`) |
| `SocranskyMapping` | Feature ID → taxonomy, complex, resolution and confidence. | `SocranskyMappingDirectoryFormat` (`mapping.tsv`) |
| `SocranskyReference` | User taxon → complex reference (`taxon`, `complex` columns). | `SocranskyReferenceDirectoryFormat` (`reference.tsv`) |

`SocranskyComplexes`, `PeriodontalIndices`, `PeriodontalIndicesDelta` and `SocranskyMapping` artifacts can be viewed as QIIME 2 Metadata.

---

## 8. Integration with Downstream Pipelines

Because the indices and complexes can be read as Metadata, they can be passed with `--m-metadata-file` to other plugins.

### Tabulating results
```bash
qiime metadata tabulate \
  --m-input-file dysbiosis-indices.qza \
  --m-input-file sample-metadata.tsv \
  --o-visualization indices-table.qzv
```

### Longitudinal analyses (`q2-longitudinal`)
```bash
qiime longitudinal linear-mixed-effects \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-file dysbiosis-indices.qza \
  --p-metric pMDI \
  --p-state-column Visit \
  --p-individual-id-column PatientID \
  --o-visualization pmdi-lme.qzv
```

The complex table is not a `FeatureTable`, so it cannot be passed directly to `q2-diversity` or `q2-sample-classifier` actions that expect `FeatureTable[Frequency]`.

---

## 9. Citations & License

### Citations
If you use `q2-periodontitis` in your research, please cite (`qiime periodontitis --citations`):

1. Socransky, S. S., Haffajee, A. D., Cugini, M. A., Smith, C., & Kent, R. L. (1998). Microbial complexes in subgingival plaque. *Journal of Clinical Periodontology*, 25(2), 134–144. https://doi.org/10.1111/j.1600-051X.1998.tb02419.x
2. Socransky, S. S., & Haffajee, A. D. (2002). Dental biofilms: difficult therapeutic targets. *Periodontology 2000*, 28(1), 12–55. https://doi.org/10.1034/j.1600-0757.2002.280102.x
3. Hajishengallis, G., Darveau, R. P., & Curtis, M. A. (2012). The keystone-pathogen hypothesis. *Nature Reviews Microbiology*, 10(10), 717–725. https://doi.org/10.1038/nrmicro2873
4. Pérez-Chaparro, P. J., Gonçalves, C., Figueiredo, L. C., Faveri, M., Lobão, E., Tamashiro, N., Duarte, P., & Feres, M. (2014). Newly identified pathogens associated with periodontitis: a systematic review. *Journal of Dental Research*, 93(9), 846–858. https://doi.org/10.1177/0022034514542468
5. Curtis, M. A., Diaz, P. I., & Van Dyke, T. E. (2020). The role of the microbiota in periodontal disease. *Periodontology 2000*, 83(1), 14–25. https://doi.org/10.1111/prd.12296

### License
This project is licensed under the [MIT License](LICENSE). The bundled Vega, Vega-Lite and Vega-Embed libraries are distributed under the BSD 3-Clause License (see `q2_periodontitis/assets/web/js/`).
