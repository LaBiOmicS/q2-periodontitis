# `q2-periodontitis`

**QIIME 2 Plugin for Subgingival Biofilm Ecological Characterization, Socransky & Extended Complexes Collapse, and Periodontal Dysbiosis Indices**

---

## 📋 Table of Contents
1. [Overview & Biological Rationale](#1-overview--biological-rationale)
2. [Key Features](#2-key-features)
3. [Ecological Complexes Table (Socransky & Extended)](#3-ecological-complexes-table-socransky--extended)
4. [Mathematical Formulation of Indices](#4-mathematical-formulation-of-indices)
   - [Periodontal Microbial Dysbiosis Index ($pMDI$ & $pMDI_{\text{clr}}$)](#41-periodontal-microbial-dysbiosis-index-pmdi--pmdi_textclr)
   - [Socransky Weighted Load ($SWL$)](#42-socransky-weighted-load-swl)
   - [Keystone Risk Score ($KRS$)](#43-keystone-risk-score-krs)
   - [Longitudinal Treatment Deltas ($\Delta$-Indices)](#44-longitudinal-treatment-deltas-delta-indices)
5. [Installation](#5-installation)
6. [Usage Guide & CLI Reference](#6-usage-guide--cli-reference)
   - [1. `collapse-complexes`](#1-collapse-complexes)
   - [2. `delta-indices`](#2-delta-indices)
   - [3. `summarize`](#3-summarize)
7. [QIIME 2 Semantic Types & Directory Formats](#7-qiime-2-semantic-types--directory-formats)
8. [Integration with Downstream Pipelines](#8-integration-with-downstream-pipelines)
9. [Citations & License](#9-citations--license)

---

## 1. Overview & Biological Rationale

Periodontitis is a multifactorial chronic inflammatory disease characterized by progressive destruction of periodontal supporting tissues (periodontal ligament and alveolar bone), driven by dysbiosis in the subgingival biofilm.

In standard amplicon (16S rRNA) or shotgun metagenomic workflows within QIIME 2, taxonomic profiling often fragments interpretation into hundreds or thousands of sparse ASVs/OTUs. While standard alpha and beta diversity metrics provide broad community insights, they fail to capture the **ecological succession** and **functional microbial consortia** well established in oral microbiology.

### Theoretical Foundation
- **Socransky et al. (1998) Model:** Groups subgingival biofilm species into color-coded microbial consortia (*Red*, *Orange*, *Yellow*, *Green*, *Purple*, *Blue*), reflecting early, intermediate, and late colonization stages.
- **Keystone-Pathogen Hypothesis (Hajishengallis et al., 2012):** Establishes that low-abundance keystone pathogens (such as *Porphyromonas gingivalis*) orchestrate tissue homeostasis breakdown and elevate destructive inflammatory responses.
- **Extended Complex (Curtis et al., 2020):** Incorporates contemporary metagenomic pathogens (e.g., *Filifactor alocis*, *Fretibacterium* spp., *Treponema maltophilum*, *Selenomonas sputigena*, *Aggregatibacter actinomycetemcomitans*).

`q2-periodontitis` bridges this analytical gap by mapping feature count tables and taxonomic strings into collapsed ecological complex matrices, quantitative dysbiosis metrics, and interactive `.qzv` visual dashboards.

---

## 2. Key Features

- 🧬 **Database-Agnostic Mapping:** Native support for **eHOMD**, **SILVA**, **Greengenes2**, **GTDB**, **NCBI**, and consensus taxonomies (from `q2-consensus-taxonomy`).
- 🎯 **Configurable Resolution Modes (`--p-assignment-mode`):**
  - `strict`: Requires exact species-level (`s__`) match.
  - `genus-majority` *(Recommended)*: Resolves at the genus level (`g__`) for conserved oral taxa.
  - `probabilistic`: Proportionally weights feature abundance based on reference database species distributions.
- 📐 **Multidimensional Dysbiosis Metrics:**
  - **$pMDI$:** Log-ratio of dysbiotic versus eubiotic complexes (with adaptive pseudocount and Compositional Data Analysis CLR support).
  - **$SWL$:** Load weighted by virulency and proteolytic potential.
  - **$KRS$:** Sigmoidal risk score penalizing the presence of low-abundance keystone pathogens.
  - **$\Delta$-Indices:** Longitudinal treatment monitoring (Pre- vs. Post-treatment / follow-up).
- 📊 **Interactive `.qzv` Visualizations:** Dashboards equipped with Vega-Lite charts (stacked bar charts, index boxplots) and automated non-parametric statistics (Kruskal-Wallis and Mann-Whitney U tests with $p$-values).
- 🔄 **QIIME 2 Ecosystem Native:** Directly outputs registered semantic types, enabling seamless chaining into `q2-diversity`, `q2-sample-classifier`, and `q2-longitudinal`.

---

## 3. Ecological Complexes Table (Socransky & Extended)

| Complex | Hex Color | Reference Target Taxa | Clinical Significance & Ecological Role |
| :--- | :--- | :--- | :--- |
| **Red** | `#E41A1C` | *Porphyromonas gingivalis*, *Tannerella forsythia*, *Treponema denticola* | Mandatory late-stage obligate anaerobes; associated with deep periodontal pockets, bleeding on probing, and active tissue destruction. |
| **Orange** | `#FF7F00` | *Fusobacterium nucleatum*, *Prevotella intermedia*, *Prevotella nigrescens*, *Parvimonas micra*, *Campylobacter rectus*, *Eubacterium nodatum*, *Streptococcus constellatus* | Bridge organisms; facilitate physical co-aggregation and metabolic succession between primary colonizers and strict anaerobes. |
| **Extended** | `#984EA3` | *Filifactor alocis*, *Aggregatibacter actinomycetemcomitans*, *Fretibacterium fastidiosum*, *Treponema maltophilum*, *Selenomonas sputigena*, *Desulfobulbus* spp. | Contemporary metagenomic pathogens; linked with rapidly progressing or refractory periodontitis. |
| **Yellow** | `#FFFF33` | *Streptococcus gordonii*, *Streptococcus mitis*, *Streptococcus oralis*, *Streptococcus sanguinis*, *Streptococcus parasanguinis* | Early pioneer colonizers; bind salivary mucins of the acquired pellicle under healthy sulcular conditions. |
| **Blue** | `#377EB8` | *Actinomyces naeslundii*, *Actinomyces viscosus*, *Actinomyces israelii*, *Actinomyces gerencseriae*, *Actinomyces odontolyticus* | Initial structural matrix formers maintaining eubiotic periodontal homeostasis. |
| **Green** | `#4DAF4A` | *Eikenella corrodens*, *Capnocytophaga gingivalis*, *Capnocytophaga ochracea*, *Capnocytophaga sputigena*, *Campylobacter concisus* | Early-to-intermediate transition colonizers; range from commensals to mild opportunistic pathogens. |
| **Purple** | `#88419D` | *Veillonella parvula*, *Veillonella atypica*, *Veillonella dispar*, *Veillonella rogosae* | Lactate fermenters; facilitate metabolic succession without causing direct tissue damage. |
| **Other** | `#999999` | Taxa not associated with periodontal complexes or unclassified | Background microflora. |

---

## 4. Mathematical Formulation of Indices

Given a relative abundance matrix $\mathbf{P} \in [0, 1]^{N \times M}$ across $N$ samples and $M$ features, and a binary projection matrix $\mathbf{W} \in \{0, 1\}^{M \times K}$, the collapsed complex matrix $\mathbf{Y} \in [0, 1]^{N \times K}$ is calculated as:
$$\mathbf{Y} = \mathbf{P} \mathbf{W}$$

### 4.1. Periodontal Microbial Dysbiosis Index ($pMDI$ & $pMDI_{\text{clr}}$)

The $pMDI$ computes the log-ratio between dysbiotic/pathogenic complexes ($\mathcal{D} = \{\text{Red}, \text{Orange}, \text{Extended}\}$) and eubiotic/commensal complexes ($\mathcal{H} = \{\text{Yellow}, \text{Blue}\}$):

$$pMDI_i = \log_{10} \left( \frac{\sum_{c \in \mathcal{D}} y_{i, c} + \epsilon_i}{\sum_{c \in \mathcal{H}} y_{i, c} + \epsilon_i} \right)$$

Where $\epsilon_i$ is a sample-wise dynamic adaptive pseudocount:
$$\epsilon_i = \begin{cases} \frac{\min(y_{i, c} > 0)}{2}, & \text{if } \exists \, y_{i, c} > 0 \\ 10^{-5}, & \text{otherwise} \end{cases}$$

#### Compositional Variant ($pMDI_{\text{clr}}$)
For robust analysis under Compositional Data Analysis (CoDA) principles:
$$pMDI_{\text{clr}, i} = \text{clr}(y_{i, \mathcal{D}}) - \text{clr}(y_{i, \mathcal{H}})$$

- **Interpretation:**
  - $pMDI > 0$: Severe dysbiotic state (high pathogen burden).
  - $pMDI \approx 0$: Transitional biofilm / mild inflammatory state.
  - $pMDI < 0$: Periodontal homeostasis (commensal dominance).

---

### 4.2. Socransky Weighted Load ($SWL$)

The $SWL$ calculates pathogenic burden weighted by proteolytic and virulence potential:

$$SWL_i = \sum_{c \in \mathcal{K}} \omega_c \cdot y_{i, c}, \quad SWL_i \in [0, 1]$$

Empirical Virulence Weights ($\vec{\omega}$):
- $\omega_{\text{Red}} = 1.00$
- $\omega_{\text{Extended}} = 0.90$
- $\omega_{\text{Orange}} = 0.60$
- $\omega_{\text{Green}} = 0.20$
- $\omega_{\text{Purple}} = 0.10$
- $\omega_{\text{Yellow}} = 0.00$
- $\omega_{\text{Blue}} = 0.00$
- $\omega_{\text{Other}} = 0.00$

---

### 4.3. Keystone Risk Score ($KRS$)

Based on the *Keystone-Pathogen Hypothesis*, $KRS$ applies a sigmoidal penalty when Red complex pathogen abundance surpasses the critical threshold $\theta_{\text{Red}} = 0.005$ (0.5%):

$$KRS_i = \left( 1 + \tanh\left( \alpha \cdot \frac{y_{i, \text{Red}}}{\theta_{\text{Red}}} \right) \right) \cdot \left( y_{i, \text{Red}} + y_{i, \text{Extended}} + 0.5 \cdot y_{i, \text{Orange}} \right)$$

Where $\alpha$ controls sigmoidal steepness (default: $\alpha = 1.0$).

---

### 4.4. Longitudinal Treatment Deltas ($\Delta$-Indices)

Allows comparison of paired samples from subjects undergoing periodontal therapy (e.g., Scaling and Root Planing - SRP, adjunctive antibiotics, or phototherapy):

$$\Delta pMDI_i = pMDI_{i, \text{Target}} - pMDI_{i, \text{Baseline}}$$
$$\Delta SWL_i = SWL_{i, \text{Target}} - SWL_{i, \text{Baseline}}$$
$$\Delta KRS_i = KRS_{i, \text{Target}} - KRS_{i, \text{Baseline}}$$

- **Clinical Interpretation:**
  - $\Delta < 0$: Favorable therapeutic response (reduction in dysbiosis and eubiotic biofilm recovery).
  - $\Delta \ge 0$: Treatment failure, recalcitrance, or disease recurrence.

---

## 5. Installation

### Prerequisites
- An active **QIIME 2** environment (e.g., QIIME 2 2024.10 or compatible).

### Direct Installation from GitHub

In your active QIIME 2 environment, install directly from source:

```bash
pip install git+https://github.com/LaBiOmicS/q2-periodontitis.git
```

### Installation for Local Development

Alternatively, clone the repository and install in editable mode:

```bash
git clone https://github.com/LaBiOmicS/q2-periodontitis.git
cd q2-periodontitis
pip install -e .
```

### Refresh QIIME 2 Plugin Cache

After installation, update the QIIME 2 command-line interface cache:

```bash
qiime dev refresh-cache
```

Verify installation:
```bash
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
  --p-keystone-threshold 0.005 \
  --p-keystone-steepness 1.0 \
  --o-complex-table socransky-complexes.qza \
  --o-indices dysbiosis-indices.qza \
  --o-mapping feature-socransky-map.qza
```

#### Inputs
- `--i-table` (`FeatureTable[Frequency]`): ASV/OTU feature count table.
- `--i-taxonomy` (`FeatureData[Taxonomy]`): Taxonomic classifications for features.
- `--i-custom-mapping` (`FeatureData[SocranskyMapping]`, *Optional*): Custom taxon-to-complex mapping table.

#### Parameters
- `--p-assignment-mode` (`Choice['strict', 'genus-majority', 'probabilistic']`): Taxonomic assignment resolution strategy. Default: `genus-majority`.
- `--p-keystone-threshold` (`Float`): Critical abundance threshold ($\theta_{\text{Red}}$) for Red complex pathogens. Default: `0.005`.
- `--p-keystone-steepness` (`Float`): Steepness factor ($\alpha$) for the sigmoidal KRS calculation. Default: `1.0`.

#### Outputs
- `--o-complex-table` (`SampleData[SocranskyComplexes]`): Relative abundance matrix of collapsed ecological complexes.
- `--o-indices` (`SampleData[PeriodontalIndices]`): Calculated dysbiosis metrics ($pMDI$, $pMDI_{\text{clr}}$, $SWL$, $KRS$).
- `--o-mapping` (`FeatureData[SocranskyMapping]`): Audit table mapping each feature ID to its assigned complex.

---

### 2. `delta-indices`
Calculates longitudinal treatment changes ($\Delta$) for dysbiosis metrics between baseline and follow-up samples.

```bash
qiime periodontitis delta-indices \
  --i-indices dysbiosis-indices.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-sample-id-col SubjectID \
  --p-group-col Timepoint \
  --p-baseline-group Pre-Treatment \
  --p-target-group Post-Treatment \
  --o-delta-table delta-indices.qza
```

#### Parameters
- `--i-indices`: Dysbiosis indices artifact generated by `collapse-complexes`.
- `--m-metadata-file`: Sample metadata containing subject identifiers and time points.
- `--p-sample-id-col`: Column name for subject/patient identifiers.
- `--p-group-col`: Column name for treatment group / time point identifiers.
- `--p-baseline-group`: Name of the baseline group (e.g., `Pre-Treatment`).
- `--p-target-group`: Name of the target follow-up group (e.g., `Post-Treatment`).

---

### 3. `summarize`
Generates an interactive visualization report artifact (`.qzv`).

```bash
qiime periodontitis summarize \
  --i-complexes socransky-complexes.qza \
  --i-indices dysbiosis-indices.qza \
  --i-mapping feature-socransky-map.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization periodontitis-report.qzv
```

To view the interactive report in a browser:
```bash
qiime tools view periodontitis-report.qzv
```

---

## 7. QIIME 2 Semantic Types & Directory Formats

`q2-periodontitis` registers three dedicated **Semantic Types** within the QIIME 2 framework:

1. **`SampleData[SocranskyComplexes]`**
   - **Description:** Relative abundance matrix ($N \times K$) containing proportions across all 8 ecological complexes per sample.
   - **Directory Format:** `SocranskyComplexesDirectoryFormat` (`matrix.tsv`).

2. **`SampleData[PeriodontalIndices]`**
   - **Description:** Quantitative table containing calculated metrics ($pMDI$, $pMDI_{\text{clr}}$, $SWL$, $KRS$) for each sample.
   - **Directory Format:** `PeriodontalIndicesDirectoryFormat` (`indices.tsv`).

3. **`FeatureData[SocranskyMapping]`**
   - **Description:** Auditable feature mapping table tracking each ASV/OTU, its assigned complex, and resolution level (`species`, `genus`, `unassigned`).
   - **Directory Format:** `SocranskyMappingDirectoryFormat` (`mapping.tsv`).

---

## 8. Integration with Downstream Pipelines

Artifacts produced by `q2-periodontitis` connect seamlessly into standard QIIME 2 plugins:

### Ordination and $\beta$-Diversity (`q2-diversity`)
Perform PCoA or UMAP ordination on collapsed complex matrices to visualize community trajectories:
```bash
qiime diversity pcoa \
  --i-distance-matrix complexes-braycurtis.qza \
  --o-pcoa complexes-pcoa.qza
```

### Machine Learning & Clinical Classification (`q2-sample-classifier`)
Use complex matrices and dysbiosis indices to train predictive machine learning models (Random Forest) for periodontitis staging:
```bash
qiime sample-classifier classify-samples \
  --i-table socransky-complexes.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column PeriodontitisStage \
  --o-execution-result rf-results.qza
```

### Longitudinal Trajectory Modeling (`q2-longitudinal`)
Evaluate temporal dynamics of dysbiosis metrics across clinical follow-up cohorts:
```bash
qiime longitudinal linear-mixed-effects \
  --i-table dysbiosis-indices.qza \
  --m-metadata-file sample-metadata.tsv \
  --p-metric pMDI \
  --p-state-column Visit \
  --p-individual-id-column PatientID \
  --o-visualization linear-effects.qzv
```

---

## 9. Citations & License

### Citations
If you use `q2-periodontitis` in your research, please cite:

1. **Socransky et al. (1998):**
   > Socransky, S. S., Haffajee, A. D., Cugini, M. A., Smith, C., & Kent, R. L. (1998). Microbial complexes in subgingival plaque. *Journal of Clinical Periodontology*, 25(2), 134-144.
2. **Hajishengallis et al. (2012):**
   > Hajishengallis, G., Darveau, R. P., & Curtis, M. A. (2012). The keystone-pathogen hypothesis. *Nature Reviews Microbiology*, 10(10), 717-725.
3. **Curtis et al. (2020):**
   > Curtis, M. A., Diaz, P. I., & Van Dyke, T. E. (2020). Biomarkers of periodontal disease status and progression. *Periodontology 2000*, 84(1), 11-25.

### License
This project is licensed under the [BSD-3-Clause License](file:///home/fabiano.menegidio/workdir/TI/bioinformatics/q2-plugins/q2-periodontitis/LICENSE).
