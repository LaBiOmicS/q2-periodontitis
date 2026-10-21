# Theoretical & Mathematical Foundations of `q2-periodontitis`

This document outlines the ecological principles and mathematical formulas implemented in `q2-periodontitis`.

---

## 1. Socransky & Extended Microbial Complexes

The subgingival biofilm community is structured into ecological complexes reflecting temporal microbial succession (*Socransky et al., 1998*; *Socransky & Haffajee, 2002*):

- **Red Complex:** *Porphyromonas gingivalis*, *Tannerella forsythia*, *Treponema denticola*. Late colonizers associated with tissue breakdown.
- **Orange Complex:** *Fusobacterium nucleatum* subspecies, *Prevotella intermedia*, *P. nigrescens*, *Parvimonas micra*, *Campylobacter rectus*, *C. showae*, *C. gracilis*, *Eubacterium nodatum*, *Streptococcus constellatus*. Bridge species enabling co-aggregation.
- **Yellow Complex:** *Streptococcus gordonii*, *S. mitis*, *S. oralis*, *S. sanguinis*, *S. intermedius*. Early colonizers.
- **Green Complex:** *Capnocytophaga* spp., *Eikenella corrodens*, *Campylobacter concisus*. Early colonizers.
- **Purple Complex:** *Veillonella parvula* and *Actinomyces odontolyticus* (now *Schaalia odontolytica*). Early colonizers.
- **Blue Complex:** *Actinomyces* species (*A. naeslundii*, *A. viscosus*, *A. israelii*, *A. gerencseriae*), described by Socransky & Haffajee (2002).
- **Extended Complex (plugin-defined):** periodontitis-associated taxa identified after 1998 by culture-independent methods, e.g. *Filifactor alocis*, *Fretibacterium* spp., *Treponema lecithinolyticum*, *T. medium*, *Selenomonas sputigena* and *Desulfobulbus* sp. (moderate evidence in Pérez-Chaparro et al., 2014).

Each entry of the built-in reference (`q2_periodontitis/data/ehomd_map.tsv`) records its `source`: `Socransky1998`, a `citations.bib` key, or `unverified` for extensions that still need a supporting reference.

### 1.1. Taxonomy parsing

For each feature, the genus (`g__`) and species (`s__`) ranks are extracted (or ranks 6 and 7 when there are no rank prefixes) and normalized:

1. Placeholder labels (`uncultured*`, `unidentified`, `unclassified`, `metagenome`, `Ambiguous_taxa`, …) are ignored.
2. Database suffixes on genus names and epithets are removed, matching `_(\d+|[A-Z]{1,2})$` — SILVA (`Prevotella_7`, `Treponema_2`) and GTDB / Greengenes2 (`Fusobacterium_A`, `Prevotella_D`, `nucleatum_D`). Identifiers such as `oral taxon 361` are kept.
3. An epithet-only species (Greengenes 13_8, `s__gingivalis`) is combined with the genus.
4. Historical names are replaced through `synonyms_map.tsv` (e.g., *Bacteroides forsythus* → *Tannerella forsythia*).

The original taxonomy string is preserved in the mapping output.

### 1.2. Assignment modes

- `strict`: the species must be in the reference.
- `genus-majority` (default): uses the species when available; otherwise, the genus-level mapping.

Features with a classifier `Confidence` below `min_confidence` are assigned to `Other` (resolution `low-confidence`). Each feature belongs to exactly one complex.

---

## 2. Mathematical Formulations

Let $y_{i,c}$ be the relative abundance of complex $c$ in sample $i$ ($\sum_c y_{i,c} = 1$, or 0 for an empty sample), $\mathcal{D} = \{\text{Red}, \text{Orange}, \text{Extended}\}$, $\mathcal{H} = \{\text{Yellow}, \text{Blue}\}$, and $\epsilon$ the pseudocount (`pseudocount`, default $10^{-5}$), identical for all samples.

### 2.1. Periodontal Microbial Dysbiosis Index ($pMDI$)
$$pMDI_i = \log_{10} \left( \frac{\sum_{c \in \mathcal{D}} y_{i, c} + \epsilon}{\sum_{c \in \mathcal{H}} y_{i, c} + \epsilon} \right)$$

### 2.2. Compositional balance ($pMDI_{\text{ilr}}$)
$$pMDI_{\text{ilr}, i} = \sqrt{\frac{6}{5}} \left( \frac{1}{3} \sum_{c \in \mathcal{D}} \ln (y_{i,c} + \epsilon) - \frac{1}{2} \sum_{c \in \mathcal{H}} \ln (y_{i,c} + \epsilon) \right)$$

This is the isometric log-ratio balance between the dysbiotic ($r = 3$) and eubiotic ($s = 2$) parts, with normalizing constant $\sqrt{rs/(r+s)}$. The difference of mean CLR coordinates of the same groups (the former `pMDI_clr`) is identical up to the factor $\sqrt{6/5}$ and is therefore no longer reported.

### 2.3. Socransky Weighted Load ($SWL$)
$$SWL_i = \sum_{c \in \mathcal{K}} \omega_c \cdot y_{i, c}$$

Virulence weights ($\vec{\omega}$):
Red (1.00), Extended (0.90), Orange (0.60), Green (0.20), Purple (0.10), Yellow (0.00), Blue (0.00), Other (0.00).

### 2.4. Keystone Risk Score ($KRS$)
$$KRS_i = \left( 1 + \tanh\left( \alpha \cdot \frac{y_{i, \text{Red}} - \theta}{\theta} \right) \right) \cdot \left( y_{i, \text{Red}} + y_{i, \text{Extended}} + 0.5 \cdot y_{i, \text{Orange}} \right)$$

Default parameters: $\theta = 0.005$ (0.5%), $\alpha = 1.0$. The sigmoid factor is 1 at $y_{\text{Red}} = \theta$, approaches 2 above the threshold and equals $1 - \tanh(\alpha)$ in the absence of Red complex taxa.

> **Change from earlier versions:** the factor was $1 + \tanh(\alpha\, y_{\text{Red}}/\theta)$, which is already ≥ 1.76 at $y_{\text{Red}} = \theta$ (for $\alpha = 1$), so $\theta$ only rescaled the curve and did not behave as a threshold.

### 2.5. Longitudinal deltas
For each subject with exactly one baseline and one target sample: $\Delta X = X_{\text{target}} - X_{\text{baseline}}$ for $X \in \{pMDI, pMDI_{\text{ilr}}, SWL, KRS\}$.

### 2.6. Group comparisons in `summarize`
For each categorical metadata column and index: Mann-Whitney U (two-sided) for two groups or Kruskal-Wallis for more (groups with ≥ 2 samples). All *p*-values in a report are adjusted with the Benjamini-Hochberg procedure; *q* < 0.05 is highlighted.

---

## 3. References

- Socransky SS, Haffajee AD, Cugini MA, Smith C, Kent RL Jr. Microbial complexes in subgingival plaque. *J Clin Periodontol* 1998;25(2):134–144. doi:10.1111/j.1600-051X.1998.tb02419.x
- Socransky SS, Haffajee AD. Dental biofilms: difficult therapeutic targets. *Periodontol 2000* 2002;28:12–55. doi:10.1034/j.1600-0757.2002.280102.x
- Hajishengallis G, Darveau RP, Curtis MA. The keystone-pathogen hypothesis. *Nat Rev Microbiol* 2012;10(10):717–725. doi:10.1038/nrmicro2873
- Pérez-Chaparro PJ, et al. Newly identified pathogens associated with periodontitis: a systematic review. *J Dent Res* 2014;93(9):846–858. doi:10.1177/0022034514542468
