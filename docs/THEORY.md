# Theoretical & Mathematical Foundations of `q2-periodontitis`

This document outlines the ecological principles and mathematical formulas implemented in `q2-periodontitis`.

---

## 1. Socransky & Extended Microbial Complexes

The subgingival biofilm community is structured into functional ecological complexes reflecting temporal microbial succession (*Socransky et al., 1998; Curtis et al., 2020*):

- **Red Complex:** *Porphyromonas gingivalis*, *Tannerella forsythia*, *Treponema denticola*. Late obligate anaerobes driving tissue breakdown.
- **Orange Complex:** *Fusobacterium nucleatum*, *Prevotella intermedia*, *Parvimonas micra*, *Campylobacter rectus*, etc. Bridge species enabling co-aggregation.
- **Extended Complex:** *Filifactor alocis*, *Aggregatibacter actinomycetemcomitans*, *Fretibacterium* spp., *Treponema maltophilum*, *Selenomonas sputigena*, etc. Contemporary metagenomic pathogens.
- **Yellow Complex:** *Streptococcus* spp. Pioneer colonizers under eubiotic health.
- **Blue Complex:** *Actinomyces* spp. Early structural matrix formers.
- **Green Complex:** *Capnocytophaga* spp., *Eikenella corrodens*. Transition colonizers.
- **Purple Complex:** *Veillonella* spp. Early lactate fermenters.

---

## 2. Mathematical Formulations

### 2.1. Periodontal Microbial Dysbiosis Index ($pMDI$)
$$pMDI_i = \log_{10} \left( \frac{\sum_{c \in \{\text{Red, Orange, Extended}\}} y_{i, c} + \epsilon_i}{\sum_{c \in \{\text{Yellow, Blue}\}} y_{i, c} + \epsilon_i} \right)$$

### 2.2. Compositional Variant ($pMDI_{\text{clr}}$)
$$pMDI_{\text{clr}, i} = \text{clr}(y_{i, \mathcal{D}}) - \text{clr}(y_{i, \mathcal{H}})$$

### 2.3. Socransky Weighted Load ($SWL$)
$$SWL_i = \sum_{c \in \mathcal{K}} \omega_c \cdot y_{i, c}$$

Virulence weights ($\vec{\omega}$):
Red (1.00), Extended (0.90), Orange (0.60), Green (0.20), Purple (0.10), Yellow (0.00), Blue (0.00), Other (0.00).

### 2.4. Keystone Risk Score ($KRS$)
$$KRS_i = \left( 1 + \tanh\left( \alpha \cdot \frac{y_{i, \text{Red}}}{\theta_{\text{Red}}} \right) \right) \cdot \left( y_{i, \text{Red}} + y_{i, \text{Extended}} + 0.5 \cdot y_{i, \text{Orange}} \right)$$

Default parameters: $\theta_{\text{Red}} = 0.005$ (0.5%), $\alpha = 1.0$.
