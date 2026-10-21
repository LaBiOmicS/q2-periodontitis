## Quality control (157 samples)

| Measure | Median | IQR |
|---|---|---|
| qc_depth | 2.86e+04 | 2.4e+04-3.39e+04 |
| qc_assigned_fraction | 0.286 | 0.196-0.397 |
| qc_species_fraction | 0.827 | 0.745-0.888 |

Feature resolution: unassigned 2645, species 119, genus 28

## Health (HP) vs periodontitis (D1P), subject level

21 healthy and 48 periodontitis subjects (mean of their sites).

| Index | Median health | Median periodontitis | AUC (95% CI) | p |
|---|---|---|---|---|
| SMDI | -0.0652 | 2.43 | 0.90 (0.81-0.97) | 1.6e-07 |
| SMDI_simplified | -3.68 | 1.26 | 0.86 (0.76-0.95) | 1.8e-06 |
| PC2014_moderate | 0.0147 | 0.0953 | 0.84 (0.71-0.94) | 1.1e-05 |
| Red_complex | 0.00225 | 0.0988 | 0.86 (0.75-0.95) | 2.3e-06 |
| Orange_complex | 0.121 | 0.162 | 0.71 (0.56-0.84) | 0.0065 |
| Yellow_complex | 0.00323 | 0 | 0.18 (0.08-0.31) | 1.7e-05 |
| Green_complex | 0.00788 | 0.00466 | 0.38 (0.25-0.52) | 0.13 |
| Purple_complex | 0.0122 | 0.00217 | 0.21 (0.09-0.35) | 0.00015 |
| Blue_complex | 0.0465 | 0.00871 | 0.16 (0.07-0.25) | 6.6e-06 |

SMDI > 0.0 (subjects): sensitivity 0.90, specificity 0.52, accuracy 0.78
SMDI > 0.0 (samples): sensitivity 0.82, specificity 0.62, accuracy 0.76
SMDI > 0.4 (subjects): sensitivity 0.88, specificity 0.76, accuracy 0.84
SMDI > 0.4 (samples): sensitivity 0.79, specificity 0.71, accuracy 0.77
SMDI > 0.5 (subjects): sensitivity 0.85, specificity 0.90, accuracy 0.87
SMDI > 0.5 (samples): sensitivity 0.79, specificity 0.74, accuracy 0.78
SMDI_simplified > 0.0 (subjects): sensitivity 0.60, specificity 0.90, accuracy 0.70
SMDI_simplified > 0.0 (samples): sensitivity 0.62, specificity 0.88, accuracy 0.70

SMDI Youden-optimal cut-off (in-sample) 0.54: sensitivity 0.85, specificity 0.90
SMDI_simplified Youden-optimal cut-off (in-sample) -1.46: sensitivity 0.85, specificity 0.81
Spearman SMDI vs SMDI_simplified: rho 0.91 (p 5e-27)
Spearman SMDI vs PC2014_moderate: rho 0.73 (p 1e-12)

## Non-surgical therapy: same site, D2P minus D1P

19 sites (one per subject).

| Index | Median change | Sites decreased | Wilcoxon p |
|---|---|---|---|
| SMDI | -0.251 | 13/19 | 0.096 |
| SMDI_simplified | -0.566 | 12/19 | 0.29 |
| PC2014_moderate | -0.0188 | 13/19 | 0.13 |

