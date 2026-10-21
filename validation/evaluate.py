"""Evaluates the q2-periodontitis indices on Chen et al. (2018) plaque data.

    python evaluate.py <export_dir> <metadata.tsv> <results_dir>

Analyses (subjects, not samples, are the unit of analysis):

1. Discrimination of healthy (HP) from periodontitis (D1P) subjects: AUC
   with a 95% stratified bootstrap confidence interval (2000 resamples of
   subjects), Mann-Whitney U p-value, and the sensitivity / specificity of
   SMDI_simplified at its published cut-off of 0 (Chen et al. 2022).
2. Response to non-surgical therapy: per-site change (D2P minus D1P, same
   tooth), Wilcoxon signed-rank test and the fraction of sites whose index
   decreased.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

N_BOOT = 2000
SEED = 20181

export_dir, metadata_path, results_dir = sys.argv[1:4]
os.makedirs(results_dir, exist_ok=True)

meta = pd.read_csv(metadata_path, sep="\t", index_col=0)
indices = pd.read_csv(os.path.join(export_dir, "indices", "indices.tsv"),
                      sep="\t", index_col=0)
complexes = pd.read_csv(os.path.join(export_dir, "complexes", "complexes.tsv"),
                        sep="\t", index_col=0)
mapping = pd.read_csv(os.path.join(export_dir, "mapping", "mapping.tsv"),
                      sep="\t", index_col=0)
deltas = pd.read_csv(os.path.join(export_dir, "deltas", "deltas.tsv"),
                     sep="\t", index_col=0)

complex_cols = ["Red", "Orange", "Yellow", "Green", "Purple", "Blue"]
data = indices.join(complexes[complex_cols].add_suffix("_complex")).join(meta)
metrics = ["SMDI", "SMDI_simplified", "PC2014_moderate"] + [
    c + "_complex" for c in complex_cols]

lines = []


def out(text=""):
    print(text)
    lines.append(text)


# ---------------------------------------------------------------- QC
qc = data[["qc_depth", "qc_assigned_fraction", "qc_species_fraction"]]
out("## Quality control (%d samples)" % len(data))
out()
out("| Measure | Median | IQR |")
out("|---|---|---|")
for col in qc.columns:
    q1, med, q3 = qc[col].quantile([0.25, 0.5, 0.75])
    out("| %s | %.3g | %.3g-%.3g |" % (col, med, q1, q3))
out()
res = mapping["Resolution"].value_counts()
out("Feature resolution: " + ", ".join("%s %d" % kv for kv in res.items()))
out()

# ------------------------------------------------- 1. health vs disease
baseline = data[data["group"].isin(["HP", "D1P"])]
subjects = baseline.groupby(["subject", "status"])[metrics].mean().reset_index()
health = subjects[subjects["status"] == "health"]
disease = subjects[subjects["status"] == "periodontitis"]
subjects.to_csv(os.path.join(results_dir, "subject_indices.tsv"), sep="\t",
                index=False)


def auc(pos, neg):
    """P(disease > health) + 0.5 P(tie)."""
    u, _ = mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (len(pos) * len(neg))


rng = np.random.default_rng(SEED)
rows = []
for metric in metrics:
    pos = disease[metric].dropna().to_numpy()
    neg = health[metric].dropna().to_numpy()
    point = auc(pos, neg)
    boot = [auc(rng.choice(pos, len(pos)), rng.choice(neg, len(neg)))
            for _ in range(N_BOOT)]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    _, p = mannwhitneyu(pos, neg, alternative="two-sided")
    rows.append({"index": metric, "n_health": len(neg),
                 "n_periodontitis": len(pos),
                 "median_health": np.median(neg),
                 "median_periodontitis": np.median(pos),
                 "AUC": point, "AUC_low": lo, "AUC_high": hi, "p": p})
disc = pd.DataFrame(rows)
disc.to_csv(os.path.join(results_dir, "discrimination.tsv"), sep="\t",
            index=False)

out("## Health (HP) vs periodontitis (D1P), subject level")
out()
out("%d healthy and %d periodontitis subjects (mean of their sites)."
    % (len(health), len(disease)))
out()
out("| Index | Median health | Median periodontitis | AUC (95% CI) | p |")
out("|---|---|---|---|---|")
for r in rows:
    out("| %s | %.3g | %.3g | %.2f (%.2f-%.2f) | %.2g |"
        % (r["index"], r["median_health"], r["median_periodontitis"],
           r["AUC"], r["AUC_low"], r["AUC_high"], r["p"]))
out()

# Published cut-offs (Chen et al. 2022): 0 for both indices; 0.4 (training)
# and 0.5 (test) were the best cut-offs of the 49-species SMDI.
for metric, cutoff in [("SMDI", 0.0), ("SMDI", 0.4), ("SMDI", 0.5),
                       ("SMDI_simplified", 0.0)]:
    for level, frame in [("subjects", subjects), ("samples", baseline)]:
        pred = frame[metric] > cutoff
        truth = frame["status"] == "periodontitis"
        sens = (pred & truth).sum() / truth.sum()
        spec = (~pred & ~truth).sum() / (~truth).sum()
        acc = (pred == truth).mean()
        out("%s > %.1f (%s): sensitivity %.2f, specificity %.2f, "
            "accuracy %.2f" % (metric, cutoff, level, sens, spec, acc))
out()

# In-sample optimum (descriptive only: chosen on the same data)
truth = (subjects["status"] == "periodontitis").to_numpy()
for metric in ["SMDI", "SMDI_simplified"]:
    values = subjects[metric].to_numpy()
    best = max(((c, ((values >= c) & truth).sum() / truth.sum(),
                 ((values < c) & ~truth).sum() / (~truth).sum())
                for c in np.unique(values)), key=lambda r: r[1] + r[2])
    out("%s Youden-optimal cut-off (in-sample) %.2f: sensitivity %.2f, "
        "specificity %.2f" % (metric, *best))
for a, b in [("SMDI", "SMDI_simplified"), ("SMDI", "PC2014_moderate")]:
    rho, p = spearmanr(subjects[a], subjects[b])
    out("Spearman %s vs %s: rho %.2f (p %.1g)" % (a, b, rho, p))
out()

# --------------------------------------------------- 2. therapy response
out("## Non-surgical therapy: same site, D2P minus D1P")
out()
out("%d sites (one per subject)." % len(deltas))
out()
out("| Index | Median change | Sites decreased | Wilcoxon p |")
out("|---|---|---|---|")
rows = []
for metric in deltas.columns:
    d = deltas[metric].dropna()
    stat, p = wilcoxon(d)
    rows.append({"index": metric, "n": len(d), "median_change": d.median(),
                 "fraction_decreased": (d < 0).mean(), "wilcoxon_p": p})
    out("| %s | %.3g | %d/%d | %.2g |"
        % (metric, d.median(), (d < 0).sum(), len(d), p))
pd.DataFrame(rows).to_csv(os.path.join(results_dir, "therapy.tsv"), sep="\t",
                          index=False)
out()

with open(os.path.join(results_dir, "summary.md"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
