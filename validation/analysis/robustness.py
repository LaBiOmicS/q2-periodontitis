"""Depth and pseudocount robustness (docs/ANALYSIS_PLAN.md, section 7.3-7.4).

    python -I robustness.py <data_dir> <results_dir>

For every processed case-control cohort and classifier, recomputes the
indices with q2-periodontitis from the feature table and taxonomy:

* pseudocount: the original SMDI against fixed epsilon = 1e-6, 1e-5, 1e-4
  (collapse_complexes(smdi_pseudocount=...));
* rarefaction: 10 draws (fixed seeds) at 1,000, 5,000 and 10,000 reads,
  samples below the depth dropped, original pseudocount.

Subject-level AUCs (health vs periodontitis) are written to
<results_dir>/tables/robustness_epsilon.tsv and
<results_dir>/tables/robustness_rarefaction.tsv.
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
import qiime2
from scipy import stats as sps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stats as st  # noqa: E402
from common import (CLASSIFIERS, COMPLEXES, CUTOFFS, SEED,  # noqa: E402
                    load_cohorts, subject_case_control)
from q2_periodontitis._methods import collapse_complexes  # noqa: E402

EPSILONS = [1e-6, 1e-5, 1e-4]
DEPTHS = [1000, 5000, 10000]
DRAWS = 10
MEASURES = ["SMDI", "SMDI_simplified", "PC2014_moderate", "complex_Red"]

data_dir, results = sys.argv[1:3]
tables = os.path.join(results, "tables")
os.makedirs(tables, exist_ok=True)
warnings.simplefilter("ignore")


def indices(table, taxonomy, **kwargs):
    complexes, ind, _ = collapse_complexes(table, taxonomy, **kwargs)
    return ind.join(complexes[COMPLEXES].add_prefix("complex_"))


def case_control(ind, meta):
    data = ind.join(meta, how="inner")
    return subject_case_control(data.reset_index(), MEASURES)


def aucs(subj):
    out = {}
    for m in MEASURES:
        pos = subj.loc[subj["status"] == "periodontitis", m].dropna()
        neg = subj.loc[subj["status"] == "health", m].dropna()
        out[m] = (st.auc(pos, neg), len(neg), len(pos)) if (
            len(pos) >= 2 and len(neg) >= 2) else (np.nan, len(neg), len(pos))
    return out


eps_rows, rar_rows = [], []
cohorts = load_cohorts()
for c in cohorts[cohorts["case_control"] == "yes"].itertuples():
    out = os.path.join(data_dir, c.cohort, "out")
    meta_path = os.path.join(data_dir, c.cohort, "metadata.tsv")
    if not os.path.exists(os.path.join(out, "table.qza")):
        continue
    meta = pd.read_csv(meta_path, sep="\t", dtype=str).set_index("sample-id")
    import biom  # noqa: F401  (view type)
    table = qiime2.Artifact.load(os.path.join(out, "table.qza")).view(
        biom.Table)
    for clf in CLASSIFIERS:
        tax_path = os.path.join(out, "taxonomy-%s.qza" % clf)
        if not os.path.exists(tax_path):
            continue
        taxonomy = qiime2.Artifact.load(tax_path).view(pd.DataFrame)
        print(c.cohort, clf, flush=True)
        original = indices(table, taxonomy)
        subj0 = case_control(original, meta)
        for label, ind in [("original", original)] + [
                ("%g" % e, indices(table, taxonomy, smdi_pseudocount=e))
                for e in EPSILONS]:
            subj = case_control(ind, meta)
            for m in ["SMDI", "SMDI_simplified"]:
                pos = subj.loc[subj["status"] == "periodontitis", m].dropna()
                neg = subj.loc[subj["status"] == "health", m].dropna()
                both = pd.concat([original[m], ind[m]], axis=1).dropna()
                row = {"cohort": c.cohort, "classifier": clf,
                       "meta": c.meta, "epsilon": label, "measure": m,
                       "n_health": len(neg), "n_periodontitis": len(pos),
                       "auc": st.auc(pos, neg),
                       "spearman_vs_original": sps.spearmanr(
                           both.iloc[:, 0], both.iloc[:, 1])[0],
                       "median_shift_vs_original": (
                           both.iloc[:, 1] - both.iloc[:, 0]).median()}
                for mm, cut in CUTOFFS:
                    if mm != m:
                        continue
                    truth = subj["status"] == "periodontitis"
                    pred = subj[m] > cut
                    row["sens_%g" % cut] = (pred & truth).sum() / truth.sum()
                    row["spec_%g" % cut] = (~pred & ~truth).sum() / (
                        ~truth).sum()
                eps_rows.append(row)
        base = aucs(subj0)
        for depth in DEPTHS:
            keep = [s for s, n in zip(table.ids(), table.sum(axis="sample"))
                    if n >= depth]
            sub = table.filter(keep, inplace=False)
            for draw in range(DRAWS):
                rar = sub.subsample(depth, axis="sample", by_id=False,
                                    seed=SEED + draw)
                res = aucs(case_control(indices(rar, taxonomy), meta))
                for m, (a, nh, npos) in res.items():
                    rar_rows.append({"cohort": c.cohort, "classifier": clf,
                                     "meta": c.meta, "depth": depth,
                                     "draw": draw, "measure": m,
                                     "n_health": nh, "n_periodontitis": npos,
                                     "auc": a, "auc_unrarefied": base[m][0]})
    pd.DataFrame(eps_rows).to_csv(os.path.join(tables,
                                               "robustness_epsilon.tsv"),
                                  sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(rar_rows).to_csv(os.path.join(tables,
                                               "robustness_rarefaction.tsv"),
                                  sep="\t", index=False, float_format="%.6g")
print("done")
