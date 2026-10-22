"""Collects the q2-periodontitis outputs of every processed cohort.

    python -I collect.py <data_dir> <results_dir>

For each cohort of cohorts.tsv and each classifier, joins the indices, the
complex fractions and the standardised metadata into
<results_dir>/samples/<cohort>-<classifier>.tsv, and writes the quality
control summary <results_dir>/qc.tsv with the QC flags of the analysis
plan (section 8). Cohorts without outputs are listed as not processed.
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import (CLASSIFIERS, COMPLEXES, load_cohorts,  # noqa: E402
                    absent_smdi_taxa)

data_dir, results = sys.argv[1:3]
os.makedirs(os.path.join(results, "samples"), exist_ok=True)
cohorts = load_cohorts()
qc_rows = []

for _, c in cohorts.iterrows():
    base = os.path.join(data_dir, c["cohort"])
    meta_path = os.path.join(base, "metadata.tsv")
    for clf in CLASSIFIERS:
        exp = os.path.join(base, "out", "export-" + clf)
        ind_path = os.path.join(exp, "indices", "indices.tsv")
        if not (os.path.exists(ind_path) and os.path.exists(meta_path)):
            qc_rows.append({"cohort": c["cohort"], "classifier": clf,
                            "processed": "no"})
            continue
        meta = pd.read_csv(meta_path, sep="\t", dtype=str).set_index(
            "sample-id")
        ind = pd.read_csv(ind_path, sep="\t", index_col=0)
        com = pd.read_csv(os.path.join(exp, "complexes", "complexes.tsv"),
                          sep="\t", index_col=0)
        ind.index = ind.index.astype(str)
        com.index = com.index.astype(str)
        data = ind.join(com[COMPLEXES].add_prefix("complex_"))
        data = data.join(meta, how="inner")
        dropped = len(meta) - len(data)
        data.index.name = "sample-id"
        data.to_csv(os.path.join(results, "samples",
                                 "%s-%s.tsv" % (c["cohort"], clf)), sep="\t")
        stats_path = os.path.join(base, "out", "export-denoising-stats",
                                  "stats.tsv")
        retained = float("nan")
        if os.path.exists(stats_path):
            st = pd.read_csv(stats_path, sep="\t", comment="#", index_col=0)
            st = st[st.index != "q2:types"]
            if "reads-raw" in st.columns:   # Deblur sample stats
                first, last = st["reads-raw"], st["reads-hit-reference"]
            else:                           # DADA2
                first, last = st["input"], st["non-chimeric"]
            retained = (pd.to_numeric(last, errors="coerce")
                        / pd.to_numeric(first, errors="coerce")).median()
        absent = absent_smdi_taxa(os.path.join(exp, "mapping",
                                               "mapping.tsv"))
        depth = data["qc_depth"].median()
        assigned = data["qc_assigned_fraction"].median()
        flags = []
        if depth < 5000:
            flags.append("depth<5000")
        if assigned < 0.1:
            flags.append("assigned<0.1")
        if absent > 20:
            flags.append("SMDI_absent>20")
        qc_rows.append({
            "cohort": c["cohort"], "classifier": clf, "processed": "yes",
            "samples": len(data), "samples_without_output": dropped,
            "subjects": data["subject"].nunique(),
            "median_depth": depth,
            "median_reads_retained": retained,
            "median_assigned_fraction": assigned,
            "median_species_fraction": data["qc_species_fraction"].median(),
            "smdi_species_absent": absent,
            "qc_flags": ";".join(flags)})

pd.DataFrame(qc_rows).to_csv(os.path.join(results, "qc.tsv"), sep="\t",
                             index=False)
print(pd.DataFrame(qc_rows).to_string())
