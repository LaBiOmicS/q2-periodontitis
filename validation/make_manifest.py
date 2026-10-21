"""Builds the QIIME 2 manifest and sample metadata for SRP075100 plaque.

Sample titles are GROUP_SUBJECT_SITE: HP (healthy), D1P (periodontitis,
before therapy) and D2P (after therapy); SITE is the tooth number plus a
suffix ("R" marks the re-sampled site after therapy).
"""
import os
import re
import sys

import pandas as pd

runs_path, fastq_dir, manifest_path, metadata_path = sys.argv[1:5]
cols = ["run_accession", "sample_accession", "sample_alias", "sample_title",
        "library_layout", "instrument_model", "read_count", "base_count",
        "fastq_ftp", "fastq_bytes"]
runs = pd.read_csv(runs_path, sep="\t", header=None, names=cols)

rows = []
for run, title in zip(runs["run_accession"], runs["sample_title"]):
    group, subject, site = title.split("_")
    tooth = re.match(r"\d+", site).group(0)
    rows.append({
        "sample-id": title,
        "run": run,
        "group": group,
        "status": "health" if group == "HP" else "periodontitis",
        "visit": "after-therapy" if group == "D2P" else "baseline",
        "subject": f"{group[0]}_{subject}",
        "site": f"{group[0]}_{subject}_{tooth}",
        "absolute-filepath": os.path.abspath(
            os.path.join(fastq_dir, f"{run}.fastq.gz")),
    })
df = pd.DataFrame(rows)
df[["sample-id", "absolute-filepath"]].to_csv(manifest_path, sep="\t",
                                              index=False)
df.drop(columns=["absolute-filepath"]).to_csv(metadata_path, sep="\t",
                                              index=False)
print(df.groupby("group").agg(samples=("sample-id", "size"),
                              subjects=("subject", "nunique")))
