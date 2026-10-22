"""Shared definitions of the validation analysis."""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
COHORTS_TSV = os.path.join(HERE, "..", "cohorts.tsv")

CLASSIFIERS = ["nb", "vsearch"]
PRIMARY_CLASSIFIER = "nb"
COMPLEXES = ["Red", "Orange", "Yellow", "Green", "Purple", "Blue"]
PRIMARY = "SMDI"
COMPARATOR = "complex_Red"
MEASURES = (["SMDI", "SMDI_simplified", "PC2014_moderate"]
            + ["complex_" + c for c in COMPLEXES])
# Published cut-offs (Chen et al. 2022); positive if index > cut-off
CUTOFFS = [("SMDI", 0.0), ("SMDI", 0.4), ("SMDI", 0.5),
           ("SMDI_simplified", 0.0)]
CLINICAL = ["probing_depth", "clinical_attachment_loss",
            "bleeding_on_probing", "plaque_index", "stage"]
N_BOOT = 2000
SEED = 20261022


def load_cohorts():
    return pd.read_csv(COHORTS_TSV, sep="\t", dtype=str).fillna("")


def load_samples(results, cohort, clf):
    path = os.path.join(results, "samples", "%s-%s.tsv" % (cohort, clf))
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", dtype={"sample-id": str, "subject": str,
                                            "unit": str, "time": str,
                                            "status": str, "arm": str,
                                            "response": str})
    for col in CLINICAL:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def subject_case_control(df, measures=MEASURES):
    """Baseline samples of health and periodontitis subjects, averaged per
    subject (analysis plan, section 4)."""
    d = df[(df["time"] == "baseline")
           & df["status"].isin(["health", "periodontitis"])]
    return d.groupby(["subject", "status"])[measures].mean().reset_index()


def absent_smdi_taxa(mapping_path):
    """Number of the 49 SMDI species without any feature in the table."""
    from q2_periodontitis._methods import (SMDI_SPECIES, load_signature,
                                           load_synonyms, taxon_key,
                                           _canonical)
    synonyms = load_synonyms()
    mapping = pd.read_csv(mapping_path, sep="\t", dtype=str).fillna("")
    present = {taxon_key(s) for s in mapping["Species"] if s and "/" not in s}
    signature = load_signature(SMDI_SPECIES)
    keys = [taxon_key(_canonical(t, synonyms)) for t in signature["taxon"]]
    return sum(k not in present for k in keys)
