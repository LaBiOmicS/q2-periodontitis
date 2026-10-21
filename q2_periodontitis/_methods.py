import os
import numpy as np
import pandas as pd
import biom
from scipy.stats import gmean

COMPLEXES = [
    "Red",
    "Orange",
    "Extended",
    "Yellow",
    "Blue",
    "Green",
    "Purple",
    "Other"
]

WEIGHTS = {
    "Red": 1.00,
    "Extended": 0.90,
    "Orange": 0.60,
    "Green": 0.20,
    "Purple": 0.10,
    "Yellow": 0.00,
    "Blue": 0.00,
    "Other": 0.00
}

ASSIGNMENT_MODES = ["strict", "genus-majority"]

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def load_default_mappings():
    """Loads the canonical tax-complex mapping and synonyms table."""
    ehomd_path = os.path.join(DATA_DIR, "ehomd_map.tsv")
    syn_path = os.path.join(DATA_DIR, "synonyms_map.tsv")

    tax_map = {}
    if os.path.exists(ehomd_path):
        df_ehomd = pd.read_csv(ehomd_path, sep="\t")
        for _, row in df_ehomd.iterrows():
            taxon = str(row["taxon"]).strip()
            complex_name = str(row["complex"]).strip()
            tax_map[taxon] = complex_name

    synonyms = {}
    if os.path.exists(syn_path):
        df_syn = pd.read_csv(syn_path, sep="\t")
        for _, row in df_syn.iterrows():
            old_name = str(row["old_name"]).strip()
            canon_name = str(row["canonical_name"]).strip()
            synonyms[old_name] = canon_name

    return tax_map, synonyms


def parse_taxonomy_string(tax_str, synonyms=None):
    """
    Parses QIIME 2 / Greengenes / SILVA taxonomy string into (genus, species).
    Applies synonym cross-walk if available.
    """
    if synonyms is None:
        synonyms = {}

    if not isinstance(tax_str, str) or not tax_str.strip():
        return None, None

    # Split by semicolon
    ranks = [r.strip() for r in tax_str.split(";")]
    genus = None
    species = None

    for r in ranks:
        clean_r = r
        # Remove rank prefixes like d__, k__, p__, c__, o__, f__, g__, s__
        for prefix in ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]:
            if clean_r.startswith(prefix):
                clean_r = clean_r[len(prefix):].strip()
                break

        if not clean_r or clean_r.lower() in ["unassigned", "uncultured", "ambiguous_taxa"]:
            continue

        if r.startswith("g__") or (genus is None and len(ranks) >= 6 and r == ranks[5]):
            genus = clean_r
        elif r.startswith("s__") or (species is None and len(ranks) >= 7 and r == ranks[6]):
            species = clean_r

    # If species contains genus (e.g. s__Porphyromonas_gingivalis or s__Porphyromonas gingivalis), format binomial
    if species:
        species_clean = species.replace("_", " ")
        if genus and not species_clean.startswith(genus):
            species_full = f"{genus} {species_clean}"
        else:
            species_full = species_clean
    else:
        species_full = None

    if genus:
        genus = genus.replace("_", " ")

    # Apply synonym lookup
    if species_full in synonyms:
        species_full = synonyms[species_full]
    if genus in synonyms:
        genus = synonyms[genus]

    return genus, species_full


def build_lookup(tax_map):
    """
    Validates a taxon -> complex map and returns a case-insensitive lookup
    (lower-cased taxon -> complex). Raises ValueError on unknown complexes.
    """
    invalid = sorted({c for c in tax_map.values() if c not in COMPLEXES})
    if invalid:
        offenders = sorted(t for t, c in tax_map.items() if c in invalid)
        raise ValueError(
            "Invalid complex name(s) %s (taxa: %s). Valid complexes are: %s."
            % (", ".join(repr(c) for c in invalid),
               ", ".join(offenders[:10]) + (" ..." if len(offenders) > 10 else ""),
               ", ".join(COMPLEXES))
        )
    return {k.lower(): v for k, v in tax_map.items()}


def assign_feature_to_complex(tax_str, lookup, synonyms, assignment_mode="genus-majority"):
    """
    Assigns a single taxonomic string to a Socransky / Extended complex.
    `lookup` is the case-insensitive map returned by build_lookup().

    assignment_mode:
      - "strict": species-level match only;
      - "genus-majority": uses the species when available; otherwise, the
        genus-level mapping.
    """
    if assignment_mode not in ASSIGNMENT_MODES:
        raise ValueError("assignment_mode must be one of %s, got %r."
                         % (", ".join(ASSIGNMENT_MODES), assignment_mode))
    genus, species = parse_taxonomy_string(tax_str, synonyms)

    # 1. Species level match
    if species and species.lower() in lookup:
        return lookup[species.lower()], "species"

    # 2. Genus level match (genus-majority only)
    if assignment_mode == "genus-majority":
        if genus and genus.lower() in lookup:
            return lookup[genus.lower()], "genus"

    return "Other", "unassigned"


def collapse_complexes(
    table: biom.Table,
    taxonomy: pd.DataFrame,
    assignment_mode: str = "genus-majority",
    custom_mapping: pd.DataFrame = None,
    min_confidence: float = 0.0,
    keystone_threshold: float = 0.005,
    keystone_steepness: float = 1.0
) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Collapses ASV/OTU feature table into Socransky complexes and computes dysbiosis metrics.
    """
    tax_map, synonyms = load_default_mappings()

    if custom_mapping is not None and not custom_mapping.empty:
        missing = [c for c in ("taxon", "complex") if c not in custom_mapping.columns]
        if missing:
            raise ValueError(
                "custom_mapping must have the columns 'taxon' and 'complex'; "
                "missing: %s." % ", ".join(missing))
        for t_name, c_name in zip(custom_mapping["taxon"], custom_mapping["complex"]):
            tax_map[str(t_name).strip()] = str(c_name).strip()

    lookup = build_lookup(tax_map)

    # Taxonomy arrives as the q2-types DataFrame (Taxon + optional Confidence)
    if isinstance(taxonomy, pd.Series):
        taxonomy = taxonomy.to_frame("Taxon")
    taxonomy_series = taxonomy["Taxon"] if "Taxon" in taxonomy.columns \
        else taxonomy.iloc[:, 0]
    confidence_series = None
    if "Confidence" in taxonomy.columns:
        # Missing / non-numeric confidences carry no information and are
        # never filtered out.
        confidence_series = pd.to_numeric(taxonomy["Confidence"],
                                          errors="coerce")

    sample_ids = list(table.ids(axis="sample"))
    feature_ids = list(table.ids(axis="observation"))

    # Convert biom table to matrix: features x samples
    matrix = table.matrix_data.toarray()  # (M x N)

    # Relative abundance matrix P: samples x features
    col_sums = matrix.sum(axis=0)
    col_sums[col_sums == 0] = 1.0  # Prevent division by zero
    rel_abundance = (matrix / col_sums).T  # (N x M)

    # Map each feature to exactly one complex. Features whose classifier
    # confidence is below min_confidence are sent to "Other" so that each
    # sample's complex profile still sums to 1.
    feature_assignments = []
    w_matrix = np.zeros((len(feature_ids), len(COMPLEXES)))

    for idx, f_id in enumerate(feature_ids):
        tax_str = str(taxonomy_series.get(f_id, ""))
        confidence = np.nan
        if confidence_series is not None:
            confidence = float(confidence_series.get(f_id, np.nan))

        if not np.isnan(confidence) and confidence < min_confidence:
            complex_name, resolution = "Other", "low-confidence"
        else:
            complex_name, resolution = assign_feature_to_complex(
                tax_str, lookup, synonyms, assignment_mode
            )
        c_idx = COMPLEXES.index(complex_name)
        w_matrix[idx, c_idx] = 1.0
        feature_assignments.append({
            "Feature ID": f_id,
            "Taxonomy": tax_str,
            "Complex": complex_name,
            "Resolution": resolution,
            "Confidence": confidence
        })

    # Collapsed complexes matrix Y = P x W (N x K)
    collapsed_matrix = rel_abundance @ w_matrix
    df_complexes = pd.DataFrame(
        collapsed_matrix,
        index=sample_ids,
        columns=COMPLEXES
    )

    # Compute Dysbiosis Indices (pMDI, pMDI_clr, pMDI_ilr, SWL, KRS)
    dysbiotic = df_complexes["Red"] + df_complexes["Orange"] + df_complexes["Extended"]
    health = df_complexes["Yellow"] + df_complexes["Blue"]

    pMDI_vals = []
    pMDI_clr_vals = []
    pMDI_ilr_vals = []
    SWL_vals = []
    KRS_vals = []

    # ILR scaling factor: sqrt( (3 * 2) / (3 + 2) ) = sqrt(6 / 5)
    ilr_factor = np.sqrt(6.0 / 5.0)

    for i in range(len(sample_ids)):
        row = df_complexes.iloc[i]
        d_val = dysbiotic.iloc[i]
        h_val = health.iloc[i]

        # Dynamic pseudocount epsilon_i
        non_zero = row[row > 0]
        epsilon_i = (non_zero.min() / 2.0) if len(non_zero) > 0 else 1e-5

        # Standard pMDI
        pmdi = np.log10((d_val + epsilon_i) / (h_val + epsilon_i))
        pMDI_vals.append(pmdi)

        # pMDI_clr (Centered Log-Ratio)
        row_eps = row + epsilon_i
        g_mean = np.exp(np.mean(np.log(row_eps)))
        clr_row = np.log(row_eps / g_mean)
        clr_d = np.mean([clr_row["Red"], clr_row["Orange"], clr_row["Extended"]])
        clr_h = np.mean([clr_row["Yellow"], clr_row["Blue"]])
        pMDI_clr_vals.append(clr_d - clr_h)

        # pMDI_ilr (Isometric Log-Ratio Balance)
        d_gmean = np.exp(np.mean(np.log([row["Red"] + epsilon_i, row["Orange"] + epsilon_i, row["Extended"] + epsilon_i])))
        h_gmean = np.exp(np.mean(np.log([row["Yellow"] + epsilon_i, row["Blue"] + epsilon_i])))
        ilr_val = ilr_factor * np.log(d_gmean / h_gmean)
        pMDI_ilr_vals.append(ilr_val)

        # SWL (Socransky Weighted Load)
        swl = sum(WEIGHTS[c] * row[c] for c in COMPLEXES)
        SWL_vals.append(swl)

        # KRS (Keystone Risk Score)
        red_val = row["Red"]
        krs = (1.0 + np.tanh(keystone_steepness * (red_val / keystone_threshold))) * (red_val + row["Extended"] + 0.5 * row["Orange"])
        KRS_vals.append(krs)

    df_indices = pd.DataFrame({
        "pMDI": pMDI_vals,
        "pMDI_clr": pMDI_clr_vals,
        "pMDI_ilr": pMDI_ilr_vals,
        "SWL": SWL_vals,
        "KRS": KRS_vals
    }, index=sample_ids)

    df_mapping = pd.DataFrame(feature_assignments).set_index("Feature ID")

    return df_complexes, df_indices, df_mapping


def delta_indices(
    indices: pd.DataFrame,
    metadata: object,
    sample_id_col: str,
    group_col: str,
    baseline_group: str,
    target_group: str
) -> pd.DataFrame:
    """
    Computes longitudinal delta (Pre vs Post treatment) for pMDI, SWL, KRS.
    """
    if hasattr(metadata, "to_dataframe"):
        metadata_df = metadata.to_dataframe()
    else:
        metadata_df = metadata

    merged = indices.join(metadata_df, how="inner")
    if sample_id_col not in merged.columns or group_col not in merged.columns:
        raise ValueError(f"Columns '{sample_id_col}' and '{group_col}' must be present in metadata.")

    deltas = []
    grouped = merged.groupby(sample_id_col)

    for subject_id, group_df in grouped:
        baseline = group_df[group_df[group_col] == baseline_group]
        target = group_df[group_df[group_col] == target_group]

        if not baseline.empty and not target.empty:
            b_row = baseline.iloc[0]
            t_row = target.iloc[0]
            deltas.append({
                "Subject_ID": subject_id,
                "pMDI": t_row["pMDI"] - b_row["pMDI"],
                "pMDI_clr": t_row["pMDI_clr"] - b_row["pMDI_clr"],
                "SWL": t_row["SWL"] - b_row["SWL"],
                "KRS": t_row["KRS"] - b_row["KRS"]
            })

    if not deltas:
        df_res = pd.DataFrame(columns=["pMDI", "pMDI_clr", "SWL", "KRS"])
        df_res.index.name = "sample-id"
        return df_res

    df_res = pd.DataFrame(deltas).set_index("Subject_ID")
    df_res.index.name = "sample-id"
    return df_res

