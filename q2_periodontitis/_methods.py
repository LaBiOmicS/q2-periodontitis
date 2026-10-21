import os
import re
import numpy as np
import pandas as pd
import biom
import qiime2

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

DYSBIOTIC = ["Red", "Orange", "Extended"]
EUBIOTIC = ["Yellow", "Blue"]

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


RANK_PREFIX_RE = re.compile(r"^[a-z]__")

# Database suffixes appended to genus names or epithets: SILVA
# (Prevotella_7, Treponema_2) and GTDB / Greengenes2 (Fusobacterium_A,
# Prevotella_D, nucleatum_D).
NAME_SUFFIX_RE = re.compile(r"_(\d+|[A-Z]{1,2})$")
SUFFIX_TOKEN_RE = re.compile(r"\d+|[A-Z]{1,2}")
# Tokens after which a number is an identifier, not a suffix ("oral taxon 361")
ID_PREFIX_TOKENS = {"taxon", "hmt", "hot", "ot", "str.", "strain"}

UNINFORMATIVE_PREFIXES = (
    "uncultured", "unidentified", "unclassified", "unknown", "unassigned",
    "ambiguous", "metagenome", "environmental", "incertae"
)


def _is_uninformative(name):
    """True for placeholder labels such as 'uncultured_bacterium',
    'metagenome', 'unidentified' or 'Ambiguous_taxa'."""
    lower = name.strip().lower()
    return (not lower
            or lower.startswith(UNINFORMATIVE_PREFIXES)
            or "metagenome" in lower
            or lower in {"bacterium", "sp", "sp."})


def _clean_name(raw):
    """
    Normalizes a genus or species label: splits on spaces/underscores and
    drops database suffix tokens (Prevotella_7 -> Prevotella,
    Fusobacterium_A nucleatum_D -> Fusobacterium nucleatum).
    """
    tokens = [t for t in re.split(r"[\s_]+", raw.strip()) if t]
    kept = []
    for tok in tokens:
        if (kept and SUFFIX_TOKEN_RE.fullmatch(tok)
                and kept[-1].lower() not in ID_PREFIX_TOKENS):
            continue
        kept.append(tok)
    return " ".join(kept)


def parse_taxonomy_string(tax_str, synonyms=None):
    """
    Parses QIIME 2 / Greengenes / Greengenes2 / GTDB / SILVA taxonomy strings
    into (genus, species). Database suffixes (SILVA "_7", GTDB "_A") are
    removed and placeholder names ("uncultured", "metagenome", ...) ignored.
    Applies synonym cross-walk if available.
    """
    if synonyms is None:
        synonyms = {}

    if not isinstance(tax_str, str) or not tax_str.strip():
        return None, None

    ranks = [r.strip() for r in tax_str.split(";")]
    has_prefixes = any(RANK_PREFIX_RE.match(r) for r in ranks)
    genus_raw = None
    species_raw = None

    for pos, r in enumerate(ranks):
        clean_r = RANK_PREFIX_RE.sub("", r).strip()
        if _is_uninformative(clean_r):
            continue
        if has_prefixes:
            if r.startswith("g__"):
                genus_raw = clean_r
            elif r.startswith("s__"):
                species_raw = clean_r
        elif pos == 5:
            genus_raw = clean_r
        elif pos == 6:
            species_raw = clean_r

    genus = _clean_name(genus_raw).split(" ")[0] if genus_raw else None

    species_full = None
    if species_raw:
        species_clean = _clean_name(species_raw)
        if species_clean[:1].isupper():
            # Binomial already present (s__Porphyromonas_gingivalis)
            species_full = species_clean
        elif genus:
            # Epithet only (Greengenes 13_8: g__Porphyromonas; s__gingivalis)
            species_full = f"{genus} {species_clean}"
        if species_full and len(species_full.split(" ")) < 2:
            species_full = None

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
    pseudocount: float = 1e-5,
    keystone_threshold: float = 0.005,
    keystone_steepness: float = 1.0
) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Collapses ASV/OTU feature table into Socransky complexes and computes dysbiosis metrics.
    """
    if not 0.0 <= min_confidence <= 1.0:
        raise ValueError("min_confidence must be within [0, 1].")
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

    # Collapsed complexes matrix Y = P W (N x K), computed sparsely as
    # (X^T W) / depth so the feature table is never densified.
    counts = table.matrix_data  # sparse (M x N)
    depth = np.asarray(counts.sum(axis=0)).ravel()
    depth[depth == 0] = 1.0  # Empty samples keep an all-zero profile
    collapsed_matrix = np.asarray(counts.T @ w_matrix) / depth[:, None]
    df_complexes = pd.DataFrame(
        collapsed_matrix,
        index=sample_ids,
        columns=COMPLEXES
    )

    df_indices = compute_indices(
        df_complexes, pseudocount, keystone_threshold, keystone_steepness
    )

    df_mapping = pd.DataFrame(feature_assignments).set_index("Feature ID")

    return df_complexes, df_indices, df_mapping


def compute_indices(complexes, pseudocount=1e-5, keystone_threshold=0.005,
                    keystone_steepness=1.0):
    """
    Computes pMDI, pMDI_ilr, SWL and KRS from a samples x complexes table of
    relative abundances. `pseudocount` is added to every complex abundance of
    every sample (relative-abundance scale) before taking logs.
    """
    if pseudocount <= 0:
        raise ValueError("pseudocount must be > 0.")
    if keystone_threshold <= 0:
        raise ValueError("keystone_threshold must be > 0.")

    y = complexes
    dysbiotic = y[DYSBIOTIC].sum(axis=1)
    health = y[EUBIOTIC].sum(axis=1)

    # pMDI: log10 ratio of dysbiotic vs eubiotic complex totals
    pmdi = np.log10((dysbiotic + pseudocount) / (health + pseudocount))

    # pMDI_ilr: isometric log-ratio balance between the geometric means of
    # the dysbiotic (r=3) and eubiotic (s=2) parts; sqrt(r*s/(r+s)) = sqrt(6/5)
    log_y = np.log(y[DYSBIOTIC + EUBIOTIC] + pseudocount)
    ilr_factor = np.sqrt(len(DYSBIOTIC) * len(EUBIOTIC)
                         / (len(DYSBIOTIC) + len(EUBIOTIC)))
    pmdi_ilr = ilr_factor * (log_y[DYSBIOTIC].mean(axis=1)
                             - log_y[EUBIOTIC].mean(axis=1))

    # SWL: virulence-weighted load
    swl = y[COMPLEXES] @ pd.Series(WEIGHTS)[COMPLEXES]

    # KRS: sigmoid centred on the keystone threshold theta, so the penalty
    # factor is 1 at Red == theta, -> 2 above it and -> 1 - tanh(alpha) at 0.
    red = y["Red"]
    sigmoid = 1.0 + np.tanh(keystone_steepness
                            * (red - keystone_threshold) / keystone_threshold)
    krs = sigmoid * (red + y["Extended"] + 0.5 * y["Orange"])

    return pd.DataFrame({
        "pMDI": pmdi,
        "pMDI_ilr": pmdi_ilr,
        "SWL": swl,
        "KRS": krs
    }, index=complexes.index)


def _normalize_label(value):
    """String form of a metadata value; integral floats lose their '.0'
    so numeric time points (e.g. 0.0) match labels such as '0'."""
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return str(int(value))
    return str(value)


def delta_indices(
    indices: pd.DataFrame,
    metadata: qiime2.Metadata,
    subject_col: str,
    group_col: str,
    baseline_group: str,
    target_group: str
) -> pd.DataFrame:
    """
    Computes paired longitudinal deltas (target minus baseline) of every
    index in `indices`, one row per subject. The output index is the subject
    ID taken from `subject_col`.
    """
    metadata_df = metadata.to_dataframe()
    missing = [c for c in (subject_col, group_col) if c not in metadata_df.columns]
    if missing:
        raise ValueError("Metadata column(s) not found: %s." % ", ".join(missing))

    # Join only the two needed columns, under private names, to avoid any
    # collision with index or metadata column names.
    design = pd.DataFrame({
        "__subject__": metadata_df[subject_col],
        "__group__": metadata_df[group_col],
    }).dropna()
    design["__subject__"] = design["__subject__"].map(_normalize_label)
    design["__group__"] = design["__group__"].map(_normalize_label)

    baseline_label = str(baseline_group)
    target_label = str(target_group)
    if metadata.columns[group_col].type == "numeric":
        # Let "1" match a numeric column value of 1.0
        try:
            baseline_label = _normalize_label(float(baseline_group))
            target_label = _normalize_label(float(target_group))
        except ValueError:
            raise ValueError(
                "Column '%s' is numeric, so baseline_group and target_group "
                "must be numbers (got %r and %r)."
                % (group_col, baseline_group, target_group)) from None

    metrics = list(indices.columns)
    merged = indices.join(design, how="inner")
    merged = merged[merged["__group__"].isin([baseline_label, target_label])]

    counts = merged.groupby(["__subject__", "__group__"]).size()
    duplicated = sorted(set(counts[counts > 1].index.get_level_values(0)))
    if duplicated:
        raise ValueError(
            "More than one sample per subject and group for subject(s): %s. "
            "Each subject must have at most one '%s' and one '%s' sample."
            % (", ".join(duplicated), baseline_group, target_group))

    baseline = merged[merged["__group__"] == baseline_label].set_index("__subject__")[metrics]
    target = merged[merged["__group__"] == target_label].set_index("__subject__")[metrics]
    paired = baseline.index.intersection(target.index)
    if len(paired) == 0:
        raise ValueError(
            "No subject in '%s' has samples in both groups '%s' and '%s' of "
            "column '%s'." % (subject_col, baseline_group, target_group, group_col))

    df_res = target.loc[paired] - baseline.loc[paired]
    df_res = df_res.sort_index()
    df_res.index.name = "subject-id"
    return df_res
