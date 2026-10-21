# Technical specification

The scientific definitions are in [THEORY.md](THEORY.md). This document
describes how the code is organised.

## Modules

| File | Content |
| :--- | :--- |
| `q2_periodontitis/_methods.py` | Taxonomy parsing, reference handling, complex assignment, CLR / SMDI, `collapse_complexes`, `delta_indices`. |
| `q2_periodontitis/_visualizer.py` | `summarize`: report rendering, analysis units, group tests, BH correction. |
| `q2_periodontitis/_type.py`, `_format.py`, `_transformer.py` | Semantic types, file formats and their validation, DataFrame / Metadata transformers. |
| `q2_periodontitis/plugin_setup.py` | Action and type registration. |
| `q2_periodontitis/data/socransky_complexes.tsv` | Built-in reference: `taxon`, `complex`, `level`, `source`, `evidence`. |
| `q2_periodontitis/data/synonyms_map.tsv` | Historical → current names. |
| `q2_periodontitis/data/perezchaparro2014_moderate.tsv` | Moderate-evidence taxa with HMT numbers. |
| `q2_periodontitis/data/smdi_chen2022_species.tsv` | Built-in SMDI signature (49 species, MDG ≥ 0.4), with original labels and MDG. |
| `q2_periodontitis/data/smdi_simplified_genus.tsv` | Built-in simplified SMDI signature (3 genera). |
| `q2_periodontitis/assets/web/` | Jinja2 template, CSS and bundled Vega libraries. |
| `validation/` | Reproducible validation on SRA SRP075100 (not installed with the package). |

## `collapse_complexes` algorithm

1. Merge the built-in reference with `custom_mapping` (custom rows win,
   case-insensitively). Build a species lookup and a genus lookup; only
   rows with `level = genus` go to the genus lookup.
2. For each feature:
   1. Parse genus and species from the taxonomy string
      (`parse_taxonomy_string`). This strips rank prefixes, database
      suffixes and bracketed genus markers, skips placeholders, expands
      epithet-only species, and applies synonyms.
   2. If `Confidence < min_confidence`, assign `Other` / `low-confidence`.
   3. Otherwise look up the species (all candidates, for `a/b` labels),
      then the genus (`assign_feature_to_complex`).
   4. Flag Pérez-Chaparro taxa by name or HMT number.
3. Collapse to complexes with a sparse product `countsᵀ · W`, where `W` is
   the binary feature × complex matrix. Divide by depth to get proportions.
4. Compute SMDI, SMDI_simplified and, given a signature, SMDI_custom with
   `compute_smdi`. This aggregates counts by name at the signature level
   (species are keyed by `taxon_key`: HMT number or lower-case name) with
   a sparse indicator product, applies `clr_microbiome` and takes the
   difference of the group means.
5. Add `PC2014_moderate` and the `qc_*` columns.

Empty samples give all-zero complex rows and NaN indices.

## Conventions

* Index columns whose name starts with `qc_` are quality control. They are
  excluded from `delta_indices` and from the tests in `summarize`
  (`index_columns`).
* Numeric table formats accept empty fields as missing values.
* Every built-in reference row must cite a key of `citations.bib`. This is
  enforced by `TestReference.test_reference_is_documented`.
