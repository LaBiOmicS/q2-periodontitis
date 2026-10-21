# Changelog

All notable changes to `q2-periodontitis` are documented in this file.

## [Unreleased] — next version

### ⚠️ Breaking changes
- **`pMDI_clr` removed.** It was exactly `pMDI_ilr / sqrt(6/5)`. Use `pMDI_ilr`, the canonical compositional balance.
- **Pseudocount.** The per-sample dynamic pseudocount (half of the smallest non-zero complex abundance) was replaced by the `--p-pseudocount` parameter (default `1e-5`), applied identically to all samples. `pMDI` and `pMDI_ilr` values change.
- **KRS formula.** Now `(1 + tanh(alpha * (Red - theta) / theta)) * (Red + Extended + 0.5 * Orange)`, so `theta` (`--p-keystone-threshold`) is a real threshold (factor = 1 at `Red == theta`). KRS values change. `--p-keystone-threshold` must be > 0.
- **Confidence handling.** Abundances are no longer multiplied by the taxonomy `Confidence` (which lost mass, so complex profiles did not sum to 1). Use the new `--p-min-confidence` to send low-confidence features to `Other`.
- **`--p-assignment-mode probabilistic` removed** (it was never implemented and behaved like `genus-majority`).
- **`--i-custom-mapping` now takes a `SocranskyReference` artifact** (TSV with a mandatory `taxon`/`complex` header) instead of `FeatureData[SocranskyMapping]`. Invalid complex names raise a clear error.
- **`delta-indices`:**
  - `--p-sample-id-col` renamed to `--p-subject-col`.
  - Output type is now `SampleData[PeriodontalIndicesDelta]` (indexed by subject ID, column `subject-id`) instead of `PeriodontalIndices`.
  - Deltas are computed for every index in the input (pMDI, pMDI_ilr, SWL, KRS).
  - Subjects with more than one sample per group, or no paired subjects at all, now raise an error instead of silently using the first sample or returning an empty table.
- **Reference map changes** (`ehomd_map.tsv` / `synonyms_map.tsv`), which change complex assignments:
  - *Actinomyces odontolyticus* / *Schaalia odontolytica*: Blue → **Purple** (Socransky et al., 1998).
  - Genus-level entries for *Porphyromonas* and *Prevotella* removed (species entries kept); unmatched species of these genera are now `Other`.
  - Wrong synonyms removed (*Campylobacter sputorum* → *C. rectus*, *Eubacterium brachy* → *Mogibacterium timidum*, *Streptococcus mitis oralis*, phylum names); added *Wolinella recta* → *C. rectus* and *Eubacterium timidum* → *M. timidum*.
- **Taxonomy parsing.** SILVA/GTDB/Greengenes2 suffixes (`Prevotella_7`, `Fusobacterium_A`, `nucleatum_D`) are stripped and placeholder names (`uncultured_bacterium`, `metagenome`, `unidentified`, …) ignored, so more features are assigned than before.
- **Packaging.** `setup.py` removed (everything is in `pyproject.toml`); the `scikit-bio` dependency was dropped.

### Added
- `--p-min-confidence` and `--p-pseudocount` parameters for `collapse-complexes`.
- `SocranskyReference` and `PeriodontalIndicesDelta` semantic types with formats and transformers.
- `source` column in `ehomd_map.tsv` (`Socransky1998`, a `citations.bib` key, or `unverified`).
- Citations: Socransky & Haffajee (2002); Pérez-Chaparro et al. (2014).
- Benjamini-Hochberg *q*-values in `summarize`.
- Format validation for complexes, indices, deltas and references.
- Test suite based on `qiime2.plugin.testing.TestPluginBase` with hand-computed expectations.
- `CHANGELOG.md`.

### Fixed
- `summarize` never displayed statistics: `numpy` was not imported and the resulting `NameError` was swallowed by a bare `except`. Failures now emit a warning with the reason.
- `summarize` only tests categorical metadata columns.
- `summarize` works offline: Vega 5.30.0, Vega-Lite 5.21.0 and Vega-Embed 6.26.0 are bundled instead of loaded from a CDN.
- Sample IDs containing `</script>` can no longer break the report.
- `delta-indices`: metadata columns named like an index no longer collide; numeric group columns work.
- Incorrect metadata of the Curtis et al. (2020), Chen et al. (2010) and Olsen & Yilmaz (2019) entries in `citations.bib`.
- Collapsing no longer densifies the feature table; indices are computed vectorized.

### Changed
- CI and `environment.yml` use the official QIIME 2 amplicon 2024.10 environment.
- `espec.md` moved to `docs/SPEC.md`.
- Author e-mail updated.
