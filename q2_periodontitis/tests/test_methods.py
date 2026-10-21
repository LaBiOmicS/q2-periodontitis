import math
import os
import unittest

import numpy as np
import pandas as pd
import biom
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_periodontitis._methods import (
    COMPLEXES,
    DATA_DIR,
    parse_taxonomy_string,
    assign_feature_to_complex,
    build_lookup,
    load_default_mappings,
    collapse_complexes,
    compute_indices,
    delta_indices
)


def _taxonomy(taxa, confidence=None):
    df = pd.DataFrame({"Taxon": list(taxa.values())},
                      index=pd.Index(list(taxa.keys()), name="Feature ID"))
    if confidence is not None:
        df["Confidence"] = [confidence[k] for k in taxa]
    return df


class TestParser(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        _, self.synonyms = load_default_mappings()

    def test_silva_binomial(self):
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
                "f__Porphyromonadaceae; g__Porphyromonas; s__Porphyromonas_gingivalis"),
            ("Porphyromonas", "Porphyromonas gingivalis"))

    def test_silva_numeric_genus_suffix(self):
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
                "f__Prevotellaceae; g__Prevotella_7; s__Prevotella_melaninogenica"),
            ("Prevotella", "Prevotella melaninogenica"))
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Spirochaetota; c__Spirochaetia; o__Spirochaetales; "
                "f__Spirochaetaceae; g__Treponema_2; s__uncultured_bacterium"),
            ("Treponema", None))

    def test_gtdb_letter_suffixes(self):
        # Genus and epithet suffixes (GTDB / Greengenes2)
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Fusobacteriota; c__Fusobacteriia; o__Fusobacteriales; "
                "f__Fusobacteriaceae; g__Fusobacterium_A; s__Fusobacterium_A nucleatum_D"),
            ("Fusobacterium", "Fusobacterium nucleatum"))
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
                "f__Bacteroidaceae; g__Prevotella_D; s__Prevotella_D intermedia"),
            ("Prevotella", "Prevotella intermedia"))

    def test_greengenes_epithet_only(self):
        self.assertEqual(
            parse_taxonomy_string(
                "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; "
                "f__Porphyromonadaceae; g__Porphyromonas; s__gingivalis"),
            ("Porphyromonas", "Porphyromonas gingivalis"))

    def test_no_rank_prefixes(self):
        self.assertEqual(
            parse_taxonomy_string(
                "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;"
                "Porphyromonadaceae;Porphyromonas;gingivalis"),
            ("Porphyromonas", "Porphyromonas gingivalis"))

    def test_oral_taxon_number_is_not_a_suffix(self):
        _, species = parse_taxonomy_string(
            "d__Bacteria; g__Fretibacterium; s__Fretibacterium_human_oral_taxon_361")
        self.assertEqual(species, "Fretibacterium human oral taxon 361")

    def test_uninformative_labels(self):
        for label in ["uncultured_bacterium", "metagenome", "unidentified",
                      "human_oral_metagenome", "Ambiguous_taxa", "uncultured"]:
            genus, species = parse_taxonomy_string(
                f"d__Bacteria; g__Porphyromonas; s__{label}")
            self.assertEqual((genus, species), ("Porphyromonas", None), label)
        self.assertEqual(
            parse_taxonomy_string("d__Bacteria; g__uncultured; s__metagenome"),
            (None, None))
        self.assertEqual(parse_taxonomy_string("Unassigned"), (None, None))
        self.assertEqual(parse_taxonomy_string(""), (None, None))

    def test_synonyms(self):
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; g__Tannerella; s__Bacteroides_forsythus", self.synonyms),
            ("Tannerella", "Tannerella forsythia"))
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; g__Wolinella; s__Wolinella_recta", self.synonyms)[1],
            "Campylobacter rectus")
        self.assertEqual(
            parse_taxonomy_string(
                "d__Bacteria; g__Peptostreptococcus; s__Peptostreptococcus_micros",
                self.synonyms)[1],
            "Parvimonas micra")


class TestReferenceMap(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        tax_map, self.synonyms = load_default_mappings()
        self.lookup = build_lookup(tax_map)

    def assign(self, tax_str, mode="genus-majority"):
        return assign_feature_to_complex(tax_str, self.lookup, self.synonyms, mode)

    def test_reference_map_is_valid(self):
        ref = pd.read_csv(os.path.join(DATA_DIR, "ehomd_map.tsv"), sep="\t")
        self.assertEqual(list(ref.columns), ["taxon", "complex", "level", "source"])
        self.assertTrue(ref["complex"].isin(COMPLEXES).all())
        self.assertFalse(ref["taxon"].duplicated().any())
        self.assertFalse(ref["source"].isna().any())

    def test_odontolyticus_is_purple(self):
        self.assertEqual(self.assign("g__Actinomyces; s__Actinomyces_odontolyticus"),
                         ("Purple", "species"))
        self.assertEqual(self.assign("g__Schaalia; s__Schaalia_odontolytica"),
                         ("Purple", "species"))

    def test_heterogeneous_genera_not_mapped(self):
        # Genus-level Prevotella / Porphyromonas were removed from the map
        self.assertEqual(self.assign("g__Prevotella_7; s__Prevotella_melaninogenica"),
                         ("Other", "unassigned"))
        self.assertEqual(self.assign("g__Porphyromonas; s__Porphyromonas_catoniae"),
                         ("Other", "unassigned"))
        self.assertEqual(self.assign("g__Prevotella; s__Prevotella_intermedia"),
                         ("Orange", "species"))

    def test_synonym_resolves_to_complex(self):
        self.assertEqual(self.assign("g__Tannerella; s__Bacteroides_forsythus", "strict"),
                         ("Red", "species"))

    def test_strict_vs_genus_majority(self):
        tax = "g__Fusobacterium; s__Fusobacterium_periodonticum"
        self.assertEqual(self.assign(tax, "strict"), ("Other", "unassigned"))
        self.assertEqual(self.assign(tax, "genus-majority"), ("Orange", "genus"))

    def test_invalid_mode(self):
        with self.assertRaisesRegex(ValueError, "assignment_mode"):
            self.assign("g__Fusobacterium", "probabilistic")

    def test_build_lookup_rejects_invalid_complex(self):
        with self.assertRaisesRegex(ValueError, "Invalid complex name.*'Pink'"):
            build_lookup({"Selenomonas noxia": "Pink"})


class TestCollapseComplexes(TestPluginBase):
    package = "q2_periodontitis.tests"

    # Counts (features x samples). Relative abundances per sample:
    #   S1: Red .1, Orange .2, Yellow .3, Blue .2, Purple .1, Other .1
    #   S2: Red .005 (= theta), Yellow .6, Blue .395
    #   S3: Orange .5, Yellow .5
    #   S4: empty sample
    COUNTS = np.array([
        [10, 5, 0, 0],      # F0 P. gingivalis -> Red
        [20, 0, 50, 0],     # F1 F. nucleatum -> Orange
        [30, 600, 50, 0],   # F2 S. gordonii -> Yellow
        [20, 395, 0, 0],    # F3 A. naeslundii -> Blue
        [10, 0, 0, 0],      # F4 V. parvula -> Purple
        [10, 0, 0, 0],      # F5 L. casei -> Other
    ])
    TAXA = {
        "F0": "d__Bacteria; g__Porphyromonas; s__Porphyromonas_gingivalis",
        "F1": "d__Bacteria; g__Fusobacterium; s__Fusobacterium_nucleatum",
        "F2": "d__Bacteria; g__Streptococcus; s__Streptococcus_gordonii",
        "F3": "d__Bacteria; g__Actinomyces; s__Actinomyces_naeslundii",
        "F4": "d__Bacteria; g__Veillonella; s__Veillonella_parvula",
        "F5": "d__Bacteria; g__Lactobacillus; s__Lactobacillus_casei",
    }

    def setUp(self):
        super().setUp()
        self.table = biom.Table(self.COUNTS, list(self.TAXA),
                                ["S1", "S2", "S3", "S4"])
        self.taxonomy = _taxonomy(self.TAXA)

    def test_complexes_hand_computed(self):
        complexes, _, mapping = collapse_complexes(self.table, self.taxonomy)
        self.assertEqual(list(complexes.columns), COMPLEXES)
        exp_s1 = [0.1, 0.2, 0.0, 0.3, 0.2, 0.0, 0.1, 0.1]
        np.testing.assert_allclose(complexes.loc["S1"].values, exp_s1)
        np.testing.assert_allclose(complexes.sum(axis=1).values, [1, 1, 1, 0])
        self.assertEqual(mapping.loc["F5", "Complex"], "Other")
        self.assertEqual(mapping.loc["F5", "Resolution"], "unassigned")

    def test_indices_hand_computed(self):
        _, indices, _ = collapse_complexes(self.table, self.taxonomy)
        self.assertEqual(list(indices.columns), ["pMDI", "pMDI_ilr", "SWL", "KRS"])
        eps = 1e-5
        r = math.sqrt(6 / 5)
        ln = math.log

        # pMDI = log10((D + eps) / (H + eps)); D = Red+Orange+Extended, H = Yellow+Blue
        self.assertAlmostEqual(indices.loc["S1", "pMDI"], -0.22184295917767582, places=12)
        self.assertAlmostEqual(indices.loc["S1", "pMDI"],
                               math.log10((0.3 + eps) / (0.5 + eps)), places=12)
        self.assertAlmostEqual(indices.loc["S2", "pMDI"], -2.2979897196252086, places=12)
        self.assertAlmostEqual(indices.loc["S3", "pMDI"], 0.0, places=12)
        self.assertAlmostEqual(indices.loc["S4", "pMDI"], 0.0, places=12)

        # pMDI_ilr = sqrt(6/5) * (mean ln(D parts + eps) - mean ln(H parts + eps))
        self.assertAlmostEqual(indices.loc["S1", "pMDI_ilr"], -4.091417176548917, places=10)
        self.assertAlmostEqual(
            indices.loc["S1", "pMDI_ilr"],
            r * ((ln(0.1 + eps) + ln(0.2 + eps) + ln(eps)) / 3
                 - (ln(0.3 + eps) + ln(0.2 + eps)) / 2), places=10)
        self.assertAlmostEqual(indices.loc["S2", "pMDI_ilr"], -9.553263856876324, places=10)
        self.assertAlmostEqual(indices.loc["S3", "pMDI_ilr"], -1.9754158626392213, places=10)
        self.assertAlmostEqual(indices.loc["S4", "pMDI_ilr"], 0.0, places=12)

        # SWL = 1.0 Red + 0.9 Extended + 0.6 Orange + 0.2 Green + 0.1 Purple
        self.assertAlmostEqual(indices.loc["S1", "SWL"], 0.1 + 0.6 * 0.2 + 0.1 * 0.1)
        self.assertAlmostEqual(indices.loc["S1", "SWL"], 0.23)
        self.assertAlmostEqual(indices.loc["S2", "SWL"], 0.005)
        self.assertAlmostEqual(indices.loc["S3", "SWL"], 0.3)
        self.assertAlmostEqual(indices.loc["S4", "SWL"], 0.0)

        # KRS = (1 + tanh(alpha (Red - theta) / theta)) (Red + Extended + 0.5 Orange)
        # S1: tanh(19) ~ 1 -> 2 * (0.1 + 0.1) = 0.4
        self.assertAlmostEqual(indices.loc["S1", "KRS"], 0.4, places=12)
        # S2: Red == theta -> sigmoid factor exactly 1
        self.assertAlmostEqual(indices.loc["S2", "KRS"], 0.005, places=12)
        # S3: Red == 0 -> (1 - tanh(1)) * 0.25
        self.assertAlmostEqual(indices.loc["S3", "KRS"], 0.05960146101105879, places=12)
        self.assertAlmostEqual(indices.loc["S4", "KRS"], 0.0, places=12)

    def test_pmdi_ilr_is_not_a_rescaled_clr(self):
        # pMDI_clr (removed) was pMDI_ilr / sqrt(6/5); make sure only ilr is output
        _, indices, _ = collapse_complexes(self.table, self.taxonomy)
        self.assertNotIn("pMDI_clr", indices.columns)

    def test_pseudocount_parameter(self):
        _, indices, _ = collapse_complexes(self.table, self.taxonomy, pseudocount=0.01)
        self.assertAlmostEqual(indices.loc["S1", "pMDI"],
                               math.log10(0.31 / 0.51), places=12)
        with self.assertRaisesRegex(ValueError, "pseudocount"):
            collapse_complexes(self.table, self.taxonomy, pseudocount=0)

    def test_keystone_parameters(self):
        complexes, _, _ = collapse_complexes(self.table, self.taxonomy)
        krs = compute_indices(complexes, keystone_threshold=0.1,
                              keystone_steepness=2.0)["KRS"]
        # S1 Red == theta -> factor 1; S2 Red = 0.005 -> 1 + tanh(2 * (0.005 - 0.1) / 0.1)
        self.assertAlmostEqual(krs["S1"], 0.2, places=12)
        self.assertAlmostEqual(krs["S2"], (1 + math.tanh(-1.9)) * 0.005, places=12)

    def test_min_confidence(self):
        conf = {"F0": 0.99, "F1": 0.5, "F2": 0.95, "F3": 0.95, "F4": 0.95, "F5": 0.95}
        taxonomy = _taxonomy(self.TAXA, conf)

        complexes, _, mapping = collapse_complexes(self.table, taxonomy,
                                                   min_confidence=0.7)
        self.assertEqual(mapping.loc["F1", "Complex"], "Other")
        self.assertEqual(mapping.loc["F1", "Resolution"], "low-confidence")
        self.assertAlmostEqual(mapping.loc["F1", "Confidence"], 0.5)
        self.assertAlmostEqual(complexes.loc["S1", "Orange"], 0.0)
        self.assertAlmostEqual(complexes.loc["S1", "Other"], 0.3)
        np.testing.assert_allclose(complexes.sum(axis=1).values, [1, 1, 1, 0])

        # Default (0.0) does not filter and does not weight abundances
        complexes, _, mapping = collapse_complexes(self.table, taxonomy)
        self.assertAlmostEqual(complexes.loc["S1", "Orange"], 0.2)
        self.assertAlmostEqual(complexes.loc["S1", "Red"], 0.1)
        self.assertIn("Confidence", mapping.columns)

    def test_missing_confidence_is_not_filtered(self):
        conf = {k: "-" for k in self.TAXA}
        taxonomy = _taxonomy(self.TAXA, conf)
        _, _, mapping = collapse_complexes(self.table, taxonomy, min_confidence=0.9)
        self.assertEqual(mapping.loc["F0", "Complex"], "Red")

    def test_custom_mapping_valid(self):
        custom = pd.DataFrame({"taxon": ["Lactobacillus casei", "Streptococcus gordonii"],
                               "complex": ["Extended", "Green"],
                               "note": ["x", "y"]})
        complexes, _, mapping = collapse_complexes(self.table, self.taxonomy,
                                                   custom_mapping=custom)
        self.assertEqual(mapping.loc["F5", "Complex"], "Extended")
        self.assertEqual(mapping.loc["F2", "Complex"], "Green")
        self.assertAlmostEqual(complexes.loc["S1", "Extended"], 0.1)

    def test_custom_mapping_reads_columns_by_name(self):
        custom = pd.DataFrame({"complex": ["Extended"],
                               "taxon": ["Lactobacillus casei"]})
        _, _, mapping = collapse_complexes(self.table, self.taxonomy,
                                           custom_mapping=custom)
        self.assertEqual(mapping.loc["F5", "Complex"], "Extended")

    def test_custom_mapping_invalid_complex(self):
        custom = pd.DataFrame({"taxon": ["Lactobacillus casei"], "complex": ["red"]})
        with self.assertRaisesRegex(ValueError, "Invalid complex name.*'red'"):
            collapse_complexes(self.table, self.taxonomy, custom_mapping=custom)

    def test_custom_mapping_missing_columns(self):
        custom = pd.DataFrame({"name": ["Lactobacillus casei"], "complex": ["Red"]})
        with self.assertRaisesRegex(ValueError, "taxon"):
            collapse_complexes(self.table, self.taxonomy, custom_mapping=custom)

    def test_action_with_custom_reference_artifact(self):
        table = qiime2.Artifact.import_data("FeatureTable[Frequency]", self.table)
        taxonomy = qiime2.Artifact.import_data("FeatureData[Taxonomy]", self.taxonomy)
        reference = qiime2.Artifact.import_data(
            "SocranskyReference", self.get_data_path("reference-valid.tsv"))
        action = self.plugin.actions["collapse_complexes"]
        complexes, indices, mapping = action(table=table, taxonomy=taxonomy,
                                             custom_mapping=reference,
                                             min_confidence=0.0)
        self.assertEqual(str(complexes.type), "SocranskyComplexes")
        self.assertEqual(str(indices.type), "PeriodontalIndices")
        mapping_df = mapping.view(pd.DataFrame)
        self.assertEqual(mapping_df.loc["F5", "Complex"], "Yellow")
        complexes_df = complexes.view(pd.DataFrame)
        np.testing.assert_allclose(complexes_df.sum(axis=1).values, [1, 1, 1, 0])


class TestDeltaIndices(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        self.indices = pd.DataFrame({
            "pMDI": [1.5, -0.2, 0.5, 0.0],
            "pMDI_ilr": [2.0, 1.0, 0.0, 0.5],
            "SWL": [0.8, 0.2, 0.4, 0.4],
            "KRS": [0.9, 0.1, 0.3, 0.6]
        }, index=["a", "b", "c", "d"])

    def metadata(self, subject, group, extra=None):
        df = pd.DataFrame({"subject": subject, "visit": group},
                          index=pd.Index(["a", "b", "c", "d"][:len(subject)], name="id"))
        if extra:
            for k, v in extra.items():
                df[k] = v
        return qiime2.Metadata(df)

    def test_success(self):
        md = self.metadata(["P1", "P1", "P2", "P2"], ["Pre", "Post", "Pre", "Post"],
                           extra={"pMDI": ["x"] * 4})  # colliding column name
        res = delta_indices(self.indices, md, subject_col="subject", group_col="visit",
                            baseline_group="Pre", target_group="Post")
        self.assertEqual(list(res.index), ["P1", "P2"])
        self.assertEqual(res.index.name, "subject-id")
        self.assertEqual(list(res.columns), ["pMDI", "pMDI_ilr", "SWL", "KRS"])
        self.assertAlmostEqual(res.loc["P1", "pMDI"], -1.7)
        self.assertAlmostEqual(res.loc["P1", "pMDI_ilr"], -1.0)
        self.assertAlmostEqual(res.loc["P2", "KRS"], 0.3)

    def test_numeric_group_column(self):
        md = self.metadata(["P1", "P1", "P2", "P2"], [0, 1, 0, 1])
        self.assertEqual(md.columns["visit"].type, "numeric")
        res = delta_indices(self.indices, md, "subject", "visit", "0", "1")
        self.assertEqual(list(res.index), ["P1", "P2"])
        self.assertAlmostEqual(res.loc["P2", "pMDI"], -0.5)

    def test_unpaired_subjects_are_skipped(self):
        md = self.metadata(["P1", "P1", "P2", "P3"], ["Pre", "Post", "Pre", "Post"])
        res = delta_indices(self.indices, md, "subject", "visit", "Pre", "Post")
        self.assertEqual(list(res.index), ["P1"])

    def test_duplicates_raise(self):
        md = self.metadata(["P1", "P1", "P1", "P2"], ["Pre", "Post", "Post", "Pre"])
        with self.assertRaisesRegex(ValueError, "More than one sample.*P1"):
            delta_indices(self.indices, md, "subject", "visit", "Pre", "Post")

    def test_no_pairs_raise(self):
        md = self.metadata(["P1", "P2", "P3", "P4"], ["Pre", "Post", "Pre", "Post"])
        with self.assertRaisesRegex(ValueError, "No subject"):
            delta_indices(self.indices, md, "subject", "visit", "Pre", "Post")

    def test_missing_column_raise(self):
        md = self.metadata(["P1", "P1"], ["Pre", "Post"])
        with self.assertRaisesRegex(ValueError, "not found: patient"):
            delta_indices(self.indices, md, "patient", "visit", "Pre", "Post")

    def test_action_output_type(self):
        indices = qiime2.Artifact.import_data("PeriodontalIndices", self.indices)
        md = self.metadata(["P1", "P1", "P2", "P2"], ["Pre", "Post", "Pre", "Post"])
        delta, = self.plugin.actions["delta_indices"](
            indices=indices, metadata=md, subject_col="subject", group_col="visit",
            baseline_group="Pre", target_group="Post")
        self.assertEqual(str(delta.type), "SampleData[PeriodontalIndicesDelta]")
        self.assertEqual(list(delta.view(pd.DataFrame).index), ["P1", "P2"])


if __name__ == "__main__":
    unittest.main()
