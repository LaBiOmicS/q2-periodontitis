import unittest
import numpy as np
import pandas as pd
import biom
import qiime2
from q2_periodontitis._methods import (
    parse_taxonomy_string,
    assign_feature_to_complex,
    build_lookup,
    collapse_complexes,
    delta_indices
)

class TestMethods(unittest.TestCase):

    def setUp(self):
        self.tax_map = {
            "Porphyromonas gingivalis": "Red",
            "Porphyromonas": "Red",
            "Fusobacterium nucleatum": "Orange",
            "Streptococcus gordonii": "Yellow",
            "Actinomyces naeslundii": "Blue"
        }
        self.synonyms = {
            "Bacteroides forsythus": "Tannerella forsythia",
            "Peptostreptococcus micros": "Parvimonas micra"
        }

    def test_parse_taxonomy_string(self):
        genus, species = parse_taxonomy_string(
            "d__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae;g__Porphyromonas;s__Porphyromonas_gingivalis",
            self.synonyms
        )
        self.assertEqual(genus, "Porphyromonas")
        self.assertEqual(species, "Porphyromonas gingivalis")

    def test_assign_feature_to_complex_species(self):
        complex_name, res = assign_feature_to_complex(
            "d__Bacteria;...;g__Porphyromonas;s__Porphyromonas_gingivalis",
            build_lookup(self.tax_map),
            self.synonyms,
            assignment_mode="strict"
        )
        self.assertEqual(complex_name, "Red")
        self.assertEqual(res, "species")

    def test_assign_feature_to_complex_genus_majority(self):
        complex_name, res = assign_feature_to_complex(
            "d__Bacteria;...;g__Porphyromonas;s__uncultured",
            build_lookup(self.tax_map),
            self.synonyms,
            assignment_mode="genus-majority"
        )
        self.assertEqual(complex_name, "Red")
        self.assertEqual(res, "genus")

    def test_collapse_complexes(self):
        # Create dummy BIOM table: 3 features x 2 samples
        data = np.array([
            [100, 10],   # ASV1: P. gingivalis (Red)
            [50, 100],   # ASV2: F. nucleatum (Orange)
            [10, 200]    # ASV3: S. gordonii (Yellow)
        ])
        table = biom.Table(data, observation_ids=["ASV1", "ASV2", "ASV3"], sample_ids=["S1", "S2"])
        taxonomy = pd.Series({
            "ASV1": "d__Bacteria;...;g__Porphyromonas;s__Porphyromonas_gingivalis",
            "ASV2": "d__Bacteria;...;g__Fusobacterium;s__Fusobacterium_nucleatum",
            "ASV3": "d__Bacteria;...;g__Streptococcus;s__Streptococcus_gordonii"
        })

        df_complexes, df_indices, df_mapping = collapse_complexes(
            table, taxonomy, assignment_mode="genus-majority"
        )

        self.assertIn("Red", df_complexes.columns)
        self.assertIn("pMDI", df_indices.columns)
        self.assertIn("pMDI_ilr", df_indices.columns)
        self.assertIn("SWL", df_indices.columns)
        self.assertIn("KRS", df_indices.columns)
        self.assertEqual(len(df_complexes), 2)
        self.assertEqual(len(df_indices), 2)

    def test_collapse_complexes_with_confidence_dataframe(self):
        data = np.array([
            [100, 10],
            [50, 100]
        ])
        table = biom.Table(data, observation_ids=["ASV1", "ASV2"], sample_ids=["S1", "S2"])
        taxonomy = pd.DataFrame({
            "Taxon": [
                "d__Bacteria;...;g__Porphyromonas;s__Porphyromonas_gingivalis",
                "d__Bacteria;...;g__Fusobacterium;s__Fusobacterium_nucleatum"
            ],
            "Confidence": [0.98, 0.85]
        }, index=["ASV1", "ASV2"])

        df_complexes, df_indices, df_mapping = collapse_complexes(
            table, taxonomy, assignment_mode="genus-majority", keystone_threshold=0.005, keystone_steepness=1.0
        )

        self.assertIn("pMDI_ilr", df_indices.columns)
        self.assertIn("Confidence", df_mapping.columns)
        self.assertEqual(df_mapping.loc["ASV1", "Confidence"], 0.98)

    def test_delta_indices(self):
        indices = pd.DataFrame({
            "pMDI": [1.5, -0.2],
            "pMDI_ilr": [1.2, -0.1],
            "SWL": [0.8, 0.2],
            "KRS": [0.9, 0.1]
        }, index=["S1_Pre", "S1_Post"])

        metadata = qiime2.Metadata(pd.DataFrame({
            "Subject_ID": ["Pat1", "Pat1"],
            "TimePoint": ["Pre", "Post"]
        }, index=pd.Index(["S1_Pre", "S1_Post"], name="id")))

        df_delta = delta_indices(
            indices, metadata,
            subject_col="Subject_ID",
            group_col="TimePoint",
            baseline_group="Pre",
            target_group="Post"
        )

        self.assertEqual(len(df_delta), 1)
        self.assertIn("pMDI", df_delta.columns)
        self.assertAlmostEqual(df_delta.loc["Pat1", "pMDI"], -0.2 - 1.5)

if __name__ == "__main__":
    unittest.main()
