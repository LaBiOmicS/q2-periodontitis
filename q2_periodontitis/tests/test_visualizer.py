import os
import re
import unittest
import warnings

import numpy as np
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_periodontitis._methods import COMPLEXES
from q2_periodontitis._visualizer import summarize, benjamini_hochberg


class TestBenjaminiHochberg(unittest.TestCase):

    def test_hand_computed(self):
        # p sorted: 0.01, 0.02, 0.03, 0.04 (n=4) -> p*n/rank = 0.04, 0.04, 0.04, 0.04
        np.testing.assert_allclose(benjamini_hochberg([0.04, 0.01, 0.03, 0.02]),
                                   [0.04, 0.04, 0.04, 0.04])
        # 0.01*3/1=0.03, 0.04*3/2=0.06, 0.5*3/3=0.5
        np.testing.assert_allclose(benjamini_hochberg([0.5, 0.01, 0.04]),
                                   [0.5, 0.03, 0.06])
        # monotonicity: 0.03*2/1=0.06 > 0.04*2/2=0.04 -> both 0.04
        np.testing.assert_allclose(benjamini_hochberg([0.03, 0.04]), [0.04, 0.04])
        self.assertEqual(len(benjamini_hochberg([])), 0)


class TestSummarize(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        rng = np.random.default_rng(42)
        # One sample ID tries to close the <script> tag
        self.ids = ["S0</script><b>x"] + [f"S{i}" for i in range(1, 12)]
        n = len(self.ids)
        self.indices = pd.DataFrame(rng.random((n, 4)), index=self.ids,
                                    columns=["pMDI", "pMDI_ilr", "SWL", "KRS"])
        self.indices.iloc[6:, 0] += 10  # strong group effect on pMDI
        complexes = rng.random((n, len(COMPLEXES)))
        self.complexes = pd.DataFrame(complexes / complexes.sum(axis=1, keepdims=True),
                                      index=self.ids, columns=COMPLEXES)
        self.mapping = pd.DataFrame({"Complex": ["Red", "Other"],
                                     "Resolution": ["species", "unassigned"]},
                                    index=pd.Index(["F1", "F2"], name="Feature ID"))
        self.metadata = qiime2.Metadata(pd.DataFrame({
            "status": ["health"] * 6 + ["periodontitis"] * 6,
            "age": np.arange(n, dtype=float),
        }, index=pd.Index(self.ids, name="id")))

    def render(self, metadata=None, mapping=None):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summarize(self.temp_dir.name, self.complexes, self.indices,
                      mapping, metadata)
        with open(os.path.join(self.temp_dir.name, "index.html")) as fh:
            return fh.read()

    def stats_table(self, html):
        match = re.search(r"Non-Parametric Statistical Group Comparisons(.*?)</table>",
                          html, re.S)
        return match.group(1) if match else None

    def test_generates_index_html_without_metadata(self):
        html = self.render(mapping=self.mapping)
        self.assertIsNone(self.stats_table(html))
        for f in ["index.html", "css/style.css", "data/complexes.tsv", "data/indices.tsv"]:
            self.assertTrue(os.path.exists(os.path.join(self.temp_dir.name, f)), f)

    def test_stats_table_with_categorical_metadata(self):
        html = self.render(self.metadata)
        table = self.stats_table(html)
        self.assertIsNotNone(table)
        self.assertIn("q-value", table)
        self.assertIn("status", table)
        self.assertIn("Mann-Whitney U", table)
        # Numeric metadata columns are not tested
        self.assertNotIn("<td>age</td>", table)
        # 4 metrics x 1 categorical column
        self.assertEqual(table.count("<td>status</td>"), 4)
        # the pMDI difference is significant after FDR correction
        self.assertIn("#27ae60", table)

    def test_kruskal_for_more_groups(self):
        md = qiime2.Metadata(pd.DataFrame(
            {"site": ["a", "b", "c"] * 4}, index=pd.Index(self.ids, name="id")))
        table = self.stats_table(self.render(md))
        self.assertIn("Kruskal-Wallis", table)

    def test_skipped_test_warns(self):
        self.indices["SWL"] = 1.0  # identical values -> Kruskal raises
        md = qiime2.Metadata(pd.DataFrame(
            {"site": ["a", "b", "c"] * 4}, index=pd.Index(self.ids, name="id")))
        with self.assertWarnsRegex(UserWarning, "Skipping test of 'SWL'"):
            summarize(self.temp_dir.name, self.complexes, self.indices, None, md)

    def test_javascript_is_local(self):
        html = self.render(self.metadata)
        self.assertNotIn("cdn.jsdelivr.net", html)
        self.assertNotRegex(html, r'<script[^>]+src="https?://')
        for js in ["vega.min.js", "vega-lite.min.js", "vega-embed.min.js"]:
            self.assertIn(f'src="js/{js}"', html)
            self.assertTrue(os.path.getsize(os.path.join(self.temp_dir.name, "js", js)) > 0)

    def test_script_json_is_escaped(self):
        html = self.render()
        script = html[html.rindex("<script>"):]
        self.assertNotIn("S0</script>", script)
        self.assertIn("S0<\\/script>", script)
        # The HTML preview table escapes the ID
        self.assertIn("S0&lt;/script&gt;&lt;b&gt;x", html)
        # Only the closing tags of the real <script> elements remain
        self.assertEqual(html.count("</script>"), 4)

    def test_color_map_from_python(self):
        html = self.render()
        self.assertIn('"Red": "#E41A1C"', html)
        self.assertIn("Object.keys(colorMap)", html)


if __name__ == "__main__":
    unittest.main()
