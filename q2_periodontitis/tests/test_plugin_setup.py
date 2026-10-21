import unittest

from qiime2.plugin.testing import TestPluginBase


class TestPluginSetup(TestPluginBase):
    package = "q2_periodontitis.tests"

    def test_plugin_name(self):
        self.assertEqual(self.plugin.name, 'periodontitis')

    def test_actions_registered(self):
        self.assertIn('collapse_complexes', self.plugin.methods)
        self.assertIn('delta_indices', self.plugin.methods)
        self.assertIn('summarize', self.plugin.visualizers)

    def test_formats_registered(self):
        for fmt in ['SocranskyComplexesDirectoryFormat',
                    'PeriodontalIndicesDirectoryFormat',
                    'PeriodontalIndicesDeltaDirectoryFormat',
                    'SocranskyMappingDirectoryFormat',
                    'SocranskyReferenceDirectoryFormat']:
            self.assertIn(fmt, self.plugin.formats)

    def test_collapse_signature(self):
        sig = self.plugin.methods['collapse_complexes'].signature
        self.assertIn('min_confidence', sig.parameters)
        self.assertIn('pseudocount', sig.parameters)
        self.assertEqual(str(sig.inputs['custom_mapping'].qiime_type),
                         'SocranskyReference')
        modes = sig.parameters['assignment_mode'].qiime_type
        self.assertIn('strict', str(modes))
        self.assertNotIn('probabilistic', str(modes))

    def test_delta_signature(self):
        sig = self.plugin.methods['delta_indices'].signature
        self.assertIn('subject_col', sig.parameters)
        self.assertNotIn('sample_id_col', sig.parameters)
        self.assertEqual(str(sig.outputs['delta_table'].qiime_type),
                         'SampleData[PeriodontalIndicesDelta]')


if __name__ == "__main__":
    unittest.main()
