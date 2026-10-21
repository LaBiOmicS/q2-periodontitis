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
                    'SocranskyReferenceDirectoryFormat',
                    'DysbiosisSignatureDirectoryFormat']:
            self.assertIn(fmt, self.plugin.formats)

    def test_collapse_signature(self):
        sig = self.plugin.methods['collapse_complexes'].signature
        self.assertIn('min_confidence', sig.parameters)
        self.assertEqual(str(sig.inputs['custom_mapping'].qiime_type),
                         'SocranskyReference')
        self.assertEqual(str(sig.inputs['smdi_signature'].qiime_type),
                         'DysbiosisSignature')
        self.assertNotIn('assignment_mode', sig.parameters)
        self.assertNotIn('keystone_threshold', sig.parameters)
        self.assertNotIn('pseudocount', sig.parameters)

    def test_summarize_signature(self):
        sig = self.plugin.visualizers['summarize'].signature
        self.assertIn('group_column', sig.parameters)
        self.assertIn('subject_column', sig.parameters)

    def test_delta_signature(self):
        sig = self.plugin.methods['delta_indices'].signature
        self.assertIn('subject_col', sig.parameters)
        self.assertNotIn('sample_id_col', sig.parameters)
        self.assertEqual(str(sig.outputs['delta_table'].qiime_type),
                         'SampleData[PeriodontalIndicesDelta]')


if __name__ == "__main__":
    unittest.main()
