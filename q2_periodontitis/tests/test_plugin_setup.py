import unittest
from q2_periodontitis.plugin_setup import plugin

class TestPluginSetup(unittest.TestCase):

    def test_plugin_name(self):
        self.assertEqual(plugin.name, 'periodontitis')

    def test_actions_registered(self):
        self.assertIn('collapse_complexes', plugin.methods)
        self.assertIn('delta_indices', plugin.methods)
        self.assertIn('summarize', plugin.visualizers)

    def test_formats_registered(self):
        self.assertIn('SocranskyComplexesDirectoryFormat', plugin.formats)
        self.assertIn('PeriodontalIndicesDirectoryFormat', plugin.formats)
        self.assertIn('SocranskyMappingDirectoryFormat', plugin.formats)

if __name__ == "__main__":
    unittest.main()
