from qiime2.plugin import (
    Plugin,
    Citations,
    Str,
    Float,
    Choices,
    Metadata
)
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy
from ._version import __version__
from ._type import SocranskyComplexes, PeriodontalIndices, SocranskyMapping
from ._format import (
    SocranskyComplexesDirectoryFormat,
    PeriodontalIndicesDirectoryFormat,
    SocranskyMappingDirectoryFormat
)
from ._methods import collapse_complexes, delta_indices
from ._visualizer import summarize

citations = Citations.load('citations.bib', package='q2_periodontitis')

plugin = Plugin(
    name='periodontitis',
    version=__version__,
    website='https://github.com/LaBiOmicS/q2-periodontitis',
    package='q2_periodontitis',
    citations=[citations['socransky1998microbial'], citations['hajishengallis2012keystone'], citations['curtis2020periodontal']],
    description='QIIME 2 plugin for subgingival biofilm ecological characterization, Socransky complexes collapse, and dysbiosis indices.',
    short_description='Plugin for Periodontal Biofilm & Dysbiosis Characterization.'
)

# Register Directory Formats
plugin.register_formats(
    SocranskyComplexesDirectoryFormat,
    PeriodontalIndicesDirectoryFormat,
    SocranskyMappingDirectoryFormat
)

# Register Semantic Type to Format Mappings
plugin.register_semantic_type_to_format(
    SocranskyComplexes,
    artifact_format=SocranskyComplexesDirectoryFormat
)

plugin.register_semantic_type_to_format(
    PeriodontalIndices,
    artifact_format=PeriodontalIndicesDirectoryFormat
)

plugin.register_semantic_type_to_format(
    SocranskyMapping,
    artifact_format=SocranskyMappingDirectoryFormat
)

# Register Methods
plugin.methods.register_function(
    function=collapse_complexes,
    inputs={
        'table': FeatureTable[Frequency],
        'taxonomy': FeatureData[Taxonomy],
        'custom_mapping': FeatureData[SocranskyMapping]
    },
    parameters={
        'assignment_mode': Str % Choices(['strict', 'genus-majority', 'probabilistic']),
        'keystone_threshold': Float,
        'keystone_steepness': Float
    },
    outputs=[
        ('complex_table', SocranskyComplexes),
        ('indices', PeriodontalIndices),
        ('mapping', SocranskyMapping)
    ],
    input_descriptions={
        'table': 'Feature table of ASV/OTU counts.',
        'taxonomy': 'Taxonomic classifications for features.',
        'custom_mapping': 'Optional custom taxonomy-to-complex mapping table.'
    },
    parameter_descriptions={
        'assignment_mode': 'Taxonomic assignment mode for resolving genus-level ambiguities.',
        'keystone_threshold': 'Critical abundance threshold (theta) for Red complex keystone pathogens (default: 0.005).',
        'keystone_steepness': 'Steepness factor (alpha) for the sigmoidal Keystone Risk Score (default: 1.0).'
    },
    output_descriptions={
        'complex_table': 'Collapsed Socransky & Extended complexes relative abundance matrix.',
        'indices': 'Calculated dysbiosis indices (pMDI, pMDI_clr, SWL, KRS).',
        'mapping': 'Feature ID to Socransky complex mapping table.'
    },
    name='Collapse ASVs into Socransky Complexes & Compute Dysbiosis Indices',
    description='Maps taxonomic strings to Socransky/Extended complexes, collapses feature table, and calculates dysbiosis metrics.'
)

plugin.methods.register_function(
    function=delta_indices,
    inputs={
        'indices': PeriodontalIndices
    },
    parameters={
        'metadata': Metadata,
        'sample_id_col': Str,
        'group_col': Str,
        'baseline_group': Str,
        'target_group': Str
    },
    outputs=[
        ('delta_table', PeriodontalIndices)
    ],
    input_descriptions={
        'indices': 'Calculated dysbiosis indices table.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata containing subject IDs and treatment groups.',
        'sample_id_col': 'Metadata column containing subject/patient identifiers.',
        'group_col': 'Metadata column containing treatment/time-point group identifiers.',
        'baseline_group': 'Baseline group name (e.g. Pre-treatment / Baseline).',
        'target_group': 'Target group name (e.g. Post-treatment / Follow-up).'
    },
    output_descriptions={
        'delta_table': 'Table of longitudinal deltas (Target minus Baseline) for pMDI, SWL, KRS.'
    },
    name='Longitudinal Treatment Delta Calculation',
    description='Calculates paired subject deltas for dysbiosis metrics between baseline and follow-up groups.'
)

# Register Visualizers
plugin.visualizers.register_function(
    function=summarize,
    inputs={
        'complexes': SocranskyComplexes,
        'indices': PeriodontalIndices,
        'mapping': SocranskyMapping
    },
    parameters={
        'metadata': Metadata
    },
    input_descriptions={
        'complexes': 'Collapsed Socransky complexes table.',
        'indices': 'Calculated dysbiosis indices table.',
        'mapping': 'Feature resolution mapping table.'
    },
    parameter_descriptions={
        'metadata': 'Sample metadata for non-parametric group comparisons.'
    },
    name='Summarize Periodontal Complexes & Dysbiosis Metrics',
    description='Generates interactive visualization dashboard (.qzv) with statistical group comparisons.'
)
