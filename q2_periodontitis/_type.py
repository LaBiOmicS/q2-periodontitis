from qiime2.plugin import SemanticType
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData

SocranskyComplexes = SemanticType(
    'SocranskyComplexes',
    variant_of=SampleData.field['type']
)

PeriodontalIndices = SemanticType(
    'PeriodontalIndices',
    variant_of=SampleData.field['type']
)

# Paired deltas (target - baseline); rows are subjects, not samples
PeriodontalIndicesDelta = SemanticType(
    'PeriodontalIndicesDelta',
    variant_of=SampleData.field['type']
)

SocranskyMapping = SemanticType(
    'SocranskyMapping',
    variant_of=FeatureData.field['type']
)

# User-supplied taxon -> complex reference table (not tied to feature IDs)
SocranskyReference = SemanticType('SocranskyReference')

# Taxa whose CLR abundances define a dysbiosis index (SMDI): taxon,
# direction (dysbiotic / normobiotic) and level (species / genus)
DysbiosisSignature = SemanticType('DysbiosisSignature')
