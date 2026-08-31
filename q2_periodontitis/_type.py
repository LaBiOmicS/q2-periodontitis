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

SocranskyMapping = SemanticType(
    'SocranskyMapping',
    variant_of=FeatureData.field['type']
)
