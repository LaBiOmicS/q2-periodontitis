import pandas as pd
import biom
from q2_types.feature_table import BIOMV210Format, BIOMV210DirFmt
from .plugin_setup import plugin
from ._format import (
    SocranskyComplexesFormat,
    PeriodontalIndicesFormat,
    SocranskyMappingFormat,
)

@plugin.register_transformer
def _1(format: SocranskyComplexesFormat) -> pd.DataFrame:
    with format.open() as fh:
        return pd.read_csv(fh, sep='\t', index_col=0)

@plugin.register_transformer
def _2(df: pd.DataFrame) -> SocranskyComplexesFormat:
    format = SocranskyComplexesFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='sample-id')
    return format

@plugin.register_transformer
def _3(format: PeriodontalIndicesFormat) -> pd.DataFrame:
    with format.open() as fh:
        return pd.read_csv(fh, sep='\t', index_col=0)

@plugin.register_transformer
def _4(df: pd.DataFrame) -> PeriodontalIndicesFormat:
    format = PeriodontalIndicesFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='sample-id')
    return format

@plugin.register_transformer
def _5(format: SocranskyMappingFormat) -> pd.DataFrame:
    with format.open() as fh:
        return pd.read_csv(fh, sep='\t', index_col=0)

@plugin.register_transformer
def _6(df: pd.DataFrame) -> SocranskyMappingFormat:
    format = SocranskyMappingFormat()
    with format.open() as fh:
        df.to_csv(fh, sep='\t', index_label='Feature ID')
    return format

import qiime2

@plugin.register_transformer
def _8(format: PeriodontalIndicesFormat) -> qiime2.Metadata:
    with format.open() as fh:
        return qiime2.Metadata(pd.read_csv(fh, sep='\t', index_col=0))

@plugin.register_transformer
def _9(format: SocranskyComplexesFormat) -> qiime2.Metadata:
    with format.open() as fh:
        return qiime2.Metadata(pd.read_csv(fh, sep='\t', index_col=0))

@plugin.register_transformer
def _10(format: SocranskyMappingFormat) -> qiime2.Metadata:
    with format.open() as fh:
        return qiime2.Metadata(pd.read_csv(fh, sep='\t', index_col=0))

