import qiime2.plugin.model as model

class SocranskyComplexesFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
            if not header:
                raise model.ValidationError("File is empty.")
            # Verify TSV structure
            if "\t" not in header and not header.startswith("#"):
                raise model.ValidationError("Header line must be tab-separated.")

SocranskyComplexesDirectoryFormat = model.SingleFileDirectoryFormat(
    'SocranskyComplexesDirectoryFormat',
    'complexes.tsv',
    SocranskyComplexesFormat
)

class PeriodontalIndicesFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
            if not header:
                raise model.ValidationError("File is empty.")
            if "\t" not in header:
                raise model.ValidationError("Header line must be tab-separated.")

PeriodontalIndicesDirectoryFormat = model.SingleFileDirectoryFormat(
    'PeriodontalIndicesDirectoryFormat',
    'indices.tsv',
    PeriodontalIndicesFormat
)

class SocranskyMappingFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.open() as fh:
            header = fh.readline()
            if not header:
                raise model.ValidationError("File is empty.")
            if "\t" not in header:
                raise model.ValidationError("Header line must be tab-separated.")

SocranskyMappingDirectoryFormat = model.SingleFileDirectoryFormat(
    'SocranskyMappingDirectoryFormat',
    'mapping.tsv',
    SocranskyMappingFormat
)
