import qiime2.plugin.model as model

from ._methods import COMPLEXES

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


class SocranskyReferenceFormat(model.TextFileFormat):
    """Tab-separated taxon -> complex table with a mandatory header that
    contains (at least) the columns ``taxon`` and ``complex``."""

    def _validate_(self, level):
        max_lines = {'min': 100, 'max': None}[level]
        with self.open() as fh:
            header = fh.readline().rstrip('\r\n')
            if not header:
                raise model.ValidationError("File is empty.")
            columns = [c.strip() for c in header.split('\t')]
            missing = [c for c in ('taxon', 'complex') if c not in columns]
            if missing:
                raise model.ValidationError(
                    "Header must be tab-separated and contain the columns "
                    "'taxon' and 'complex'; missing: %s." % ', '.join(missing))
            taxon_idx = columns.index('taxon')
            complex_idx = columns.index('complex')

            for line_no, line in enumerate(fh, start=2):
                if max_lines is not None and line_no > max_lines + 1:
                    break
                line = line.rstrip('\r\n')
                if not line.strip():
                    continue
                fields = line.split('\t')
                if len(fields) != len(columns):
                    raise model.ValidationError(
                        "Line %d has %d fields, expected %d."
                        % (line_no, len(fields), len(columns)))
                if not fields[taxon_idx].strip():
                    raise model.ValidationError(
                        "Line %d has an empty taxon." % line_no)
                complex_name = fields[complex_idx].strip()
                if complex_name not in COMPLEXES:
                    raise model.ValidationError(
                        "Line %d: invalid complex %r. Valid complexes: %s."
                        % (line_no, complex_name, ', '.join(COMPLEXES)))

SocranskyReferenceDirectoryFormat = model.SingleFileDirectoryFormat(
    'SocranskyReferenceDirectoryFormat',
    'reference.tsv',
    SocranskyReferenceFormat
)
