# Thai cohort, NovaSeq 2x250 trimmed by the submitters to 225 (R1) and
# 222 nt (R2); primers already removed: R1 starts just after 341F
# (TGGGGAATATTG...) and R2 just after 806R (CCTGTTCGCTCC...), so V3-V4.
# Binned quality (4 values), flat. Expected overlap ~18 nt: no truncation.
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
DADA2_ARGS="--p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 0 --p-trunc-len-r 0"
PRIMER_F=CCTACGGGNGGCWGCAG
PRIMER_R=GGACTACHVGGGTWTCTAAT
