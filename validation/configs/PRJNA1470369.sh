# Xi'an cohort, NovaSeq 2x250, V3-V4: 341F (17 nt) at the start of R1 and
# 805R (GACTACHVGGGTATCTAATCC, 21 nt) at the start of R2. Quality flat
# (median 40) along the reads; expected overlap ~30 nt: no truncation.
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
DADA2_ARGS="--p-trim-left-f 17 --p-trim-left-r 21 --p-trunc-len-f 0 --p-trunc-len-r 0"
PRIMER_F=CCTACGGGNGGCWGCAG
PRIMER_R=GACTACHVGGGTATCTAATCC
