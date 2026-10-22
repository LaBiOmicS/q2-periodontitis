# Elderly periodontitis, NovaSeq X 2x250, V3-V4 (341F/806R), primers at
# the read starts (ending at position 20 in R1 and R2). Quality is binned
# (4 values) and flat along the reads; the expected overlap after primer
# removal is ~30 nt, so the reads are not truncated (read_qc.py).
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
DADA2_ARGS="--p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 0 --p-trunc-len-r 0"
PRIMER_F=CCTACGGGNGGCWGCAG
PRIMER_R=GGACTACHVGGGTWTCTAAT
