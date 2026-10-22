# PERIDYSSEQ (France): MiSeq 2x300, V3-V4, 341F (17 nt) at the start of R1
# and 805R (GACTACHVGGGTATCTAATCC, 21 nt) at the start of R2. The lower
# quartile of quality drops below 25 after ~240 (R1) and ~200 (R2);
# 270 + 210 leaves ~16 nt of overlap for the ~426-nt insert.
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
DADA2_ARGS="--p-trim-left-f 17 --p-trim-left-r 21 --p-trunc-len-f 270 --p-trunc-len-r 210"
PRIMER_F=CCTACGGGNGGCWGCAG
PRIMER_R=GACTACHVGGGTATCTAATCC
