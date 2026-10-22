# Santiago/Oporto, MiSeq 2x300, V3-V4 (341F / 805R). Primers follow
# variable-length spacers, so they are removed with cutadapt (reads
# without both primers discarded). The lower quartile of quality drops
# below ~25 after ~260 (R1) and ~230 (R2) raw positions; 255 + 190 trimmed
# nt leave ~19 nt of overlap for the ~426-nt insert.
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
CUTADAPT_ARGS="--p-front-f CCTACGGGNGGCWGCAG --p-front-r GACTACHVGGGTATCTAATCC --p-discard-untrimmed"
DADA2_ARGS="--p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 255 --p-trunc-len-r 190"
PRIMER_F=CCTACGGGNGGCWGCAG
PRIMER_R=GACTACHVGGGTATCTAATCC
