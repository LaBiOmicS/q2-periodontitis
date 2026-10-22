# Spanish stage I-IV cohort: MiSeq 2x300, V3-V4 (341F/805R), primers in reads
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
# R2 quality drops after ~200 nt; ~20 nt overlap remains
DADA2_ARGS="--p-trim-left-f 17 --p-trim-left-r 21 --p-trunc-len-f 280 --p-trunc-len-r 210"
PRIMER_F=CCTACGGGNGGCWGCAG
PRIMER_R=GACTACHVGGGTATCTAATCC
