# Pei et al. 2020: MiSeq 2x250, V3-V4 (338F/806R), primers in the reads
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
# ~25 nt overlap after primer removal: R1 is not truncated
DADA2_ARGS="--p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 0 --p-trunc-len-r 245"
PRIMER_F=ACTCCTACGGGAGGCAGCAG
PRIMER_R=GGACTACHVGGGTWTCTAAT
