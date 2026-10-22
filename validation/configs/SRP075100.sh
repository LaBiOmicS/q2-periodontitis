# Chen et al. 2018: MiSeq, merged primer-free V4 reads (515F/806R)
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=single
DADA2_ARGS="--p-trim-left 0 --p-trunc-len 250"
PRIMER_F=GTGYCAGCMGCCGCGGTAA
PRIMER_R=GGACTACNVGGGTWTCTAAT
