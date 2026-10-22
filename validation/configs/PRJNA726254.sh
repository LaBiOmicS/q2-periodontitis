# Houston cohort: MiSeq, only R1 is public (301 nt). Reads start with a
# 14-nt linker and 806R (V4 read from the reverse primer); the lower
# quartile of quality drops below 30 after position 240 (read_qc.py).
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=single
DADA2_ARGS="--p-trim-left 34 --p-trunc-len 240"
PRIMER_F=GTGYCAGCMGCCGCGGTAA
PRIMER_R=GGACTACNVGGGTWTCTAAT
