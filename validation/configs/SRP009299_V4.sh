# Griffen et al. 2012: 454 Titanium, V4; primer-free reads (~238 nt).
# Constant SRA quality strings: Deblur instead of DADA2 (see SRP009299_V12.sh).
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=deblur
DEBLUR_ARGS="--p-left-trim-len 0 --p-trim-length 220"
PRIMER_F=GTGYCAGCMGCCGCGGTAA
PRIMER_R=GGACTACNVGGGTWTCTAAT
