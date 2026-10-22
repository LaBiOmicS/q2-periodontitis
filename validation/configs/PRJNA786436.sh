# Byrne et al. 2022: Ion Torrent PGM, V4 from 515F (primer in the reads,
# ~292 nt); DADA2 pyro settings are the recommended ones for Ion Torrent
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=pyro
DADA2_ARGS="--p-trim-left 19 --p-trunc-len 250"
PRIMER_F=GTGYCAGCMGCCGCGGTAA
PRIMER_R=GGACTACNVGGGTWTCTAAT
