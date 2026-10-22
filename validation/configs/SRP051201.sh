# Kirst et al. 2015: 454 Titanium, V1-V3 read from 534R (reverse
# orientation); reads start with an 8-nt barcode + 17-nt 534R primer
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=pyro
DADA2_ARGS="--p-trim-left 25 --p-trunc-len 400"
PRIMER_F=AGAGTTTGATCMTGGCTCAG
PRIMER_R=ATTACCGCGGCTGCTGG
