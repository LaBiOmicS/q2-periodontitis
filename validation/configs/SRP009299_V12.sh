# Griffen et al. 2012: 454 Titanium, V1-V2; reads keep 17 nt of 27F.
# The SRA quality strings are constant (Q30 with Q3 tails: the reads were
# submitted without real base qualities), so DADA2 cannot learn an error
# model ("Error matrix is NULL"). Deblur, whose error model is static,
# is used instead (see docs/VALIDATION.md, Deviations).
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=deblur
DEBLUR_ARGS="--p-left-trim-len 17 --p-trim-length 220"
PRIMER_F=AGAGTTTGATCMTGGCTCAG
PRIMER_R=TGCTGCCTCCCGTAGGAGT
