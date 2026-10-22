# Bizzarro et al. 2016: 454 Titanium, V5-V7. Reads start with the forward
# primer GGATTAGATACCCBDGTAGTC (21 nt) and the reverse complement of the
# reverse primer (ACGTCATCCCCDCCTTCC) appears at ~390 nt; real (variable)
# base qualities, lower quartile >= 34 up to position 370 (read_qc.py).
MANIFEST_FORMAT=SingleEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[SequencesWithQuality]'
DENOISE=pyro
DADA2_ARGS="--p-trim-left 21 --p-trunc-len 370"
PRIMER_F=GGATTAGATACCCBDGTAGTC
PRIMER_R=ACGTCATCCCCDCCTTCC
