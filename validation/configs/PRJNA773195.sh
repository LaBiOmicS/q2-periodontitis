# Santiago/Oporto, MiSeq 2x300, V1-V3 (27F / 533R TTACCGCGGCKGCTGGCAC) on
# the same samples as PRJNA773202. Primers follow variable-length spacers:
# removed with cutadapt. The ~487-nt insert needs >= ~500 nt of the two
# trimmed reads, where quality is already low (lower quartile ~16 at R1
# 270 and ~15 at R2 235 raw positions).
# Pre-specified fallback (written before any index was computed): if the
# median fraction of reads retained by DADA2 is below 10 %, the forward
# reads alone are used (V1-V2 part of the amplicon), as a deviation.
MANIFEST_FORMAT=PairedEndFastqManifestPhred33V2
IMPORT_TYPE='SampleData[PairedEndSequencesWithQuality]'
DENOISE=paired
CUTADAPT_ARGS="--p-front-f AGAGTTTGATCMTGGCTCAG --p-front-r TTACCGCGGCKGCTGGCAC --p-discard-untrimmed"
DADA2_ARGS="--p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 270 --p-trunc-len-r 235"
PRIMER_F=AGAGTTTGATCMTGGCTCAG
PRIMER_R=TTACCGCGGCKGCTGGCAC
