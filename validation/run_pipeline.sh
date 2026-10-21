#!/usr/bin/env bash
# Reprocesses the subgingival plaque samples of Chen et al. (2018, ISME J,
# SRA SRP075100) from raw reads to q2-periodontitis outputs.
#
#   bash validation/run_pipeline.sh <work_dir>
#
# Requires an activated QIIME 2 amplicon environment (2024.10 was used) with
# q2-periodontitis installed, plus curl. Downloads ~540 MB of FASTQ.
set -euo pipefail

WORK=${1:?usage: run_pipeline.sh <work_dir>}
HERE=$(cd "$(dirname "$0")" && pwd)
THREADS=${THREADS:-4}
EHOMD=https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.23
mkdir -p "$WORK"/fastq "$WORK"/ehomd
cd "$WORK"

# 1. Run list: subgingival plaque of healthy (HP), periodontitis before (D1P)
#    and after (D2P) non-surgical therapy.
if [ ! -s runs.tsv ]; then
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP075100&result=read_run&fields=run_accession,sample_accession,sample_alias,sample_title,library_layout,instrument_model,read_count,base_count,fastq_ftp,fastq_bytes&format=tsv" -o runs.tsv
fi
awk -F'\t' 'NR>1 && $4 ~ /^(HP|D1P|D2P)_/' runs.tsv > plaque.tsv
(cd fastq && cut -f9 ../plaque.tsv | sed 's#^#https://#' \
   | xargs -P 8 -n 1 curl -s -C - -O --retry 3)

python "$HERE"/make_manifest.py plaque.tsv fastq manifest.tsv metadata.tsv

# 2. Reads are already merged, primer-free V4 (515F/806R, ~253 nt).
qiime tools import --type 'SampleData[SequencesWithQuality]' \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path manifest.tsv --output-path demux.qza
qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 --p-trunc-len 250 --p-n-threads "$THREADS" \
  --o-table table.qza --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats dada2-stats.qza

# 3. Naive Bayes classifier on eHOMD 16S RefSeq V15.23, trimmed to V4.
curl -s -o ehomd/refseq.fasta "$EHOMD"/HOMD_16S_rRNA_RefSeq_V15.23.p9.fasta
curl -s -o ehomd/taxonomy.tsv "$EHOMD"/HOMD_16S_rRNA_RefSeq_V15.23.qiime.taxonomy
qiime tools import --type 'FeatureData[Sequence]' \
  --input-path ehomd/refseq.fasta --output-path ehomd/refseq.qza
qiime tools import --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ehomd/taxonomy.tsv --output-path ehomd/taxonomy.qza
qiime feature-classifier extract-reads --i-sequences ehomd/refseq.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-n-jobs "$THREADS" --o-reads ehomd/refseq-v4.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ehomd/refseq-v4.qza \
  --i-reference-taxonomy ehomd/taxonomy.qza \
  --o-classifier ehomd/classifier-v4.qza
qiime feature-classifier classify-sklearn --i-classifier ehomd/classifier-v4.qza \
  --i-reads rep-seqs.qza --p-n-jobs "$THREADS" --o-classification taxonomy.qza

# 4. q2-periodontitis
qiime periodontitis collapse-complexes --i-table table.qza \
  --i-taxonomy taxonomy.qza --o-complex-table complexes.qza \
  --o-indices indices.qza --o-mapping mapping.qza
qiime periodontitis summarize --i-complexes complexes.qza \
  --i-indices indices.qza --i-mapping mapping.qza \
  --m-metadata-file metadata.tsv --p-group-column status \
  --p-subject-column subject --o-visualization summary.qzv
qiime periodontitis delta-indices --i-indices indices.qza \
  --m-metadata-file metadata.tsv --p-subject-col site \
  --p-group-col group --p-baseline-group D1P --p-target-group D2P \
  --o-delta-table deltas.qza

for a in complexes indices mapping deltas; do
  qiime tools export --input-path $a.qza --output-path export/$a
done
qiime tools export --input-path dada2-stats.qza --output-path export/dada2-stats

python "$HERE"/evaluate.py export metadata.tsv "$HERE"/results
