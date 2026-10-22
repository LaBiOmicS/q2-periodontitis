#!/usr/bin/env bash
# Processes one 16S cohort from its QIIME 2 manifest to q2-periodontitis
# outputs, with two taxonomic classifiers:
#
#   nb       naive Bayes on eHOMD V15.23 trimmed to the cohort's primers
#            (primary; the standard QIIME 2 workflow)
#   vsearch  consensus of the top global-alignment hits at >= 98 % identity
#            against full-length eHOMD (closest to the BLAST >= 98 % rule
#            used to derive the SMDI; ties between species fall back to the
#            genus)
#
#   bash process_dataset.sh <config> <cohort_dir>
#
# <cohort_dir> holds manifest.tsv (written by prepare.py); outputs go to
# <cohort_dir>/out. The config is a shell file that sets:
#   MANIFEST_FORMAT, IMPORT_TYPE  (qiime tools import)
#   DENOISE         paired | single | pyro | deblur
#   DADA2_ARGS      trimming/truncation for DADA2 (paired, single, pyro)
#   CUTADAPT_ARGS   optional: primers removed with cutadapt trim-paired or
#                   trim-single before denoising (for primers preceded by
#                   variable-length spacers), e.g. "--p-front-f ... --p-front-r
#                   ... --p-discard-untrimmed"
#   DEBLUR_ARGS     trimming for Deblur (deblur): reads are first filtered
#                   with quality-filter q-score, then denoised with
#                   deblur denoise-16S
#   PRIMER_F, PRIMER_R  primers used to trim the reference for naive Bayes
# Environment: THREADS (default 4), REF_DIR (default <cohort_dir>/../ehomd).
# Requires an activated QIIME 2 amplicon environment with q2-periodontitis.
# Every step is skipped when its output exists, except q2-periodontitis,
# which always runs so that the outputs match the installed plugin.
set -euo pipefail

CONFIG=$(readlink -f "${1:?usage: process_dataset.sh <config> <cohort_dir>}")
COHORT=$(readlink -f "${2:?usage: process_dataset.sh <config> <cohort_dir>}")
THREADS=${THREADS:-4}
DADA2_ARGS="" DEBLUR_ARGS="" CUTADAPT_ARGS=""
# shellcheck disable=SC1090
source "$CONFIG"
EHOMD=https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/V15.23
REF=$(readlink -f "${REF_DIR:-$COHORT/../ehomd}")
WORK=$COHORT/out
mkdir -p "$WORK" "$REF"
cd "$WORK"

log() { echo "$(date +%H:%M:%S) [$(basename "$COHORT")] $*"; }

# Shared full-length eHOMD reference
if [ ! -s "$REF/refseq.qza" ]; then
  log "eHOMD V15.23"
  curl -sf -o "$REF/refseq.fasta" "$EHOMD/HOMD_16S_rRNA_RefSeq_V15.23.p9.fasta"
  curl -sf -o "$REF/taxonomy.tsv" "$EHOMD/HOMD_16S_rRNA_RefSeq_V15.23.qiime.taxonomy"
  qiime tools import --type 'FeatureData[Sequence]' \
    --input-path "$REF/refseq.fasta" --output-path "$REF/refseq.qza"
  qiime tools import --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat \
    --input-path "$REF/taxonomy.tsv" --output-path "$REF/taxonomy.qza"
fi

if [ ! -s demux.qza ]; then
  log "import"
  qiime tools import --type "$IMPORT_TYPE" --input-format "$MANIFEST_FORMAT" \
    --input-path "$COHORT/manifest.tsv" --output-path demux.qza
fi

READS=demux.qza
if [ -n "$CUTADAPT_ARGS" ]; then
  READS=demux-trimmed.qza
  if [ ! -s "$READS" ]; then
    log "cutadapt"
    mode=single
    [ "$DENOISE" = paired ] && mode=paired
    # shellcheck disable=SC2086
    qiime cutadapt trim-$mode --i-demultiplexed-sequences demux.qza \
      $CUTADAPT_ARGS --p-cores "$THREADS" \
      --o-trimmed-sequences "$READS" > /dev/null
  fi
fi

if [ ! -s table.qza ]; then
  log "denoise ($DENOISE)"
  if [ "$DENOISE" = deblur ]; then
    qiime quality-filter q-score --i-demux "$READS" \
      --o-filtered-sequences demux-filtered.qza \
      --o-filter-stats filter-stats.qza
    # shellcheck disable=SC2086
    qiime deblur denoise-16S --i-demultiplexed-seqs demux-filtered.qza \
      $DEBLUR_ARGS --p-sample-stats --p-jobs-to-start "$THREADS" \
      --o-table table.qza --o-representative-sequences rep-seqs.qza \
      --o-stats denoising-stats.qza
  else
    # shellcheck disable=SC2086
    qiime dada2 denoise-$DENOISE --i-demultiplexed-seqs "$READS" \
      $DADA2_ARGS --p-n-threads "$THREADS" \
      --o-table table.qza --o-representative-sequences rep-seqs.qza \
      --o-denoising-stats denoising-stats.qza
  fi
fi

if [ ! -s taxonomy-nb.qza ]; then
  log "naive Bayes"
  qiime feature-classifier extract-reads --i-sequences "$REF/refseq.qza" \
    --p-f-primer "$PRIMER_F" --p-r-primer "$PRIMER_R" \
    --p-n-jobs "$THREADS" --o-reads ref-trimmed.qza
  qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ref-trimmed.qza \
    --i-reference-taxonomy "$REF/taxonomy.qza" --o-classifier classifier.qza
  qiime feature-classifier classify-sklearn --i-classifier classifier.qza \
    --i-reads rep-seqs.qza --p-n-jobs "$THREADS" \
    --o-classification taxonomy-nb.qza
fi

if [ ! -s taxonomy-vsearch.qza ]; then
  log "vsearch"
  qiime feature-classifier classify-consensus-vsearch --i-query rep-seqs.qza \
    --i-reference-reads "$REF/refseq.qza" \
    --i-reference-taxonomy "$REF/taxonomy.qza" \
    --p-perc-identity 0.98 --p-query-cov 0.9 --p-maxaccepts all \
    --p-top-hits-only --p-min-consensus 0.51 --p-threads "$THREADS" \
    --o-classification taxonomy-vsearch.qza --o-search-results vsearch-hits.qza
fi

for clf in nb vsearch; do
  log "q2-periodontitis ($clf)"
  qiime periodontitis collapse-complexes --i-table table.qza \
    --i-taxonomy "taxonomy-$clf.qza" \
    --o-complex-table "complexes-$clf.qza" --o-indices "indices-$clf.qza" \
    --o-mapping "mapping-$clf.qza" > /dev/null 2> "warnings-$clf.txt" || {
      cat "warnings-$clf.txt"; exit 1; }
  rm -rf "export-$clf"
  for a in complexes indices mapping taxonomy; do
    qiime tools export --input-path "$a-$clf.qza" \
      --output-path "export-$clf/$a" > /dev/null
  done
done
for a in table denoising-stats; do
  [ -s "export-$a" ] || qiime tools export --input-path "$a.qza" \
    --output-path "export-$a" > /dev/null
done
qiime --version | head -1 > versions.txt
qiime info | grep -E "periodontitis|dada2|deblur|feature-classifier" >> versions.txt
log "DONE"
