#!/bin/bash
set -euo pipefail

QUERY="../data/explore_pirna_cluster/top20_pirna_clusters_sequences.fasta"
GENOME_DIR="/data/genomes"
OUT_DIR="../data/explore_pirna_cluster/blastn_vs_genomes"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"

mkdir -p "$OUT_DIR"

for GENOME in "$GENOME_DIR"/*.fasta; do
    BASENAME=$(basename "$GENOME" .fasta)
    DB="$OUT_DIR/${BASENAME}_db"
    OUT="$OUT_DIR/blastn_top20pirnas_vs_${BASENAME}.tsv"

    echo "Building BLAST db for $BASENAME ..."
    makeblastdb \
        -in "$GENOME" \
        -dbtype nucl \
        -out "$DB" \
        -parse_seqids

    echo "Running BLASTn against $BASENAME ..."
    blastn \
        -query "$QUERY" \
        -db "$DB" \
        -evalue 1e-5 \
        -perc_identity 60 \
        -max_target_seqs 30 \
        -num_threads 8 \
        -outfmt "$FMT" \
        -out "$OUT"

    echo "Done: $OUT"
done

echo "All genomes processed."
