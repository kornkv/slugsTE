#!/bin/bash
set -euo pipefail

DIR="/data/scripts/slugsTE/data/blasting_transposase_unknown_tes"
PROT="$DIR/maverick_orf_protein.fasta"
NUCL="$DIR/20260517_DL_yahs_all_pb_final_unknowns_only.fasta"
PROT_DB="$DIR/maverick_orf_protein_db"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"

# --- BLASTx: nucleotide queries vs protein database (all 6 frames, best hit per query) ---
makeblastdb \
    -in "$PROT" \
    -dbtype prot \
    -out "$PROT_DB" \
    -parse_seqids

blastx \
    -query "$NUCL" \
    -db "$PROT_DB" \
    -evalue 1e-5 \
    -word_size 3 \
    -matrix BLOSUM62 \
    -seg no \
    -num_threads 1 \
    -max_target_seqs 5 \
    -max_hsps 1 \
    -outfmt "$FMT frame" \
    -out "$DIR/blastx_unknowns_vs_maverick_orf.tsv"

echo "BLASTx done. Results in $DIR/blastx_unknowns_vs_maverick_orf.tsv"
