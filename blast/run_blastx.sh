#!/bin/bash
set -euo pipefail

DIR="/data/scripts/slugsTE/data/blasting_transposase_unknown_tes"
PROT="$DIR/KAK7092566.1.fasta"
NUCL="$DIR/20260517_DL_yahs_all_pb_final_unknowns_only.fasta"
NUCL_DB="$DIR/unknowns_db"
PROT_DB="$DIR/transposase_prot_db"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"

# --- tBLASTn: protein query vs nucleotide database ---
makeblastdb \
    -in "$NUCL" \
    -dbtype nucl \
    -out "$NUCL_DB" \
    -parse_seqids

tblastn \
    -query "$PROT" \
    -db "$NUCL_DB" \
    -evalue 1e-5 \
    -word_size 3 \
    -matrix BLOSUM62 \
    -seg no \
    -num_threads 8 \
    -outfmt "$FMT sstrand" \
    -out "$DIR/tblastn_KAK7092566_vs_unknowns.tsv"

echo "tBLASTn done."

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
    -out "$DIR/blastx_unknowns_vs_KAK7092566.tsv"

echo "BLASTx done. Results in $DIR/blastx_unknowns_vs_KAK7092566.tsv"
