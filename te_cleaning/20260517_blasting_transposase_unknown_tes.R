library(data.table)

# tBLASTn: protein (KAK7092566.1) vs nucleotide unknowns
tblastn_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                  "qlen", "slen", "qcovs", "sstrand")
tbn <- fread("../data/blasting_transposase_unknown_tes/tblastn_KAK7092566_vs_unknowns.tsv")
setnames(tbn, tblastn_cols)
tbn[, c("repclass", "repname", "coords") := tstrsplit(sseqid, "\\|")]

# BLASTx: nucleotide unknowns vs protein (KAK7092566.1) — best hit per query, all 6 frames
blastx_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                 "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                 "qlen", "slen", "qcovs")
bx <- fread("../data/blasting_transposase_unknown_tes/blastx_unknowns_vs_KAK7092566.tsv")
setnames(bx, blastx_cols)
bx[, c("repclass", "repname", "coords") := tstrsplit(qseqid, "\\|")]

# frame is +1/+2/+3 (forward) or -1/-2/-3 (reverse)
bx[, .N, frame]
summary_results <- bx[, .N, repname]
fwrite(summary_results, "../data/blasting_transposase_unknown_tes/blastx_summary.csv")
