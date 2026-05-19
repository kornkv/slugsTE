library(data.table)
library(ggplot2)


pirna_clusters <- fread("../data/explore_pirna_cluster/piRNA_clusters_DL_CZ_fixed.sorted.bed")

setnames(pirna_clusters, c("seqnames", "start", "end", "cluster_id", "V5","V6"))
setkey(pirna_clusters, "seqnames", "start", "end")
pirna_clusters[,width:=end-start]
# check if i have data for pirna expression in the pirna clusters:
smallrna <- fread("../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv")
ovls <- foverlaps(smallrna, pirna_clusters, by.x=c("seqnames", "start", "end"), by.y=c("seqnames", "start", "end"), type="within", nomatch=0L)[order(seqnames, start)]

# ovls represent only the transposable elements within pirna clusters.
ovls
ovls[,repclass2:=stringr::str_remove(repclass, "\\/.*")]
ovls[repclass2=="SINE?",repclass2:="SINE_"]
ovls[ovotestis>0][,.N,.(repclass2, repname)][order(-N)]
# now lets visualize cluster of highest total expression :
ovls[,total_exp_repname:=sum(ovotestis), by=.(cluster_id, repclass2, repclass, repname)]
ovls[,elementlocation:=paste0(seqnames, ":", i.start, "-", i.end)]
ovls[repclass2=="DNA"][order(-total_exp_repname)][order(-orf_length)][1:20,]

#OK calculate per pirnacluster total pirna expression and get top 20
top20pirna <- ovls[,.(total_exp_cluster=sum(ovotestis)), by=.(cluster_id, seqnames, start, end, width)][order(-total_exp_cluster)][1:20]

library(GenomicRanges)
library(Biostrings)
library(BSgenome)
gr <- GRanges(top20pirna$seqnames, IRanges(top20pirna$start, top20pirna$end), cluster_id=top20pirna$cluster_id, total_exp_cluster=top20pirna$total_exp_cluster)

dna_seq <- readDNAStringSet("../data/Deroceras_laeve_CZ.1.0.fasta")

# Extract sequence at those coordinates
seq_at_coords <- getSeq(dna_seq, gr)
names(seq_at_coords) <- paste0(gr$cluster_id,"_",seqnames(gr),"_",start(gr),"_",end(gr),"_",gr$total_exp_cluster)
writeXStringSet(seq_at_coords, "../data/explore_pirna_cluster/top20_pirna_clusters_sequences.fasta")
