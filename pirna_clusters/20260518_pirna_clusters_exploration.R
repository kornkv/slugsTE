library(data.table)
library(ggplot2)


pirna_clusters <- fread("../data/explore_pirna_cluster/piRNA_clusters_DL_CZ_fixed.sorted.bed")
setnames(pirna_clusters, c("seqnames", "start", "end", "cluster_id", "V5","V6"))
setkey(pirna_clusters, "seqnames", "start", "end")
bedgraph <- fread("../data/explore_pirna_cluster/Dlaeve_L532_Jonathan2_Euphallic.bedGraph")
setnames(bedgraph, c("seqnames", "start", "end", "score"))
setkey(bedgraph, "seqnames", "start", "end")

# annotate the bedgraph with the cluster information:
bedgraph_annotated <- foverlaps(bedgraph, pirna_clusters, by.x=c("seqnames", "start", "end"), by.y=c("seqnames", "start", "end"), type="within", nomatch=0L)

# bedgraph is 0-based half-open: 1-based positions = (start+1):end
dt_perbase <- bedgraph[, .(pos = (start + 1):end, coverage = score), by = .I]
dt_perbase[,quantile(coverage, probs=c(0.1,0.9))]

ggplot(dt_perbase[I%in%1:3], aes(x=pos, y=coverage)) +
  geom_point() + theme_bw() + 
  facet_wrap(~I, scales="free_x") + xlab("Position in cluster") + ylab("Coverage")
