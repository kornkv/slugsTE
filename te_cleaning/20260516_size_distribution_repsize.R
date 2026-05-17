library(data.table)
x <- fread("../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv")
dnas <- x[repclass %like% "DNA"]
library(ggplot2)

ggplot(dnas, aes(fill=repclass, x= width)) + 
  geom_histogram(bins=50,color="black") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("DLcz_rnd-1_family-136") + xlab("Size distribution of reads mapping to DLcz_rnd-1_family-136") + ylab("Count of reads") 
  

pdf("../data/size_distribution_repsize.pdf", width=10, height=5)

lapply(unique(dnas$repname), function(i){
  ggplot(dnas[repname==i], aes(fill=repclass, x= width)) + 
  geom_histogram(bins=50,color="black") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(i) + xlab("Size distribution of reads") + ylab("Count of reads") 

})

dev.off()
httpgd::hgd()
