library(data.table)
library(ggplot2)

# here are the heterozygous sites in the L681 genome:
hets <- fread("../data/het_L681/L681/L681.het.vcf.gz", skip="#CHROM")


# here is the mapping of the L681 genome to the reference genome:
mapping <- fread("../data/het_L681/L681/L681_to_ref_mapping.tsv")


flist <- list.files("../data/heterozygosity", pattern = "*paf", full.names = TRUE)
readinpaf <- function(file) {
  dt <- fread(file, header=FALSE)
  setnames(dt, c("query_name", "query_length", "query_start", "query_end", "strand", "target_name", "target_length", "target_start", "target_end", "num_matching_bases", "alignment_block_length", "mapping_quality"))
  dt[,file:=basename(file)]
  return(dt)
}
paf_dt <- rbindlist(lapply(flist, readinpaf))

library(GenomicRanges)
#we now make disjunkt union of all organisms;
# and annotate which are covered in which genome:
gr <- makeGRangesFromDataFrame(paf_dt, seqnames.field="target_name", start.field="target_start", end.field="target_end", keep.extra.columns=TRUE)
# we want a split granges object per species:by
gr_list <- split(gr, gr$file)

# disjoint union of all granges:
dj <- disjoin(gr)

co <- countOverlaps(dj, gr_list)

cm <- sapply(names(gr_list), function(i){
  countOverlaps(dj,gr_list[[i]])
})
as.data.table(cm)

library(pheatmap)
library(corrplot)

colnames(cm) <- stringr::str_remove(colnames(cm), "Deroceras_laeve_CZ_1_0__vs__")
colnames(cm) <- stringr::str_remove(colnames(cm), ".paf")
corrplot(cor(cm), method="color", addCoef.col = "black", number.cex=0.7, tl.cex=0.7)

dt <- as.data.table(cm)
dt[,disjoin_region:=1:.N]
mm <- melt(dt, id.vars="disjoin_region", variable.name="species", value.name="overlaps")
mm[, (c("name","location","phallus")):=tstrsplit(species, "_")]

frac <- mm[,sum(overlaps), by=.(disjoin_region,phallus)]
dd <- dcast(frac, disjoin_region~phallus, value.var="V1", fill=0)
dd[order(Aphallic-Euphallic), ]

dt[disjoin_region==386064]

dt[L451_Glasgow_Aphallic==0&L633_Wein6_Aphallic==0&L684_Jonathan_Aphallic==0]
dt[,aphalic_covered:=((L451_Glasgow_Aphallic>0)+(L633_Wein6_Aphallic>0)+(L684_Jonathan_Aphallic>0))]
dt[,euphallic_covered:=((L532_Jonathan2_Euphallic>0)+(L685_Jonathan2_Euphallic>0)+(L681_Schlossteich2_Euphallic>0)+(derLae1_Mexico_Euphallic>0))]
dt[order(euphallic_covered-aphalic_covered), .(disjoin_region, aphalic_covered, euphallic_covered)]


regions <- dt[(euphallic_covered==4&aphalic_covered==0)|(euphallic_covered==0&aphalic_covered==3), disjoin_region]
regions_pl <- dt[(euphallic_covered==4&aphalic_covered==0)|(euphallic_covered==0&aphalic_covered==3), ]




dj$coverage_aphalic <- dt[order(disjoin_region),aphalic_covered]
dj$coverage_euphallic <- dt[order(disjoin_region),euphallic_covered]
cov <- as.data.table(dj)
ggplot(cov[regions][seqnames=="scaffold_1"], aes(x=start, y=coverage_aphalic-coverage_euphallic)) + 
  geom_point() + 
  theme_bw() + 
  xlab("Genomic position") +
   ylab("Covered in aphalic - covered in euphallic")+
   facet_wrap(~seqnames, scales="free_x") 

rr <- cov[coverage_aphalic==3&coverage_euphallic==0, .(seqnames, start, end)]
regions_specific_for_aphalic <- as.data.table(reduce(GRanges(rr)))
fwrite(regions_specific_for_aphalic, "../data/heterozygosity/regions_specific_for_aphalic.bed", sep="\t", col.names=FALSE)
regions_specific_for_aphalic[,sum(width)]
