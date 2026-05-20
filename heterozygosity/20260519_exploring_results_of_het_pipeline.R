library(data.table)

flist <- list.files("../data/heterozygosity/othersamples/", pattern = "*het_per_contig.tsv", full.names = TRUE, recursive = TRUE)
readinhet <- function(file) {
  dt <- fread(file)
  dt[,file:=stringr::str_remove(stringr::str_remove(file,"../data/heterozygosity/othersamples//"),"/het_per_contig.tsv")]  
  return(dt)
}

het_dt <- rbindlist(lapply(flist, readinhet))
het_dt


mapping_to_ref <- list.files("../data/heterozygosity/", pattern = "*paf", full.names = TRUE, recursive = TRUE)
mapping_to_ref

maps <- lapply(mapping_to_ref, function(file) {
  dt <- fread(file, header=FALSE)
  setnames(dt, c("query_name", "query_length", "query_start", "query_end", "strand", "target_name", "target_length", "target_start", "target_end", "num_matching_bases", "alignment_block_length", "mapping_quality"))
  dt[,file:=stringr::str_remove(stringr::str_remove(file,"../data/heterozygosity//Deroceras_laeve_CZ_1_0__vs__"),".paf")]
  return(dt)
})
maps_dt <- rbindlist(maps)
maps_dt


top30 <- het_dt[order(-contig_length),.SD[1:30], by=file]
library(ggplot2)
httpgd::hgd()
ggplot(top30, aes(x=contig_length, y=het_per_100kb, color=file)) + geom_point() + theme_bw() + scale_x_log10()

maps_dt[target_name=="scaffold_1" & target_end<10000000 & target_start>5000000]
