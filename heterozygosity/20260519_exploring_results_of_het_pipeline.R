library(data.table)

flist <- list.files("../data/heterozygosity/othersamples/", pattern = "*het_per_contig.tsv", full.names = TRUE, recursive = TRUE)
readinhet <- function(file) {
  dt <- fread(file)
  dt[,file:=stringr::str_remove(stringr::str_remove(file,"../data/heterozygosity/othersamples//"),"/het_per_contig.tsv")]  
  return(dt)
}

flist2 <- list.files("../data/heterozygosity/othersamples/", pattern = "*het_per_window.annotated.tsv", full.names = TRUE, recursive = TRUE)
readinhet100 <- function(file) {
  dt <- fread(file)
  dt[,file:=stringr::str_remove(stringr::str_remove(file,"../data/heterozygosity/othersamples//"),"/het_per_contig.tsv")]  
  return(dt)
}

hh <- rbindlist(lapply(flist2, readinhet))
hh[,file:=stringr::str_remove(file,"/het_per_window.annotated.tsv")] 
hh[file=="het_L681/L681", file:="L681"]
mm

het_dt <- rbindlist(lapply(flist, readinhet))
het_dt[file=="het_L681/L681", file:="L681"]


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
httpgd::hgd()

# your query coords you want to lift over
hh
# maps_dt has: query_name, query_start, query_end, target_name, target_start, target_end
# find which PAF alignment each query range falls in, then project to target coords
lifted <- maps_dt[hh, on = .(??_name == contig), allow.cartesian = TRUE][
    start >= query_start & end <= query_end,
    .(target_name,
      target_lifted_start = fifelse(strand == "+",
        target_start + (start - query_start),
        target_end   - (end   - query_start)),
      target_lifted_end = fifelse(strand == "+",
        target_start + (end   - query_start),
        target_end   - (start - query_start))
)
]

lifted

library(ggplot2)

maps_dt[target_name=="scaffold_1" & target_end<10000000 & target_start>5000000]
maps_dt[,.N,file]
maps_dt[,file:=stringr::str_extract(file, "L\\d+")]
maps_dt[,contig:=stringr::str_extract(query_name, "ptg\\d+")]
maps_dt[,contig:=paste0(file,"_", contig)]
het_dt[,contig:=stringr::str_replace(paste0(file,"_", contig),"l$","")]
setkey( maps_dt, file, contig)
setkey( het_dt, file, contig)
het_dt <- unique(het_dt)
mm <- merge(maps_dt, het_dt)

# adding frankfurt
mm <- rbind(mm, hh[file=="DerLae",.(file, target_name=contig, target_start=start, target_end=end, het_per_100kb)], fill=TRUE)


mm[order(target_start)]
mm[,nm:=stringr::str_remove(stringr::str_remove(query_name, "L\\d+_"),"_ptg\\d+l")]
mm[is.na(nm), nm := "DerLae_Frankfurt"]
ggplot(mm[target_name=="scaffold_6"], aes(x=target_start, y=het_per_100kb, color=file)) + 
geom_segment(aes(x=target_start, xend=target_end, y=het_per_100kb, yend=het_per_100kb), size=2) + 
facet_wrap(target_name~nm)+
theme_bw()+ coord_cartesian(ylim=c(-10, 50))
httpgd::hgd()

# for plot limit to 50:
mm[het_per_100kb>50, het_per_100kb:=50]
mm[order(-target_start),.SD[1],target_name][,target_name]
smmal <- mm[target_name%in%paste0("scaffold_",1:30)]
lvls <- smmal[order(-target_length),.SD[1],target_name][,target_name]
smmal[,target_name:=factor(as.character(target_name), levels=lvls)]
smmal[,nm:=factor(stringr::str_extract(target_name, "\\d+"), levels=c(1:30))]
ggplot(smmal, aes(x=target_start, y=het_per_100kb, color=file)) + 
geom_segment(aes(x=target_start, xend=target_end, y=het_per_100kb, yend=het_per_100kb), size=2) + 
facet_grid(nm~target_name, scales="free_x", space="free")+
theme_bw()+ theme(
    panel.spacing.x = unit(0, "lines")  )+
 coord_cartesian(ylim=c(0, 52))


ggplot(smmal, aes(x=target_start, y=het_per_100kb, color=file)) + 
geom_point() + 
facet_grid(factor(as.character(nm),levels=as.character(1:30))~target_name, scales="free_x", space="free")+
theme_bw()+ theme(
    panel.spacing.x = unit(0, "lines")  )+
 coord_cartesian(ylim=c(0, 52))

