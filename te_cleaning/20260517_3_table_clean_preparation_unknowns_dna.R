library(data.table)

# read in the fished out annotation of unknowns with transposase:
interesting_unknowns <- fread("../data/blasting_transposase_unknown_tes/blastx_summary.csv")

# readin the original rm_table with orfs and pirna:
rm_file_with_orfs <- fread("../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv", sep="\t", quote=FALSE)

#now extract only the DNA unknowns:  
dna_only_unknowns <- rm_file_with_orfs[repname %in% interesting_unknowns$repname]
dna_only <- dna_only_unknowns[order(repname, -width)]
dna_only[is.na(orf_length),orf_length:=0]

# i will make a pdf of size distribution of these unknowns: 
library(ggplot2)
  

pdf("../data/size_distribution_repsize_unknowns_dna_transposase.pdf", width=10, height=5)

lapply(unique(dna_only$repname), function(i){
  ggplot(dna_only[repname==i], aes(fill=repclass, x= width)) + 
  geom_histogram(bins=50,color="black") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle(i) + xlab("Size distribution of reads") + ylab("Count of reads") 

})

dev.off()

dna_summary <- dna_only[,.(ovotestis_total_family=sum(ovotestis),
            juvenile_total_family=sum(juvenile)
), by=repname]
perfamily <- dna_summary[order(-ovotestis_total_family)]


# two biggest longorf elements per repname:
dna_toporf <- dna_only[order(-orf_length),.SD[1:2], by=repname]
dna_toporf[,anno_element:=paste0(seqnames,":", start, "-", end)]
dna_toporf[,n:=1:.N, by=repname]

ddorf <- dcast(dna_toporf, repname~n, value.var=c("anno_element","orf_length"), fill=0)

# and order them by: orf length and targeting in the ovotestis
dna_top <- dna_only[order(-width),.SD[1:5], by=repname]
dna_top[,anno_element:=paste0(seqnames,":", start, "-", end)]
dna_top[,n:=1:.N, by=repname]
dd <- dcast(dna_top[], repname~n, value.var="anno_element", fill=0)
dna_top <- merge(perfamily, dd, by="repname", all.x=TRUE)

dna_top <- merge(dna_top, ddorf, by="repname", all.x=TRUE)



# write it into one excel:

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, sheetName = "top5_per_repname")
writeData(wb, sheet = "top5_per_repname", dna_top)
saveWorkbook(wb, "../data/unknown_dna_transposons_top5_per_repname_file_smallrna_orf.xlsx", overwrite = TRUE)
# now this table will be manually checked by people.
# it is uploaded to: https://kuzmanconsultingdoo-my.sharepoint.com/:x:/g/personal/maja_kuzmanconsulting_com/IQDXgpcR-iMITqMNjfvyMXMfATOKAe6BoHh51AZFW75jd5o?e=PynPV4 


clean_name <- function(x) substr(gsub("[\\[\\]:*?/\\\\]", "_", x), 1, 31)

# create and save individual excel files for each repname:
for (rn in unique(dna_only$repname)) {
  wb <- createWorkbook()
    addWorksheet(wb, sheetName = clean_name(rn))
    writeData(wb, sheet = clean_name(rn), dna_only[repname == rn, ][order(-width)])
    saveWorkbook(wb, paste0("../data/full_dna_tables/unknown_", clean_name(rn), "_file_smallrna_orf.xlsx"), overwrite = TRUE)

}


dna_only[is.na(orf_length)]

# now getting the repclass for the elements to put into the table:
fwrite(rm_file_with_orfs[(repclass %like% "DNA") | (repname %in% interesting_unknowns$repname),.N, .(repname, repclass)],
file="../data/repname_repclass_dna_unknowns.csv")

# the unknowns and the annotation are added to the table on onedrive. 
# moving on to /data/scripts/slugsTE/te_cleaning/20260517_3_table_clean_preparation_other_classes.R