library(data.table)

# readin the original rm_table with orfs and pirna:
rm_file_with_orfs <- fread("../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv", sep="\t", quote=FALSE)
rm_file_with_orfs[,repclass2 := stringr::str_extract(repclass,"^[^/]+")]


### per class, make a folder and save orf table, subfolder with excels for each separate repname

#te <- rm_file_with_orfs[!repclass2 %like% "DNA"]
te <- rm_file_with_orfs
te[repclass2=="SINE?", repclass2:="SINE_"]

te[is.na(orf_length),orf_length:=0]
te <- te[order(repclass2, repclass, repname, -width)]

# i will make a pdf of size distribution of these unknowns: 
library(ggplot2)
  
makeplotforclass <- function(classname){
  print(paste0("Making plot for class: ", classname))

  # subset only to this class:
  tesubset <- te[repclass2==classname]

  te_top <- tesubset[order(-width),.SD[1], by=repname]
  a <- ggplot(te_top, aes(width))+geom_histogram(bins=50, color="black", fill="steelblue") +
    theme_bw() + 
    facet_grid(repclass2~., scales="free_y") +
    xlab("Length of longest insertion per family") + 
    ylab("Count of families") + 
    ggtitle("Distribution of longest insertion per family")
  a2 <- ggplot(te_top[width<10000], aes(width))+
    geom_histogram(bins=100, color="black", fill="steelblue") +
    theme_bw() + 
    facet_grid(repclass2~., scales="free_y") +
    xlab("Length of longest insertion per family") + 
    ylab("Count of families") + 
    ggtitle("Distribution of longest insertion per family")
  a3 <- ggplot(te_top, aes(orf_length))+
    geom_histogram(bins=100, color="black", fill="steelblue") +
    theme_bw() + 
    facet_grid(repclass2~., scales="free_y") +
    xlab("Length of longest orf per family") + 
    ylab("Count of families") +  
    ggtitle("Distribution of longest orf per family")
  a4 <- ggplot(te_top[orf_length>0], aes(orf_length))+
    geom_histogram(bins=100, color="black", fill="steelblue") +
    theme_bw() + 
    facet_grid(repclass2~., scales="free_y") +
    xlab("Length of longest orf per family") + 
    ylab("Count of families") +  
    ggtitle("Distribution of longest orf per family - nonzero families")
  
  aa <- list(a,a2,a3,a4)
  ggsave(paste0("../data/longest_insertion_size_distribution_", classname, ".pdf"), aa, width=10, height=5)  
  print(paste0("Plot saved: ../data/longest_insertion_size_distribution_", classname,".pdf"))
 
  
  # now for each repname in this class, make a histogram of size distribution:
  p <- lapply(unique(tesubset$repname), function(i){
    ggplot(tesubset[repname==i], aes(fill=repclass, x= width)) + 
    geom_histogram(bins=50,color="black") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(i) + xlab("Size distribution of elements") + ylab("Count of elements") 

  })

  ggsave(paste0("../data/size_distribution_repsize_", classname, ".pdf"), p, width=10, height=5)
  print(paste0("Plot saved: ../data/size_distribution_repsize_", classname, ".pdf"))
  
  
  # now for each repname in this class, make a histogram of size distribution:
  p2 <- lapply(unique(te_top[orf_length>0][order(-width)]$repname), function(i){
    ggplot(tesubset[width<20000][repname==i], aes(fill=repclass, x= width)) + 
    geom_histogram(bins=50,color="black") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(i) + xlab("Size distribution of elements") + ylab("Count of elements") 

  })
  

  ggsave(paste0("../data/size_distribution_repsize_nonzeroorf_", classname, ".pdf"),p2, , width=10, height=5)
  print(paste0("Plot saved: ../data/size_distribution_repsize_nonzeroorf_", classname, ".pdf"))

}

lapply(unique(te$repclass2), makeplotforclass)



maketablesforclass <- function(classname){
  print(paste0("Making tables for class: ", classname))
  # subset only to this class:
  tesubset <- te[repclass2==classname]
  
  # two biggest longorf elements per repname:
  te_toporf <- tesubset[order(-orf_length),.SD[1:2], by=repname]
  te_toporf[,anno_element:=paste0(seqnames,":", start, "-", end)]
  te_toporf[,n:=1:.N, by=repname]
  
  ddorf <- dcast(te_toporf, repname~n, value.var=c("anno_element","orf_length"), fill=0)
  
  # and order them by: orf length and targeting in the ovotestis
  te_top <- tesubset[order(-width),.SD[1:5], by=repname]
  te_top[,anno_element:=paste0(seqnames,":", start, "-", end)]
  te_top[,n:=1:.N, by=repname]
  dd <- dcast(te_top[], repname~n, value.var="anno_element", fill=0)
  
  # summary table per family:
  te_summary <- tesubset[,.(ovotestis_total_family=sum(ovotestis),
            juvenile_total_family=sum(juvenile)
  ), by=.(repclass, repname)]

  te_top <- merge(te_summary, dd, by="repname", all.x=TRUE)
  te_top <- merge(te_top, ddorf, by="repname", all.x=TRUE)

  # write it into one excel:

  library(openxlsx)
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = classname)
  writeData(wb, sheet = classname, te_top)
  saveWorkbook(wb, paste0("../data/full_length_tables/", classname, "_top5_per_repname_file_smallrna_orf.xlsx"), overwrite = TRUE)
  # now this table will be manually checked by people.
  # it is uploaded to: 

  clean_name <- function(x) substr(gsub("[\\[\\]:*?/\\\\]", "_", x), 1, 31)

  dir.create(paste0("../data/full_length_tables/full_", classname, "_tables/"), showWarnings = FALSE)
  # create and save individual excel files for each repname:
  for (rn in unique(te_top[orf_length_1>0]$repname)) {
    wb <- createWorkbook()
      addWorksheet(wb, sheetName = clean_name(rn))
      writeData(wb, sheet = clean_name(rn), te[repname == rn, ][order(-orf_length)])
      saveWorkbook(wb, paste0("../data/full_length_tables/full_", classname, "_tables/", clean_name(rn), "_file_smallrna_orf.xlsx"), overwrite = TRUE)

  }
}

te[,.N,repclass2]
lapply(unique(te$repclass2), maketablesforclass)



### extra: making table for all classes, one mastertable, probably would have been smartter aproach :)


{
    te_toporf <- rm_file_with_orfs[order(-orf_length),.SD[1:2], by=.(repclass2,repclass,repname)]
    te_toporf[,anno_element:=paste0(seqnames,":", start, "-", end)]
    te_toporf[,n:=1:.N, by=.(repclass2,repclass,repname)]
  
  ddorf <- dcast(te_toporf, repclass2+repclass+repname~n, value.var=c("anno_element","orf_length"), fill=0)
  
  # and order them by: orf length and targeting in the ovotestis
  te_top <- rm_file_with_orfs[order(-width),.SD[1:5], by=.(repclass2,repclass,repname)]
  te_top[,anno_element:=paste0(seqnames,":", start, "-", end)]
  te_top[,n:=1:.N, by=.(repclass2,repclass,repname)]
  dd <- dcast(te_top[], repclass2+repclass+repname~n, value.var="anno_element", fill=0)
  
  # summary table per family:
  te_summary <- rm_file_with_orfs[,.(ovotestis_total_family=sum(ovotestis),
            juvenile_total_family=sum(juvenile)
  ), by=.(repclass, repname)]

  te_top <- merge(te_summary, dd, by="repname", all.x=TRUE)
  te_top <- merge(te_top, ddorf, by="repname", all.x=TRUE)

 fwrite(te_top, "../data/all_classes_top5_per_repname_file_smallrna_orf.tsv", sep="\t", quote=FALSE)
}
