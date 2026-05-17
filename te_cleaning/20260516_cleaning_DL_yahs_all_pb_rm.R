library(data.table)
setDTthreads(8)

setwd("/data/scripts/slugsTE/te_cleaning")
rm_file <- fread("../data/DL_yahs_all_pb_final.fasta.out", fill=TRUE, skip=3, header=FALSE)
rm_file
setnames(rm_file, c("swscore", "percdiv", "percdel", "percins","qseqid", "qstart", "qend", "qleft", "strand", "repname", "repclass", "sstart",  "send", "sleft", "sid", "starred"))

# cleaning sstart send and sleft columns based on strand
rm_file[strand == "+", `:=`(sstart2 = as.numeric(sstart), send2 = as.numeric(send), sleft2 = as.numeric(stringr::str_extract(sleft, "\\d+")))]
rm_file[strand == "C", `:=`(sstart2 = as.numeric(sleft), send2 = as.numeric(send), sleft2 = as.numeric(stringr::str_extract(sstart, "\\d+")))]
rm_file[,':='(sstart = NULL, send = NULL, sleft = NULL)]
setnames(rm_file, c("sstart2", "send2", "sleft2"), c("sstart", "send", "sleft"))
rm_file

rm_file[,.N,sid][order(-N)]

# smallrna counts 
smallrna <- fread("../data/piRNA_counts_matrix.per_ins.with_coords.tsv")
smallrna[,':='(qseqid=chr, qstart=start+1, qend=end)]
smallrna[,':='(chr=NULL, start=NULL, end=NULL)]

# clean up experiment names
expnames <- names(smallrna)[grepl("27_30", names(smallrna))] 
newnames <- stringr::str_remove(stringr::str_remove(tolower(expnames),".cellrnaclean.27_30nt" ),"dl")
setnames(smallrna, expnames, newnames)
setnames(smallrna, "qseqid", "qseqid")

# scaling per replicate:
smallrna[,lapply(.SD, function(x) x / sum(x)), .SDcols=newnames]
# check correlations but just for those nonzero  
ll <- smallrna[,..newnames]

cors <- cor(ll[rowSums(ll) > 0,])
corrplot::corrplot(cors, method="color" )
corrplot::corrplot.mixed(cors, upper="number", lower="color",  number.cex=0.8, tl.pos="lt")
smallrna[,':='(juvenile = juv1 + juv2 + juv3 )]
smallrna[,':='(ovotestis = ovo1 + ovo2 + ovo3 + jon2ovo1 + jon2ovo2 + jon2ovo3)]
smallrna[,':='(ovotestis_2rep =  ovo1 + ovo2 + ovo3 + jon2ovo2 + jon2ovo3)]
expnames_new <- c("juvenile", "ovotestis", "ovotestis_2rep", newnames)


rm_file_smallrna <- merge(rm_file, smallrna, by.x=c("qseqid", "qstart", "qend", "strand"), by.y=c("qseqid", "qstart", "qend", "strand"), all.x=TRUE)
rm_file_smallrna[is.na(rm_file_smallrna)] <- 0
rm_file_smallrna[,repclass2:=stringr::str_remove(repclass, "\\/.*")]
rm_file_smallrna[,.N,repclass2]

rm_file_smallrna <- rm_file_smallrna[repclass2%in%c("LINE","Unknown","DNA","LTR","SINE","RC","SINE?")]

## now lets first sum the smallrna counts per id:
#rm_file_smallrna[,lapply(.SD, sum), by=sid, .SDcols=expnames_new]

# now we will collapse repeats with same repname, which are overlapping, and will put start as min strart and end as max end
library(GenomicRanges)
rm_file_smallrna_gr <- GRanges(rm_file_smallrna$qseqid, IRanges(rm_file_smallrna$qstart, rm_file_smallrna$qend), 
  strand= ifelse(rm_file_smallrna$strand == "+", "+", "-"),
  sid=rm_file_smallrna$sid,
  repname=rm_file_smallrna$repname, 
  repclass=rm_file_smallrna$repclass,
  juvenile=rm_file_smallrna$juvenile,
  ovotestis=rm_file_smallrna$ovotestis,
  ovotestis_2rep=rm_file_smallrna$ovotestis_2rep
  )

# split by repname so reduce only merges overlapping ranges of the same repeat
grl_orig <- split(rm_file_smallrna_gr, rm_file_smallrna_gr$repname)
grl_red  <- reduce(grl_orig, with.revmap=TRUE)  # revmap indices are local per repname

grl_annotated <- Map(function(red, orig) {
    revmap <- red$revmap
    red$juvenile       <- sum(extractList(orig$juvenile,       revmap))
    red$ovotestis      <- sum(extractList(orig$ovotestis,      revmap))
    red$ovotestis_2rep <- sum(extractList(orig$ovotestis_2rep, revmap))
    red$repname        <- orig$repname[1]   # all same within group
    red$repclass       <- orig$repclass[1]  # all same within group
    red$revmap         <- NULL
    red
}, as.list(grl_red), as.list(grl_orig))

rm_file_smallrna_gr <- unlist(GRangesList(grl_annotated))

ff <- as.data.table(rm_file_smallrna_gr)
fwrite(ff, "../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv", sep="\t", quote=FALSE)


dna <- fread("../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv")
dna[repclass%like%"DNA"]

clean_name <- function(x) substr(gsub("[\\[\\]:*?/\\\\]", "_", x), 1, 31)

rep_df <- as.data.frame(dna)

rep_df_dna <- rep_df[grepl("DNA", rep_df$repclass, ignore.case = TRUE), ]

wb <- createWorkbook()
for (rn in unique(rep_df_dna$repname)) {
    addWorksheet(wb, sheetName = clean_name(rn))
    writeData(wb, sheet = clean_name(rn), rep_df_dna[rep_df_dna$repname == rn, ])
}
saveWorkbook(wb, "rm_file_smallrna_DNA.xlsx", overwrite = TRUE)

