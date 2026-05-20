# Install
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ORFik", "Biostrings", "GenomicRanges", "BSgenome"))

library(ORFik)
library(Biostrings)
library(GenomicRanges)
library(data.table)

# Load assembled genome (FASTA)
genome <- readDNAStringSet("../data/Deroceras_laeve_CZ.1.0.fasta")

# Define your coordinates as a GRanges object
rm_file <- fread("../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv")
dnas <- rm_file[repclass %like% "DNA"]


coords <- GRanges(
  seqnames = dnas$seqnames,
  ranges   = IRanges(start = dnas$start, end = dnas$end)
#  strand   = dnas$strand
)
coords_all <- GRanges(
  seqnames = rm_file$seqnames,
  ranges   = IRanges(start = rm_file$start, end = rm_file$end)
#  strand   = rm_file$strand # we want all orientation orfs
)


# Extract sequence at those coordinates
seq_at_coords <- getSeq(genome, coords)
seq_at_coords_all <- getSeq(genome, coords_all)

# write fasta for sequence at coordinates:
rm_file[,namesforfasta:=paste0(repclass,"|",repname,"|",seqnames, ":", start, "-", end)]
names(seq_at_coords_all) <- rm_file$namesforfasta
writeXStringSet(seq_at_coords_all, "../data/20260517_DL_yahs_all_pb_final.fasta")

unknowns <- seq_at_coords_all[grep("Unknown", names(seq_at_coords_all))]
writeXStringSet(unknowns, "../data/20260517_DL_yahs_all_pb_final_unknowns_only.fasta")

# Find all ORFs in the extracted sequence
orfs <- findORFs(
  seq_at_coords_all,
  startCodon   = "ATG",
  stopCodon    = "TAA|TAG|TGA",
  longestORF   = TRUE,   # TRUE = only longest per region
  minimumLength = 50     # min ORF length in nt
)

orfs_multiple <- findORFs(
  seq_at_coords_all,
  startCodon   = "ATG",
  stopCodon    = "TAA|TAG|TGA",
  longestORF   = FALSE,   # FALSE = multiple ORFs per region
  minimumLength = 50     # min ORF length in nt
)

seq_at_coords_all_rc <- reverseComplement(seq_at_coords_all)
names(seq_at_coords_all_rc) <- names(seq_at_coords_all)
orfs_multiple_neg <- findORFs(
  seq_at_coords_all_rc,
  startCodon   = "ATG",
  stopCodon    = "TAA|TAG|TGA",
  longestORF   = FALSE,   # FALSE = multiple ORFs per region
  minimumLength = 50     # min ORF length in nt
)

orfs_multiple <- as.data.table(unlist(orfs_multiple))  # ORFs are returned as a GRangesList, unlist to get a single GRanges
orfs_multiple_neg <- as.data.table(unlist(orfs_multiple_neg))  # ORFs are returned as a GRangesList, unlist to get a single GRanges

orfs <- rbind(orfs_multiple, orfs_multiple_neg)

orfs <- orfs[order(names)]
orfs[order(-width)][width<5000,.N,names]


# Keep only the longest two ORFs per input sequence
longest_orfs <- orfs[order(-width),.SD[1:2], by=names]
setnames(longest_orfs,c("start","end","width"),c("orf_start","orf_end","orf_length")) 
longest_orfs[,n:=1:.N, by=names]
lo <- dcast(longest_orfs, names~n, value.var=c("orf_start","orf_end","orf_length"), fill=0)
lo[order(-orf_length_1)]
# names of the orfs are position of elements in the original list

rm_file[,names:=as.character(1:.N)]
rm_file_with_orfs <- merge(rm_file, lo, by="names", all.x=TRUE)
setnames(rm_file_with_orfs, c( "orf_start_1", "orf_end_1", "orf_length_1"), c("orf_start", "orf_end", "orf_length"))
rm_file_with_orfs[repclass %like% "DNA" ][order(-orf_length)]


fwrite(rm_file_with_orfs, "../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_two_orfs.tsv", sep="\t", quote=FALSE)
# now go to 20260517_4_table_clean_preparation_other_classes.R to make tables for all classes for two orfs



# old code uses ../data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv just longest orf
#now extract only the DNA 
dna_only <- rm_file_with_orfs[repclass %like% "DNA"]
dna_only <- dna_only[order(repname, -width)]
dna_only[is.na(orf_length),orf_length:=0]


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
saveWorkbook(wb, "../data/dna_transposons_top5_per_repname_file_smallrna_orf.xlsx", overwrite = TRUE)
# now this table will be manually checked by people.
# it is uploaded to: https://kuzmanconsultingdoo-my.sharepoint.com/:x:/g/personal/maja_kuzmanconsulting_com/IQDXgpcR-iMITqMNjfvyMXMfATOKAe6BoHh51AZFW75jd5o?e=PynPV4 


clean_name <- function(x) substr(gsub("[\\[\\]:*?/\\\\]", "_", x), 1, 31)

# create and save individual excel files for each repname:
for (rn in unique(dna_only$repname)) {
  wb <- createWorkbook()
    addWorksheet(wb, sheetName = clean_name(rn))
    writeData(wb, sheet = clean_name(rn), dna_only[repname == rn, ][order(-width)])
    saveWorkbook(wb, paste0("../data/full_dna_tables/", clean_name(rn), "_file_smallrna_orf.xlsx"), overwrite = TRUE)

}


dna_only[is.na(orf_length)]
