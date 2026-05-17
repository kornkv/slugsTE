# slugsTE

Analysis of transposable elements (TEs) in *Deroceras laeve* (CZ assembly), integrating RepeatMasker annotation with piRNA targeting data and ORF prediction to characterise transposable elements.

## Data inputs

| File | Description |
|---|---|
| `DL_yahs_all_pb_final.fasta.out` | RepeatMasker output for the *D. laeve* assembly |
| `piRNA_counts_matrix.per_ins.with_coords.tsv` | piRNA counts per TE insertion (juvenile + ovotestis replicates) |
| `Deroceras_laeve_CZ.1.0.fasta` | Genome assembly (symlink) |
| `data/blasting_transposase_unknown_tes/KAK7092566.1.fasta` | Transposase protein query (hAT superfamily) |

## Workflow

```
RepeatMasker output
       │
       ▼
20260516_cleaning_DL_yahs_all_pb_rm.R        ← merge piRNA counts, collapse overlapping insertions
       │
       ├──► DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv
       │
       ├──► 20260516_size_distribution_repsize.R          ← size distribution plots per DNA TE family
       │
       └──► 20260517_rm_orf_finding.R                     ← longest ORF per element (all classes)
                   │
                   └──► ..._with_orfs.tsv
                               │
                   ┌───────────┴──────────────────────────┐
                   ▼                                        ▼
       data/blasting_transposase_unknown_tes/      20260517_table_clean_preparation.R
       run_blastx.sh  (BLASTx)           (DNA transposon summary tables)
                   │
                   ▼
       20260517_blasting_transposase_unknown_tes.R  ← parse BLAST results
                   │
                   ▼
       20260517_3_table_clean_preparation_unknowns_dna.R   ← unknown elements with transposase hits (DNA transposon from unknown repeats summary tables)
       20260517_3_table_clean_preparation_other_classes.R  ← other TE classes (in progress)
```

## Scripts

### `te_cleaning/20260516_cleaning_DL_yahs_all_pb_rm.R`
Parses the raw RepeatMasker `.out` file, fixes strand-dependent coordinate columns, and merges piRNA count data (juvenile + ovotestis replicates) per insertion. Collapses overlapping insertions of the same repeat family using `GenomicRanges::reduce()`, summing piRNA counts across merged intervals. Filters to major TE classes (LINE, LTR, DNA, SINE, RC, Unknown). Writes the annotated table to `DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv` and an Excel workbook with one sheet per DNA TE family.

**Output:** `data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv`, `data/rm_file_smallrna_DNA.xlsx`

---

### `te_cleaning/20260516_size_distribution_repsize.R`
Plots size (width) distributions for all DNA TE families as histograms, one panel per family, saved to a single multi-page PDF. 

**Output:** `data/size_distribution_repsize.pdf`

---

### `te_cleaning/20260517_rm_orf_finding.R`
Extracts genomic sequences for all annotated TE insertions using `Biostrings::getSeq()`, then finds the longest ORF per element with `ORFik::findORFs()` (ATG start, TAA|TAG|TGA stop, minimum 50 nt). Merges ORF lengths back into the main annotation table. For DNA transposons specifically, generates a summary Excel with the top 5 longest elements and top 2 ORF-containing elements per family, ordered by ovotestis piRNA targeting.

**Output:** `data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv`, `data/dna_transposons_top5_per_repname_file_smallrna_orf.xlsx`, `data/full_dna_tables/<repname>_file_smallrna_orf.xlsx`

---

### `data/blasting_transposase_unknown_tes/run_blastx.sh`
Builds BLAST databases and runs searches of the nucleotide sequences from the transposable elements against transposase protein `KAK7092566.1`: (we are fishing out elements which have repclass "Unknown" , but are potentially DNA transposons). Later we are adding those to the table with other DNA transposons. 
- **BLASTx** — unknown-TE nucleotide queries vs protein database (all 6 frames, best hit per query)

**Output:**  `blastx_unknowns_vs_KAK7092566.tsv`

---

### `te_cleaning/20260517_blasting_transposase_unknown_tes.R`
Reads tBLASTn and BLASTx output tables, parses sequence IDs into `repclass`, `repname`, and genomic coordinates, and summarises which unknown TE families have transposase hits.

**Output:** `data/blasting_transposase_unknown_tes/blastx_summary.csv`

---

### `te_cleaning/20260517_3_table_clean_preparation_unknowns_dna.R`
Filters the ORF-annotated table to unknown elements that had transposase BLAST hits. Produces size distribution plots, per-family Excel files, and a combined summary table with the top 5 largest and top 2 ORF-containing elements per family ranked by ovotestis piRNA targeting. Also writes a `repname`/`repclass` lookup CSV covering all DNA and transposase-hit unknown elements.

**Output:** `data/size_distribution_repsize_unknowns_dna_transposase.pdf`, `data/unknown_dna_transposons_top5_per_repname_file_smallrna_orf.xlsx`, `data/full_dna_tables/unknown_<repname>_file_smallrna_orf.xlsx`, `data/repname_repclass_dna_unknowns.csv`

---

### `te_cleaning/20260517_table_clean_preparation.R`
General table preparation for DNA transposons — in progress.

### `te_cleaning/20260517_3_table_clean_preparation_other_classes.R`
Table preparation for non-DNA TE classes — in progress.
