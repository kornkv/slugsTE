# slugsTE

Analysis of transposable elements (TEs) in *Deroceras laeve* (CZ assembly), integrating RepeatMasker annotation with piRNA targeting data and ORF prediction to characterise transposable elements.

## Data inputs

| File | Description |
|---|---|
| `DL_yahs_all_pb_final.fasta.out` | RepeatMasker output for the *D. laeve* assembly |
| `piRNA_counts_matrix.per_ins.with_coords.tsv` | piRNA counts per TE insertion (juvenile + ovotestis replicates) |
| `Deroceras_laeve_CZ.1.0.fasta` | Genome assembly (symlink) |
| `data/blasting_transposase_unknown_tes/KAK7092566.1.fasta` | Transposase protein query (hAT superfamily) |
| `data/blasting_transposase_unknown_tes/maverick_orf_protein.fasta` | Maverick transposase protein query |
| `data/explore_pirna_cluster/piRNA_clusters_DL_CZ_fixed.sorted.bed` | piRNA cluster coordinates (267 clusters) |
| `data/explore_pirna_cluster/top20_pirna_clusters_sequences.fasta` | Sequences of the top 20 piRNA clusters |
| `/data/genomes/*.fasta` | *D. laeve* genome assemblies (11 individuals + reference) |

## Workflow

```
RepeatMasker output
       │
       ▼
te_cleaning/20260516_cleaning_DL_yahs_all_pb_rm.R   ← merge piRNA counts, collapse overlapping insertions
       │
       ├──► DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv
       │
       ├──► te_cleaning/20260516_size_distribution_repsize.R     ← size distribution plots per DNA TE family
       │
       └──► te_cleaning/20260517_rm_orf_finding.R                ← longest ORF per element (all classes)
                   │
                   └──► ..._with_orfs.tsv
                               │
             ┌─────────────────┼──────────────────────────────┐
             ▼                 ▼                               ▼
      blast/run_blastx.sh   blast/run_blastx_maverick.sh    te_cleaning/20260517_table_clean_preparation.R
      (KAK7092566.1)        (Maverick ORF protein)          (DNA transposon summary tables, in progress)
             │                 │
             └────────┬────────┘
                      ▼
      te_cleaning/20260517_blasting_transposase_unknown_tes.R  ← parse BLAST results
                      │
                      ▼
      te_cleaning/20260517_3_table_clean_preparation_unknowns_dna.R   ← unknown elements with transposase hits
      te_cleaning/20260517_4_table_clean_preparation_other_classes.R  ← all non-DNA TE classes

piRNA cluster analysis (independent):
      data/explore_pirna_cluster/bw_extr.sh                ← extract bedGraphs per cluster per sample
      pirna_clusters/20260518_pirna_clusters_exploration.R  ← per-base coverage and cluster annotation
             │
             ▼
      pirna_clusters/blastn_pirnas_vs_genomes.nf            ← BLASTn top 20 piRNA clusters vs all genomes (Nextflow)
             │
             ▼
      data/explore_pirna_cluster/blastn_vs_genomes/*.tsv    ← one tabular result file per genome
```

## Scripts

### `te_cleaning/20260516_cleaning_DL_yahs_all_pb_rm.R`
Parses the raw RepeatMasker `.out` file, fixes strand-dependent coordinate columns, and merges piRNA count data (juvenile + ovotestis replicates) per insertion. Collapses overlapping insertions of the same repeat family using `GenomicRanges::reduce()`, summing piRNA counts across merged intervals. Filters to major TE classes (LINE, LTR, DNA, SINE, RC, Unknown). Writes the annotated table and an Excel workbook with one sheet per DNA TE family.

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

### `blast/run_blastx.sh`
Builds BLAST databases and runs two searches against the hAT transposase protein `KAK7092566.1` to fish out unknown-class elements that are likely DNA transposons:
- **tBLASTn** — protein query vs nucleotide unknown-TE database
- **BLASTx** — unknown-TE nucleotide queries vs protein database (all 6 frames, best hit per query)

**Output:** `data/blasting_transposase_unknown_tes/tblastn_KAK7092566_vs_unknowns.tsv`, `data/blasting_transposase_unknown_tes/blastx_unknowns_vs_KAK7092566.tsv`

---

### `blast/run_blastx_maverick.sh`
BLASTx search of unknown-TE nucleotide sequences against a Maverick transposase ORF protein database to identify Maverick/Polinton-class elements among the unknowns.

**Output:** `data/blasting_transposase_unknown_tes/blastx_unknowns_vs_maverick_orf.tsv`

---

### `te_cleaning/20260517_blasting_transposase_unknown_tes.R`
Reads tBLASTn and BLASTx output tables, parses sequence IDs into `repclass`, `repname`, and genomic coordinates, and summarises which unknown TE families have transposase hits.

**Output:** `data/blasting_transposase_unknown_tes/blastx_summary.csv`

---

### `te_cleaning/20260517_3_table_clean_preparation_unknowns_dna.R`
Filters the ORF-annotated table to unknown elements that had transposase BLAST hits. Produces size distribution plots, per-family Excel files, and a combined summary table with the top 5 largest and top 2 ORF-containing elements per family ranked by ovotestis piRNA targeting. Also writes a `repname`/`repclass` lookup CSV.

**Output:** `data/size_distribution_repsize_unknowns_dna_transposase.pdf`, `data/unknown_dna_transposons_top5_per_repname_file_smallrna_orf.xlsx`, `data/full_dna_tables/unknown_<repname>_file_smallrna_orf.xlsx`, `data/repname_repclass_dna_unknowns.csv`

---

### `te_cleaning/20260517_4_table_clean_preparation_other_classes.R`
Processes all non-DNA TE classes (LINE, LTR, SINE, RC, Unknown). For each class: generates size distribution PDFs (all elements and non-zero-ORF elements only), and creates a summary Excel with the top 5 largest and top 2 ORF-containing elements per family ranked by ovotestis piRNA targeting. Saves per-family Excel files for elements with ORFs. Also produces a single master TSV covering all TE classes combined.

**Output:** `data/<class>_top5_per_repname_file_smallrna_orf.xlsx`, `data/full_<class>_tables/<repname>_file_smallrna_orf.xlsx`, `data/size_distribution_repsize_<class>.pdf`, `data/longest_insertion_size_distribution_<class>.pdf`, `data/all_classes_top5_per_repname_file_smallrna_orf.tsv`

---

### `data/explore_pirna_cluster/bw_extr.sh`
Loops over all 267 piRNA clusters in `piRNA_clusters_DL_CZ_fixed.sorted.bed` and extracts per-sample bedGraph coverage for each cluster region from the WGS BigWig files. Output is organised into per-cluster subdirectories.

**Output:** `data/explore_pirna_cluster/<cluster_id>/<sample>.bedGraph`

---

### `pirna_clusters/20260518_pirna_clusters_exploration.R`
Loads per-cluster bedGraph files and annotates them with piRNA cluster coordinates using `foverlaps`. Expands bedGraph intervals to per-base coverage and plots coverage profiles per cluster.

---

### `te_cleaning/20260517_table_clean_preparation.R`
General table preparation for DNA transposons — in progress.

---

### `pirna_clusters/blastn_pirnas_vs_genomes.nf`
Nextflow DSL2 pipeline that BLASTs the top 20 piRNA cluster sequences against all *D. laeve* genome assemblies in parallel. For each genome it builds a nucleotide BLAST database, then runs `blastn` (e-value ≤ 1e-5, ≥ 60% identity, up to 30 hits per query, 4 threads per job). With a 32-core machine all 11 genomes run 8 at a time, completing in ~1–2 hours. Supports `-resume` to restart from any failed job without rerunning completed ones.

**Run:** `nextflow run blastn_pirnas_vs_genomes.nf` (from `pirna_clusters/`)

**Config:** `pirna_clusters/nextflow.config` — set `executor.cpus` to match available cores

**Output:** `data/explore_pirna_cluster/blastn_vs_genomes/blastn_top20pirnas_vs_<genome>.tsv` — tabular format 6 with columns: `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs`
