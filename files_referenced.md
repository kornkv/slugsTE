# Files Referenced in slugsTE Scripts

## Shell Scripts

### `blast/run_blastx.sh`

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/blasting_transposase_unknown_tes/KAK7092566.1.fasta` | 1MB | `d16b92f89e11e6c510d44df11b9f7803` | |
| `data/blasting_transposase_unknown_tes/20260517_DL_yahs_all_pb_final_unknowns_only.fasta` | <1MB | `e224cacd37fea266a21bcd4863c83d57` | |
| `data/blasting_transposase_unknown_tes/unknowns_db` | | | output |
| `data/blasting_transposase_unknown_tes/transposase_prot_db` | | | output |
| `data/blasting_transposase_unknown_tes/tblastn_KAK7092566_vs_unknowns.tsv` | 1MB | `01118e2dba669f1218b746bce4fde00c` | output |
| `data/blasting_transposase_unknown_tes/blastx_unknowns_vs_KAK7092566.tsv` | 1MB | `651c579bd1c414881beea6cff79a4d33` | output |

### `blast/run_blastx_maverick.sh`

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/blasting_transposase_unknown_tes/maverick_orf_protein.fasta` | 1MB | `6752cfaff3e802ad74242aa0c0e00e53` | |
| `data/blasting_transposase_unknown_tes/blastx_unknowns_vs_maverick_orf.tsv` | 1MB | `658ac66d94866bf144c75c4852b25f88` | output |

### `heterozygosity/` run scripts

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `heterozygosity/L681_samplesheet.csv` | 1MB | `901bac0d31386140ff87934f4ca6786a` | |
| `heterozygosity/frankfurt_samplesheet.csv` | 1MB | `b8dd32b8d462d15b2c8a372d941316c4` | |
| `heterozygosity/invasive_samplesheet.csv` | 1MB | `277538cc58a8e63d77ac6e79cab37671` | |
| `/data2/othersamples_samplesheet.csv` | | | **MISSING** |
| `data/heterozygosity/frankfurt/` | | | output dir |
| `data/heterozygosity/invasive/` | | | output dir |
| `/data2/heterozygosity/othersamples/` | | | output dir |

### `pirna_clusters/` run scripts

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/explore_pirna_cluster/top20_pirna_clusters_sequences.fasta` | 2MB | `ed5001bab3859dbaba0cb7fc2ad3861f` | |
| `/data/genomes/*.fasta` | | | pattern |
| `data/explore_pirna_cluster/blastn_vs_genomes/` | | | output dir |

---

## R Scripts

### `heterozygosity/20260518_explore_hets.R`

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/het_L681/L681/L681.het.vcf.gz` | | | **MISSING** |
| `data/het_L681/L681/L681_to_ref_mapping.tsv` | | | **MISSING** |
| `data/heterozygosity/regions_specific_for_aphalic.bed` | 1MB | `25abf4821221f65a6b660ef1fd575759` | output |

### `heterozygosity/20260519_exploring_results_of_het_pipeline.R`

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/heterozygosity/othersamples/**/*het_per_contig.tsv` | | | pattern |
| `data/heterozygosity/othersamples/**/*het_per_window.annotated.tsv` | | | pattern |
| `data/heterozygosity/**/*.paf` | | | see PAF table below |

### `pirna_clusters/20260518_pirna_clusters_exploration.R`

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/explore_pirna_cluster/piRNA_clusters_DL_CZ_fixed.sorted.bed` | 1MB | `a4874919d7cf4d258507a6530a2b0885` | |
| `data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv` | 282MB | `9b2b50e781a4c31dba3ddfa5ff00005b` | |
| `data/Deroceras_laeve_CZ.1.0.fasta` | 1762MB | `4d9767c908855f447d9066695b7a040a` | |
| `data/explore_pirna_cluster/top20_pirna_clusters_sequences.fasta` | 2MB | `ed5001bab3859dbaba0cb7fc2ad3861f` | output |

### `te_cleaning/` scripts

| File | Size | MD5 | Notes |
|------|------|-----|-------|
| `data/DL_yahs_all_pb_final.fasta.out` | 633MB | `3f05d0b455c745c526d15c592f80598e` | |
| `data/piRNA_counts_matrix.per_ins.with_coords.tsv` | <1MB | `d9010532046bf74f9c10a282f0655c3b` | |
| `data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated.tsv` | 229MB | `7b8a9976f0b3a4d714dff5fdcbe7330c` | output |
| `data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_orfs.tsv` | 282MB | `9b2b50e781a4c31dba3ddfa5ff00005b` | output |
| `data/DL_yahs_all_pb_final.fasta.out.smallrna_annotated_with_two_orfs.tsv` | 462MB | `fd644e57bb30644243d375e642153f19` | output |
| `data/Deroceras_laeve_CZ.1.0.fasta` | 1762MB | `4d9767c908855f447d9066695b7a040a` | |
| `data/blasting_transposase_unknown_tes/blastx_summary.csv` | 1MB | `376145211b64baa2e03031f221e768e4` | output |
| `data/blasting_transposase_unknown_tes/tblastn_KAK7092566_vs_unknowns.tsv` | 1MB | `01118e2dba669f1218b746bce4fde00c` | |
| `data/blasting_transposase_unknown_tes/blastx_unknowns_vs_KAK7092566.tsv` | 1MB | `651c579bd1c414881beea6cff79a4d33` | |

---

## PAF Files (`data/heterozygosity/*.paf`)

| File | Size | MD5 |
|------|------|-----|
| `Deroceras_laeve_CZ_1_0__vs__L691_Schlossteich_Euphallic.paf` | 19MB | `b754291d14800518ce30f6f7fec67bcc` |
| `Deroceras_laeve_CZ_1_0__vs__L451_Glasgow_Aphallic.paf` | 6MB | `7bc3ec09d9879881024c0f3749c04041` |
| `Deroceras_laeve_CZ_1_0__vs__L684_Jonathan_Aphallic.paf` | 18MB | `11816b99efb93cec31a39d53d298b3e0` |
| `Deroceras_laeve_CZ_1_0__vs__L633_Wein6_Aphallic.paf` | 18MB | `d02c86b9502f7538ae1a7c4fe3f89bae` |
| `Deroceras_laeve_CZ_1_0__vs__L532_Jonathan2_Euphallic.paf` | 19MB | `70b340755b082fdec63873d8bcaa12dd` |
| `Deroceras_laeve_CZ_1_0__vs__L681_Schlossteich2_Euphallic.paf` | 19MB | `555deb6dad7364111c58fd1e4e2e4ece` |
| `Deroceras_laeve_CZ_1_0__vs__derLae1_Mexico_Euphallic.paf` | 16MB | `368fd849d402f6d78bd449961f870ab2` |
| `Deroceras_laeve_CZ_1_0__vs__L685_Jonathan2_Euphallic.paf` | 19MB | `4f0523bdd123d35eef37e434106c312e` |


## checksum file

All other files md5sums are in checksums.md5 file!