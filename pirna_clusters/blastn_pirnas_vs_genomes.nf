#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.query   = "/data/scripts/slugsTE/data/explore_pirna_cluster/top20_pirna_clusters_sequences.fasta"
params.genomes = "/data/genomes/*.fasta"
params.outdir  = "/data/scripts/slugsTE/data/explore_pirna_cluster/blastn_vs_genomes"
params.evalue  = "1e-5"
params.perc_id = 60
params.max_hits = 30
params.threads  = 4

process MAKEBLASTDB {
    tag "$name"

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path("${name}_db.*")

    script:
    """
    makeblastdb \
        -in $genome \
        -dbtype nucl \
        -out ${name}_db \
        -parse_seqids
    """
}

process BLASTN {
    tag "$name"
    publishDir params.outdir, mode: 'copy'

    cpus params.threads

    input:
    tuple val(name), path(db_files), path(query)

    output:
    path "blastn_top20pirnas_vs_${name}.tsv"

    script:
    def db_prefix = db_files[0].toString().replaceAll(/\.[^.]+$/, '')
    """
    blastn \
        -query $query \
        -db \$(ls *.nhr | sed 's/\\.nhr//') \
        -evalue ${params.evalue} \
        -perc_identity ${params.perc_id} \
        -max_target_seqs ${params.max_hits} \
        -num_threads ${task.cpus} \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \
        -out blastn_top20pirnas_vs_${name}.tsv
    """
}

workflow {
    query_ch = file(params.query)

    genome_ch = Channel
        .fromPath(params.genomes)
        .map { fasta -> [ fasta.baseName, fasta ] }

    db_ch = MAKEBLASTDB(genome_ch)

    BLASTN(db_ch.combine([query_ch]))
}
