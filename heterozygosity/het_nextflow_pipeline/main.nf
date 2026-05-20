#!/usr/bin/env nextflow
/*
 * Per-window heterozygosity from PacBio HiFi reads + assembly.
 * One process per sample, fanned out over rows of a CSV samplesheet.
 *
 * Input  : --samplesheet samples.csv   (columns: sample,assembly,reads)
 * Output : results/<sample>/  with BAM, VCFs, het_per_window.annotated.tsv,
 *                                  het_per_contig.tsv, flagstat
 */

nextflow.enable.dsl = 2

// -------- parameters (override on CLI: --window 50000) ----------------------
params.samplesheet = "${projectDir}/samplesheet.csv"
params.outdir      = "${projectDir}/results"
params.window      = 100000
params.min_qual    = 20
params.min_dp      = 5
params.max_dp      = 200

// =============================================================================
process INDEX_REF {
    tag "$sample"
    input:
        tuple val(sample), path(ref), path(reads)
    output:
        tuple val(sample), path(ref), path("${ref}.fai"), path(reads), emit: indexed
    script:
    """
    samtools faidx ${ref}
    """
}

// =============================================================================
process MAP_HIFI {
    tag "$sample"
    cpus   { Math.min(32, Runtime.runtime.availableProcessors()) }
    memory '24 GB'
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy', pattern: "*.{bam,bai,flagstat.txt}"

    input:
        tuple val(sample), path(ref), path(fai), path(reads)
    output:
        tuple val(sample), path(ref), path(fai),
              path("${sample}.sorted.bam"), path("${sample}.sorted.bam.bai"), emit: bam
        path "${sample}.flagstat.txt"
    script:
    def sort_threads = Math.max(2, (task.cpus.intValue() / 4) as int)
    """
    minimap2 -t ${task.cpus} -ax map-hifi --secondary=no ${ref} ${reads} \\
      | samtools sort -@ ${sort_threads} -m 2G -o ${sample}.sorted.bam -
    samtools index -@ ${task.cpus} ${sample}.sorted.bam
    samtools flagstat ${sample}.sorted.bam > ${sample}.flagstat.txt
    """
}

// =============================================================================
process CALL_VARIANTS {
    tag "$sample"
    cpus   { Math.min(8, Runtime.runtime.availableProcessors()) }
    memory '8 GB'
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy', pattern: "*.vcf.gz*"

    input:
        tuple val(sample), path(ref), path(fai), path(bam), path(bai)
    output:
        tuple val(sample), path(ref), path(fai),
              path("${sample}.vcf.gz"), path("${sample}.vcf.gz.csi"), emit: vcf
    script:
    """
    bcftools mpileup -f ${ref} -B -Q 20 -q 20 \\
        -a 'FORMAT/AD,FORMAT/DP,INFO/AD' \\
        --threads ${task.cpus} ${bam} \\
      | bcftools call -m -v --threads ${task.cpus} -Oz -o ${sample}.vcf.gz && \
    bcftools index ${sample}.vcf.gz
    """
}

// =============================================================================
process FILTER_HET {
    tag "$sample"
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy', pattern: "*.het.vcf.gz*"

    input:
        tuple val(sample), path(ref), path(fai), path(vcf), path(csi)
    output:
        tuple val(sample), path(ref), path(fai),
              path("${sample}.het.vcf.gz"), path("${sample}.het.vcf.gz.csi"), emit: het
    script:
    """
    bcftools view -v snps \\
        -i 'QUAL>=${params.min_qual} && FMT/DP>=${params.min_dp} && FMT/DP<=${params.max_dp} && (GT="0/1" || GT="0|1" || GT="1|0")' \\
        ${vcf} -Oz -o ${sample}.het.vcf.gz
    bcftools index ${sample}.het.vcf.gz
    """
}

// =============================================================================
process COUNT_HET {
    tag "$sample"
    publishDir { "${params.outdir}/${sample}" }, mode: 'copy'

    input:
        tuple val(sample), path(ref), path(fai), path(het_vcf), path(het_csi)
    output:
        tuple val(sample),
              path("het_per_window.annotated.tsv"),
              path("het_per_contig.tsv"),
              path("windows_${params.window}.bed")
    script:
    """
    awk -v w=${params.window} 'BEGIN{OFS="\\t"} {for(i=0;i<\$2;i+=w) print \$1, i, (i+w>\$2?\$2:i+w)}' \\
        ${fai} > windows_${params.window}.bed

    bedtools intersect -a windows_${params.window}.bed -b ${het_vcf} -c > raw_counts.tsv

    {
      printf "contig\\tstart\\tend\\tlen\\thet_count\\thet_per_100kb\\n"
      awk -v w=${params.window} 'BEGIN{OFS="\\t"} {len=\$3-\$2; rate=(len>0?\$4/(len/100000):0); print \$1,\$2,\$3,len,\$4,rate}' raw_counts.tsv
    } > het_per_window.annotated.tsv

    awk 'BEGIN{OFS="\\t"; print "contig","contig_length","het_count","het_per_100kb"}
         NR>1 {tot[\$1]+=\$5; len[\$1]=(len[\$1]>\$3?len[\$1]:\$3)}
         END  {for(c in tot) printf "%s\\t%d\\t%d\\t%.3f\\n", c, len[c], tot[c], tot[c]/(len[c]/100000)}' \\
        het_per_window.annotated.tsv \\
      | (read header; echo "\$header"; sort -k2,2nr) \\
      > het_per_contig.tsv
    """
}

// =============================================================================
workflow {
    samples = Channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            tuple(row.sample,
                  file(row.assembly, checkIfExists: true),
                  file(row.reads,    checkIfExists: true))
          }

    INDEX_REF(samples)
    MAP_HIFI(INDEX_REF.out.indexed)
    CALL_VARIANTS(MAP_HIFI.out.bam)
    FILTER_HET(CALL_VARIANTS.out.vcf)
    COUNT_HET(FILTER_HET.out.het)
}
