#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// -- Define pipeline parameters --
params.reads = "$baseDir/data/gg_ko_SRR1039860{8,9}_[12].fastq.gz" // Path to input FASTQ files
params.transcriptome = "$baseDir/data/c_elegans.PRJEB28388.WBPS18.transcriptome.fa" // Path to transcriptome index
params.outdir = "$baseDir/results"

log.info """
         S I N G L E - C E L L - P I P E L I N E
         =======================================
         Reads         : ${params.reads}
         Transcriptome : ${params.transcriptome}
         Output Dir    : ${params.outdir}
         """
         .stripIndent()

// -- Define the workflow --
workflow {
    // Channel to feed FASTQ files into the pipeline
    // groupTuple groups paired-end reads by their sample ID
    Channel
        .fromFilePairs(params.reads)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { ch_reads }

    // Run processes
    FASTQC(ch_reads)
    KALLISTO(ch_reads, file(params.transcriptome))
}

// -- Process 1: Quality Control --
process FASTQC {
    publishDir "$params.outdir/fastqc", mode: 'copy'
    conda "envs/qc.yaml" // Specifies the Conda environment

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_out/${sample_id}/*"

    script:
    """
    mkdir -p fastqc_out/${sample_id}
    fastqc -o fastqc_out/${sample_id} --noextract ${reads[0]} ${reads[1]}
    """
}

// -- Process 2: Quantification --
process KALLISTO {
    publishDir "$params.outdir/kallisto", mode: 'copy'
    conda "envs/quant.yaml" // Specifies the Conda environment

    input:
    tuple val(sample_id), path(reads)
    path transcriptome

    output:
    path sample_id

    script:
    """
    kallisto quant -i ${transcriptome} -o ${sample_id} --single -l 200 -s 20 ${reads[0]} ${reads[1]}
    """
}

DEEP_ANALYSIS(KALLISTO.out.collect())
process DEEP_ANALYSIS {
    publishDir "$params.outdir/deep_analysis", mode: 'copy'
    conda "envs/analysis.yaml"

    input:
    path kallisto_outputs

    output:
    path "*.h5ad"
    path "*.png"

    script:
    """
    python ${baseDir}/scripts/downstream_analysis.py ${kallisto_outputs[0]}
    """
}