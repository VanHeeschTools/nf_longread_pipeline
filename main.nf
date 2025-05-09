#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LONGREAD } from './workflows/longread'

workflow {
    if (params.input) {
        if (file(params.input).isDirectory()) {
            log.info "Processing directory: ${params.input}"
            
            input_ch = Channel
                .fromPath("${params.input}/**/*.{fastq,fastq.gz,bam}")
                .map { file -> 
                    def parent = file.parent.name
                    def sample = params.sample ?: (parent == params.input.tokenize('/')[-1] ? 'sample' : parent)
                    return tuple(sample, file)
                }
                .groupTuple()
                .map { sample, files -> 
                    if (files.isEmpty()) {
                        error "No FASTQ/BAM files found for sample: ${sample}"
                    }
                    return tuple(sample, files)
                }
                .ifEmpty { error "No input files found in directory: ${params.input}" }
        } else {            
            input_ch = Channel
                .fromPath(params.input)
                .map { file -> 
                    if (!file.exists()) {
                        error "Input file does not exist: ${file}"
                    }
                    log.info "Processing single file: ${file.simpleName}"
                    return tuple(params.sample ?: file.simpleName, file)
                }
        }
    } else {
        error "Input not specified. Please provide --input parameter."
    }

    // Sample sheet handling with error check
    if (params.sample_sheet) {
        if (!file(params.sample_sheet).exists()) {
            error "Sample sheet file does not exist: ${params.sample_sheet}"
        }
        log.info "Using sample sheet: ${params.sample_sheet}"
        sample_sheet_ch = Channel.fromPath(params.sample_sheet)
    } else {
        log.warn "No sample sheet provided, continuing without it."
        sample_sheet_ch = Channel.empty()
    }

    LONGREAD(input_ch, sample_sheet_ch)
}

workflow.onComplete {
    println "Workflow finished at: ${workflow.complete}"
    println "Duration: ${workflow.duration}"
    println "Succeeded: ${workflow.success}"
    println "Work dir: ${workflow.workDir}"
}