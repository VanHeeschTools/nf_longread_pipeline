#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LONGREAD } from './workflows/longread'

workflow {

    LONGREAD()

    workflow.onComplete {
        println "Workflow finished at: ${workflow.complete}"
        println "Duration: ${workflow.duration}"
        println "Succeeded: ${workflow.success}"
        println "Work dir: ${workflow.workDir}"
    }
}

