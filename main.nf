#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LONGREAD } from './workflows/longread'
include { printHeader } from "./modules/local/helperfunctions/main.nf"
workflow {
    printHeader()
    LONGREAD()
}

