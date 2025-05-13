include { nanoplot } from '../modules/local/nanoplot/main'
include { pychopper } from '../modules/local/pychopper/main'

workflow QC {
    take:
    reads // Input reads from fastq

    main:
    nanoplot(reads,
            params.nanoplot_extra_opts)

    //Check if pychopper custom primers are provided
    //Use the kit if no custom primers are provided
    def primer_opts = params.custom_primers_file ? "-b ${params.custom_primers_file}" : "-k ${params.cdna_kit}"
    
    pychopper(reads,
                primer_opts,
                params.pychopper_backend,
                params.pychopper_extra_opts)

    emit:
    // Main output: trimmed and oriented full length reads
    full_length_reads = pychopper.out.full_length_reads

    // Logs for MultiQC
    nanoplot_logs = nanoplot.out
    pychopper_logs = pychopper.out.stats

    // Versions
    pychopper_versions = pychopper.out.versions
}
