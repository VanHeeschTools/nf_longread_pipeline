include { nanoplot } from '../modules/local/nanoplot/main'
include { pychopper } from '../modules/local/pychopper/main'

workflow QC {
    take:
    reads // Input reads from fastq
    skip_pychopper // Whether to skip pychopper (for direct RNA-seq)

    main:

    nanoplot(reads,
            params.nanoplot_extra_opts)
    nanoplot_logs = nanoplot.out

    if (skip_pychopper) {
        pychopper_logs = Channel.empty()
        full_length_reads = reads
    } else {
        //Check if pychopper custom primers are provided
        //Use the kit if no custom primers are provided
        def primer_opts = params.custom_primers_file ? "-b ${params.custom_primers_file}" : "-k ${params.cdna_kit}"
        
        pychopper(reads,
                    primer_opts,
                    params.pychopper_backend,
                    params.pychopper_extra_opts)
        pychopper_logs = pychopper.out.stats
        full_length_reads = pychopper.out.full_length_reads
    } 

    emit:
    // Main output: trimmed and oriented full length reads
    full_length_reads = full_length_reads

    // Logs for MultiQC
    nanoplot_logs = nanoplot_logs
    pychopper_logs = pychopper_logs

    // Versions
    //pychopper_versions = pychopper.out.versions
}
