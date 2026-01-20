// Run fusion detection with JAFFAL
process jaffal {
    label 'process_extreme'

    input:
        path full_length_reads // Path, location of input files
        val jaffal_data_dir    // Path, location of directory with required annotation files
        val genome_version     // String, genome version
        val annotation_version // String, annotation version

    output:
        path "*" // Remove after testing
        path "jaffa_results.csv", emit: jaffa_results_csv
        path "jaffa_results.fasta", emit: jaffa_results_fasta

    when:
        task.ext.when == null || task.ext.when 

    script:
        """
        # Run Jaffal
        bpipe run \
            -p inputPathsAreAbsolute=true \
            -p genome=${genome_version} \
            -p annotation=${annotation_version} \
            -p refBase=${jaffal_data_dir} \
            /JAFFA/JAFFAL.groovy \
            ${full_length_reads.join(' ')}
        """
}
