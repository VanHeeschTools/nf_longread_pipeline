process pychopper {
    tag "$sample"
    label 'process_high'

    input:
    tuple val(sample), path(reads)
    val(primer_opts)
    val(backend)
    val(extra_opts)

    output:
    tuple val(sample), path("${sample}_full_length_reads.fastq.gz"), emit: full_length_reads
    path "${sample}_pychopper_report.pdf", emit: report
    path "${sample}_pychopper_stats.tsv", emit: stats
    path "versions.yml", emit: versions

    script:
    //Check if reads is a list
    //If not, make it a list for .join() to work
    def input_files = reads instanceof List ? reads : [reads]

    """
    # Concatenate input files if there are multiple
    cat ${input_files.join(' ')} > ${sample}_combined_input.fastq.gz

    pychopper ${primer_opts} \
    -m ${backend} \
    -r ${sample}_pychopper_report.pdf \
    -S ${sample}_pychopper_stats.tsv \
    -t ${task.cpus} \
    ${extra_opts} \
    ${sample}_combined_input.fastq.gz \
    ${sample}_full_length_reads.fastq

    # Compress the output
    gzip ${sample}_full_length_reads.fastq

    echo "Pychopper completed for sample: $sample"

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pychopper: \$(python -c "import pychopper; print(f'pychopper,{pychopper.__version__}')")
    END_VERSIONS
    """
}