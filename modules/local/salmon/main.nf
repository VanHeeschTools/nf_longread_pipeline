process SALMON {
    // Count transcripts using Salmon.
    // library type is specified as forward stranded (-l SF) as it should have either been through pychopper or come from direct RNA reads.
    label "process_medium"

    input:
        tuple val(sample), path(bam)
        path ref_transcriptome
        val extra_opts
    output:
         tuple val(sample), path("${sample}"), emit: quant
         path "versions.yml", emit: versions
         
    """
    salmon quant ${extra_opts} -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o ${sample}
        # Generate versions.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/^salmon //; s/Last.*\$//')
    END_VERSIONS
    """
}