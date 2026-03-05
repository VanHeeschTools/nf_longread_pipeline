process salmon {
    // Count transcripts using Salmon.
    // library type is specified as forward stranded (-l SF) as it should have either been through pychopper or come from direct RNA reads.
    tag "$sample"
    label "process_medium"

    input:
        tuple val(sample), path(bam) // Tuple, containing sample id and minimap2 output bam file
        path ref_transcriptome       // Path, location of transcriptome fasta file
        val extra_opts               // Val, potential extra parameters given in config file

    output:
        tuple val(sample), path("${sample}/quant.sf"), emit: quant
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:    
        """
        salmon quant ${extra_opts} -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o ${sample}
        
        # Generate versions.yml
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(salmon --version | sed 's/^salmon //; s/Last.*\$//')
        END_VERSIONS
        """
}

// Create satistics table using the salmon output
process salmon_tables {

    label "process_low"
    containerOptions '--entrypoint='

    input:
        tuple val(ids), path("inputs/*") // Val, string containing all paths to salmon quant output files
        path gtf        // Path, reference gtf file
        val prefix      // Val, string of output prefix
        val min_tpm     // Val, min_tpm for statistics

    output:
        path "${prefix}*"
        path "${prefix}_multiqc_summary_mqc.tsv", emit: salmon_summary

    script:
        """
        # Place all file names into txt file for the R script to parse
        find -L inputs/ -name "quant.sf" > local_quant_paths.txt
        
        salmon_cohort_tables.R \
        local_quant_paths.txt \
        ${gtf} \
        ${prefix} \
        ${min_tpm}
        """
}
