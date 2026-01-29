// Run seqkit to get statistics from minimap2 output BAM files
process seqkit_stats {
    label 'process_medium'

    input:
        val minimap2_bams // Val, string containing all paths to minimap2 output bam files

    output:
        path "minimap2_bams_seqkit_stats.tsv", emit: seqkit_stats
        path "versions.yml", emit:versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        """
        # Run seqkit bam stats and write std error to ouput file, removes file paths
        seqkit bam \
        -j ${task.cpus} \
        -s \
        ${minimap2_bams.join(' ')} \
        2>&1 | awk 'NR==1{print; next} {sub(".*/","",\$NF); print}' \
        > minimap2_bams_seqkit_stats.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            \$(seqkit version)
        END_VERSIONS
        """
}