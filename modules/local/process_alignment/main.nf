process process_alignment {
    label 'samtools'
    label 'process_high'

    input:
    tuple val(sample), path(sam)

    output:
    path "${sample}_mapping.stats", emit: stats
    path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    //TODO filter unmapped reads
    script:
        """
        samtools stats ${sample}_aligned.sorted.bam > ${sample}_mapping.stats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

//TODO add process generate alignment stats with sekqit

process process_alignment_transcriptome {
    label 'samtools'
    label 'remap'
    label 'process_medium'

    input:
        tuple val(sample), path(sam)

    output:
        tuple val(sample), path("*.bam"), emit: bam
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools view -@ ${task.cpus} -Sb ${sam} > "${sample}.bam"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}
