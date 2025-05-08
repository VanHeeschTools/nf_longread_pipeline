process PROCESS_ALIGNMENT {
    label 'samtools'
    label 'process_high'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("*.bam"), emit: bam
    tuple val(sample), path("*.bam.bai"), emit: bai
    path("*_mapping.stats"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    samtools sort -@ $task.cpus -o ${sample}.bam $sam

    samtools index ${sample}.bam

    samtools stats ${sample}.bam > ${sample}_mapping.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process PROCESS_ALIGNMENT_TRANSCRIPTOME {
    label 'samtools'
    label 'remap'
    label 'process_medium'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("*.bam"), emit: bam
    path "versions.yml", emit: versions


    script:
    """
    samtools view -@ ${task.cpus} -Sb ${sam} > "${sample}.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}