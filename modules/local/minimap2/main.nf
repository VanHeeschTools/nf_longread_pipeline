// Map reads to reference genome using Minimap2
process minimap2 {
    label 'minimap2'
    label 'process_high'

    input:
        tuple val(sample), path(reads)
        path reference
        val extra_opts

    output:
        tuple val(sample), path("${sample}_aligned.sorted.bam"), emit: minimap2_bam
        tuple val(sample), path("${sample}_aligned.sorted.bam.bai"), emit: minimap2_bam_bai
        path "${sample}_mapping.stats", emit: bam_stats
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
        """
        echo "Nextflow requested CPUs: ${task.cpus}"
        echo "Nextflow requested memory: ${task.memory}"

        minimap2 \
            -ax splice \
            -t $task.cpus \
            $extra_opts \
            $reference \
            $reads | \
        samtools sort \
            -@ ${task.cpus} \
            -o ${sample}_aligned.sorted.bam

        samtools index ${sample}_aligned.sorted.bam

        samtools stats ${sample}_aligned.sorted.bam > ${sample}_mapping.stats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}

// Build minimap index from custom transcriptome
process create_minimap2_index {

    label 'remap'
    label 'minimap2'
    label 'process_high'

    input:
        path reference
        val extra_opts
    output:
        path "transcriptome_index.mmi", emit: index
        path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
        """
        minimap2 -t "${task.cpus}" ${extra_opts} -I 1000G -d "transcriptome_index.mmi" "${reference}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
        END_VERSIONS
        """
}

// Map reads to custom transcriptome using Minimap2
process minimap2_transcriptome{
    label 'remap'
    label 'minimap2'
    label 'process_high'

    input:
        tuple val(sample), path (fastq_reads)
        path index 
        val extra_opts

    output:
        tuple val(sample), path("${sample}_transcripts_aligned.bam"), emit: minimap2_transcriptome_bam
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        minimap2 -t ${task.cpus} \
            -ax map-ont ${extra_opts} \
            -N 100 ${index} \
            ${fastq_reads} |
        samtools view \
            -@ ${task.cpus} \
            -Sb  \
            -o ${sample}_transcripts_aligned.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}
