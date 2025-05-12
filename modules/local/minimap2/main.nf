process MINIMAP2 {
    /*
    Map reads to reference genome.
    */
    label 'minimap2'
    label 'process_high'

    input:
    tuple val(sample), path(reads)
    path reference
    val extra_opts

    output:
    tuple val(sample), path("*.sam"), emit: sam
    path "versions.yml", emit: versions

    script:
    """
    echo "Nextflow requested CPUs: ${task.cpus}"
    echo "Nextflow requested memory: ${task.memory}"

    minimap2 \
        -ax splice \
        -t $task.cpus \
        $extra_opts \
        $reference \
        $reads > ${sample}_aligned.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process MINIMAP2_INDEX {
    /*
    Build minimap index from custom transcriptome
    */
    label 'remap'
    label 'minimap2'
    label 'process_high'

    input:
        path reference
        val extra_opts
    output:
        path "transcriptome_index.mmi", emit: index
        path "versions.yml", emit: versions

    script:
    """
    minimap2 -t "${task.cpus}" ${extra_opts} -I 1000G -d "transcriptome_index.mmi" "${reference}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}


process MINIMAP2_TRANSCRIPTOME{
    /*
    Map reads to custom transcriptome.
    */
    label 'remap'
    label 'minimap2'
    label 'process_high'

    input:
       tuple val(sample), path (fastq_reads)
       path index 
       val extra_opts
    output:
       tuple val(sample), path("${sample}_transcripts_aligned.sam"), emit: sam
        path "versions.yml", emit: versions


    """
    minimap2 -t ${task.cpus} -ax map-ont ${extra_opts} -N 100 ${index} ${fastq_reads} > ${sample}_transcripts_aligned.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
