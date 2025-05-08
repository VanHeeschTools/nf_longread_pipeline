process STRINGTIE {
    tag "$sample"
    label 'process_high'

    input:
    tuple val(sample), path(bam)
    path reference
    val extra_opts

    output:
    tuple val(sample), path("*.gff"), emit: gff
    path "*_stringtie.log", emit: log
    path "versions.yml", emit: versions

    script:
    def reference_command = reference.name != 'NO_FILE' ? "-G $reference" : ''
    
    """
    stringtie \\
        $reference_command \\
        --rf -L -v \\
        -p $task.cpus \\
        $extra_opts \\
        -o ${sample}.gff \\
        -l $sample \\
        $bam 2> ${sample}_stringtie.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}