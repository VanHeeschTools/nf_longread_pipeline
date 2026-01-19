// StringTie process using input BAM file and reference genome
process stringtie {
    tag "$sample"
    label 'process_high'

    input:
        tuple val(sample), path(bam) // Tuple, sample id and BAM file
        path reference_genome        // Path, reference genome
        val extra_opts               // Val, extra options defined in nextflow.config

    output:
        tuple val(sample), path("${sample}.gff"), emit: stringtie_gff
        path "${sample}_stringtie.log", emit: log
        path "versions.yml", emit: versions

    script:
        def reference_command = reference_genome.name != 'NO_FILE' ? "-G $reference_genome" : ''
        
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