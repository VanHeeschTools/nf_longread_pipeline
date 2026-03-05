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
        path("${sample}.gff"), emit: gff_paths
        path "${sample}_stringtie.log", emit: log
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def reference_command = reference_genome.name != 'NO_FILE' ? "-G $reference_genome" : ''
        
        """
        stringtie \
            $reference_command \
            --rf -L -v \
            -p $task.cpus \
            $extra_opts \
            -o ${sample}.gff \
            -l $sample \
            $bam 2> ${sample}_stringtie.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            StringTie: \$(stringtie --version 2>&1)
        END_VERSIONS
        """
}


process stringtie_summary {
    label 'process_superlow'
    label 'gffcompare'

    input:
        path gff_list       // Path, list of StringTie gff output files
        path reference_gtf   // Val, reference gtf path

    output:
        path "all_samples_stringtie_counts_mqc.tsv", emit: stringtie_multiqc

    script:
        """
        OUT="all_samples_stringtie_counts_mqc.tsv"
        echo -e "Sample\tTranscripts\tExons\tKnown_transcripts\tNovel_transcripts" >> "\$OUT"

        # Iterate over StringTie output gtfs to obtain statistics when comparing to reference gff
        for GFF in ${gff_list.join(' ')}; do
            # Extract sample_id from filename
            SAMPLE=\$(basename "\$GFF" .gff)

            transcripts=\$(awk '\$3=="transcript"' "\$GFF" | wc -l)
            exons=\$(awk '\$3=="exon"' "\$GFF" | wc -l)

            gffcompare -r "${reference_gtf}" -o "\${SAMPLE}_gffcmp" "\$GFF"
            ann="\${SAMPLE}_gffcmp.annotated.gtf"
            
            all=\$(grep -c \$'\ttranscript\t' "\$ann")
            known=\$(grep 'class_code "[=c]"' "\$ann" | grep -c \$'\ttranscript\t')
            novel=\$((all - known))

            echo -e "\$SAMPLE\t\$transcripts\t\$exons\t\$known\t\$novel" >> "\$OUT"
        done
        """
}


// Create samplesheet showing id, location of StringTie gtf in output samplesheet and data type (longread)
process write_output_samplesheet{
    label 'process_superlow'
    label 'python'

    input:
        val minimap2_meta     // Val, string containing sample id and location of Minimap2 output BAM in output directory
        val gtf_file_location // Val, string containing sample id and location of StringTie output gtf in output directory

    output:
        path "output_samplesheet.csv", emit: output_samplesheet

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    printf '%s\n' \
    'sample_id,biomaterial_id,file,file_type,disease_state,seq_type' \
    ${minimap2_meta.collect { r -> "'${r[0]},null,${r[1]},bam,null,longread'" }.join(' ')} \
    ${gtf_file_location.collect { r -> "'${r[0]},null,${r[1]},gtf,null,longread'" }.join(' ')} \
    > output_samplesheet.csv
    """


}
