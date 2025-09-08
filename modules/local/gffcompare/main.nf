process GFFCOMPARE {
    label 'gffcompare'
    label 'process_medium'

    input:
    tuple val(sample), path(gtf)
    path reference_gtf
    path masked_fasta


    output:
    tuple val(sample), path("${sample}.annotated.gtf"), emit: annotated_gtf
    path "${sample}_gffcompare.stats", emit: stats
    path "${sample}.tracking", emit: tracking
    path "versions.yml", emit: versions

    script:
    def gtf_command = reference_gtf.name != 'NO_FILE' ? "-r $reference_gtf" : ''
    def masked_fasta_command = masked_fasta.name != 'NO_FILE' ? "-s $masked_fasta" : ''

    """
    gffcompare \
        -V \
        ${gtf_command} \
        ${masked_fasta_command} \
        -o "${sample}" \
        $gtf

    mv ${sample}.stats ${sample}_gffcompare.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(gffcompare --version 2>&1 | sed 's/^gffcompare v//')
    END_VERSIONS
    """
}

process MERGE_GTFS {
    label 'gffcompare'
    label 'merge_gtfs'
    label 'process_medium'

    input:
    path gtf_list
    path reference_gtf
    path masked_fasta
    val output_prefix

    output:
    path "*.gtf", emit: merged_gtf
    path "*.stats", emit: stats
    path "*.tracking", emit: tracking
    path "versions.yml", emit: versions

    script:
    def gtf_command = reference_gtf.name != 'NO_FILE' ? "-r $reference_gtf" : ''
    def masked_fasta_command = masked_fasta.name != 'NO_FILE' ? "-s $masked_fasta" : ''

    """
    gffcompare \
        -V \
        ${gtf_command} \
        ${masked_fasta_command} \
        -o "${output_prefix}" \
        -i "${gtf_list}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(gffcompare --version 2>&1 | sed 's/^gffcompare v//')
    END_VERSIONS
    """
}

process PARSE_TRACKING {
    label 'merge_gtfs'
    label 'python'

    input:
    path tracking_file
    val output_prefix

    output:
    path "${output_prefix}_transcript_presence.tsv"

    script:
    """
    python ${projectDir}/bin/parse_tracking.py ${tracking_file} ${output_prefix}
    """
}

// Define process for transcript filtering and annotation
process FILTER_ANNOTATE {
    label "merge_gtfs"
    label "process_medium"

    input:
    val reference_gtf   // Path to the input reference gtf file
    val refseq_gtf     // Path to the refseq gtf file (optional)
    path gtf_novel      // Path to the merged gtf file
    path gtf_tracking   // Path to the tracking file created by the merge step
    val min_occurrence  // Val contatining the minimum occurence of transcripts for filtering (defaults to 1)
    val min_tpm         // Val containing the minium tpm of transcripts for filtering (defaults to 0.1)
    val output_prefix // Val containing output basename

    output:
    path "${output_prefix}_novel_filtered.gtf", emit: filtered_gtf
    path "${output_prefix}_novel_filtered.log"
    path "${output_prefix}_novel_filtered.tsv"

    script:
    """
    filter_annotate.R \
    "${reference_gtf}" \
    "${gtf_novel}" \
    "${gtf_tracking}" \
    "${min_occurrence}" \
    "${min_tpm}" \
    "${output_prefix}_novel_filtered" \
    "${projectDir}/bin/" \
    "${refseq_gtf}"
    """
}


// Creates a fasta file of the transcript sequence using the reference fasta file and the transcriptome gtf
process TRANSCRIPTOME_FASTA {
    label "merge_gtfs"
    label "process_low"

    input:
    val gtf    // Merged and filtered transcriptome file
    path fasta  // Path to input reference fasta file
    val prefix

    output:
    path "*_transcriptome.fa", emit: fasta
    path "versions.yml", emit: versions


    script:
    """
    gffread -w ${prefix}_transcriptome.fa -g ${fasta} ${gtf}

    # Generate versions.yml
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version 2>&1)
    END_VERSIONS
    """
}