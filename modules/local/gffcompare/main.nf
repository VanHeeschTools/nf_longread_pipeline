// Join StringTie output gtf paths into a single list
process make_gtf_list {
    input:
        path gtfs   //Path, contains channel of StringTie ouput gtfs

    output:
        path "gtf_list.txt"

    script:
        """
        printf "%s\n" ${gtfs.join(' ')} > gtf_list.txt
        """
}

// Run gffcompare on list of StringTie output gtfs
process merge_gtfs {
    label 'gffcompare'
    label 'merge_gtfs'
    label 'process_medium'

    input:
        path gtf_list       // Path, file containing paths to StringTie output gtfs
        path reference_gtf  // Path, reference gtf file
        path masked_fasta   // Path, masked genome fasta file
        val output_prefix   // Val, string containing prefix for output

    output:
        path "${output_prefix}.combined.gtf", emit: merged_gtf
        path "${output_prefix}.stats", emit: stats
        path "${output_prefix}.tracking", emit: tracking
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def gtf_command = reference_gtf.name != 'NO_FILE' ? "-r $reference_gtf" : ''
        def masked_fasta_command = masked_fasta.name != 'NO_FILE' ? "-s $masked_fasta" : ''

        """
        ls *.gff > gtflist.txt
        gffcompare \
            -V \
            ${gtf_command} \
            ${masked_fasta_command} \
            -o "${output_prefix}" \
            -i gtflist.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffcompare: \$(gffcompare --version 2>&1 | sed 's/^gffcompare v//')
        END_VERSIONS
        """
}

process parse_tracking {
    label 'merge_gtfs'
    label 'python'

    input:
        path tracking_file
        val output_prefix

    output:
        path "${output_prefix}_transcript_presence.tsv"

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        parse_tracking.py ${tracking_file} ${output_prefix}
        """
}

// Define process for transcript filtering and annotation
process filter_annotate {
    label "merge_gtfs"
    label "process_medium"

    input:
        path reference_gtf   // Path, input reference gtf file
        path refseq_gtf      // Path, refseq gtf file (optional)
        path gtf_novel      // Path, merged gtf file
        path gtf_tracking   // Path, tracking file created by the merge step
        val min_occurrence  // Val, minimum occurence of transcripts for filtering (defaults to 1)
        val min_tpm         // Val, minium tpm of transcripts for filtering (defaults to 0.1)
        val output_prefix   // Val, output basename
        path "filter_annotate.R"
        path "filter_annotate_functions.R"

    output:
        path "${output_prefix}.filtered.extended_reference.gtf", emit: filtered_gtf
        path "${output_prefix}.filtered.novel_transcripts.gtf", emit: novel_gtf
        path "${output_prefix}.filtered.log", emit: filtered_log
        path "${output_prefix}.filtered.tsv", emit: filtered_tsv

    when:
        task.ext.when == null || task.ext.when

    script:
        def refseq_arg = refseq_gtf ? "${refseq_gtf}" : ""
        """
        filter_annotate.R \
        "${reference_gtf}" \
        "${gtf_novel}" \
        "${gtf_tracking}" \
        "${min_occurrence}" \
        "${min_tpm}" \
        "${output_prefix}.filtered" \
        "${projectDir}/bin/" \
        "${refseq_arg}"
        """
}


// Creates a fasta file of the transcript sequences using the reference fasta file and the transcriptome gtf
process transcriptome_fasta {
    label "merge_gtfs"
    label "process_low"

    input:
        path gtf     // Merged, and filtered transcriptome file
        path fasta  // Path, to input reference fasta file
        path fai
        val prefix  // Val, string containing output prefix

    output:
        path "${prefix}_transcriptome.fa", emit: fasta
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        gffread \\
            -w ${prefix}_transcriptome.fa \\
            -g ${fasta} \\
            ${gtf}

        # Generate versions.yml
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gffread: \$(gffread --version 2>&1)
        END_VERSIONS
        """
}
