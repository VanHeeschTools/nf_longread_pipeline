include { minimap2 } from '../modules/local/minimap2/main'
include { stringtie } from '../modules/local/stringtie/main'
include { merge_gtfs; parse_tracking; filter_annotate; transcriptome_fasta } from '../modules/local/gffcompare/main'

workflow ASSEMBLY {
    take:
    reads            // Trimmed and oriented reads
    reference_genome // Reference genome
    annotation       // Reference gtf

    main: 
    minimap2(reads,
                reference_genome,
                params.minimap_extra_opts)

    stringtie(minimap2.out.minimap2_bam,
                annotation,
                params.stringtie_extra_opts)
    
    // Set masked_fasta if present
    // Point to assets/NO_FILE if not set for proper path reading inside process
    masked_fasta = params.masked_fasta
        ? params.masked_fasta
        : "${projectDir}/assets/NO_FILE"
    
    // Collect GTF files and create a list file
    ch_gtf_list = stringtie.out.stringtie_gff.map { it[1] }.collect().map { gtfs ->
        def gtf_list = file("${workDir}/gtf_list.txt")
        gtf_list.text = gtfs.join('\n')
        return gtf_list
    }

    // Merge all GTFs
    merge_gtfs(ch_gtf_list, annotation, masked_fasta, params.output_prefix)

    // Parse the tracking file into transcript presence/absence in each sample
    parse_tracking(merge_gtfs.out.tracking, params.output_prefix)

    // Filter anotation 
    // TODO require GTF
    filter_annotate(annotation,
                    params.refseq_gtf ?: "",
                    merge_gtfs.out.merged_gtf,
                    merge_gtfs.out.tracking, 
                    params.min_occurrence,
                    params.min_tpm,
                    params.output_prefix)

    transcriptome_fasta(filter_annotate.out.filtered_gtf,
                        reference_genome,
                        params.output_prefix)

    emit:
    transcriptome_gtf = filter_annotate.out.filtered_gtf
    transcriptome_fasta = transcriptome_fasta.out.fasta
    mapping_logs =  minimap2.out. bam_stats.collect()
}