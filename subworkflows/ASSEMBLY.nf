include { MINIMAP2 } from '../modules/local/minimap2/main'
include { PROCESS_ALIGNMENT } from '../modules/local/process_alignment/main'
include { STRINGTIE } from '../modules/local/stringtie/main'
include { MERGE_GTFS; GFFCOMPARE; PARSE_TRACKING; FILTER_ANNOTATE; TRANSCRIPTOME_FASTA } from '../modules/local/gffcompare/main'

workflow ASSEMBLY {
    take:
    reads // Trimmed and oriented reads
    reference // Reference genome
    annotation // Reference gtf

    main: 
    MINIMAP2(reads,
                reference,
                params.minimap_extra_opts)

    PROCESS_ALIGNMENT(MINIMAP2.out.sam)

    STRINGTIE(PROCESS_ALIGNMENT.out.bam,
                annotation,
                params.stringtie_extra_opts)
    
    // Set masked_fasta if present
    // Point to assets/NO_FILE if not set for proper path reading inside process
    masked_fasta = params.masked_fasta
        ? params.masked_fasta
        : "${projectDir}/assets/NO_FILE"
    
    // Run GFFCOMPARE on each sample's GTF
    GFFCOMPARE(STRINGTIE.out.gff, annotation, masked_fasta)

    // Collect GTF files and create a list file
    ch_gtf_list = STRINGTIE.out.gff.map { it[1] }.collect().map { gtfs ->
        def gtf_list = file("${workDir}/gtf_list.txt")
        gtf_list.text = gtfs.join('\n')
        return gtf_list
    }

    // Merge all GTFs
    MERGE_GTFS(ch_gtf_list, annotation, masked_fasta, params.output_prefix)

    // Parse the tracking file into transcript presence/absence in each sample
    PARSE_TRACKING(MERGE_GTFS.out.tracking, params.output_prefix)

    // Filter anotation 
    // TODO require GTF
    FILTER_ANNOTATE(annotation,
                    MERGE_GTFS.out.merged_gtf,
                    MERGE_GTFS.out.tracking, 
                    params.min_occurrence,
                    params.min_tpm,
                    params.output_prefix)


    TRANSCRIPTOME_FASTA(FILTER_ANNOTATE.out.filtered_gtf,
                        reference,
                        params.output_prefix)

    emit:
    transcriptome = FILTER_ANNOTATE.out.filtered_gtf
    fasta = TRANSCRIPTOME_FASTA.out.fasta
    mapping_logs = PROCESS_ALIGNMENT.out.stats.collect()
    gffcompare_logs = GFFCOMPARE.out.stats.collect()
}