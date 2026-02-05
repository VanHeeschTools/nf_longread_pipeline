include { minimap2 } from '../modules/local/minimap2/main'
include { stringtie; stringtie_summary; write_output_samplesheet } from '../modules/local/stringtie/main'
include { seqkit_stats } from '../modules/local/process_alignment/main'
include { make_gtf_list; merge_gtfs; parse_tracking; filter_annotate; transcriptome_fasta } from '../modules/local/gffcompare/main'

workflow ASSEMBLY {
    take:
    reads            // Trimmed and oriented reads
    reference_genome // Reference genome
    annotation       // Reference gtf

    main: 
    // Create empty channel for versions
    ch_versions = Channel.empty()

    minimap2(reads,
                reference_genome,
                params.minimap_extra_opts)
    ch_versions = ch_versions.mix(minimap2.out.versions)

    
    // Only keep the bam file outputs of minimap2
    bam_files = minimap2.out.minimap2_bam
        .map { _sample_id, file -> file }.collect()
    seqkit_stats(bam_files)

    ch_versions = ch_versions.mix(seqkit_stats.out.versions)

    // Run StringTie
    stringtie(minimap2.out.minimap2_bam,
                annotation,
                params.stringtie_extra_opts)
    ch_versions = ch_versions.mix(stringtie.out.versions)

    // Obtain StringTie output stats for MultiQC
    stringtie_summary(stringtie.out.gff_paths.collect(), annotation)    

    // Create tuple containing sample_id and the location of StringTie output gtfs in output directory
    minimap2_meta = minimap2.out.minimap2_bam
        .map { sample, bam -> [sample, "${params.outdir}/minimap2/${bam.name}"]}
        .collect(flat:false)
    stringtie_meta = stringtie.out.stringtie_gff
        .map { sample, gtf -> [sample, "${params.outdir}/stringtie/${gtf.name}"]}
        .collect(flat:false)

    write_output_samplesheet(minimap2_meta,stringtie_meta)

    
    // Set masked_fasta if present
    // Point to assets/NO_FILE if not set for proper path reading inside process
    masked_fasta = params.masked_fasta
        ? params.masked_fasta
        : "${projectDir}/assets/NO_FILE"
    
    // Collect GTF files and create a list file
    gtf_paths = stringtie.out.gff_paths.collect().flatten()
            .map { it -> it.toString() }

    gtf_list = gtf_paths.collectFile(
    name: 'gtflist.txt',
            newLine: true, sort: true )

    // Merge all GTFs
    merge_gtfs(gtf_list, annotation, masked_fasta, params.output_prefix)
    ch_versions = ch_versions.mix(merge_gtfs.out.versions)

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
    ch_versions = ch_versions.mix(transcriptome_fasta.out.versions)


    emit:
    stringtie_mqc = stringtie_summary.out.stringtie_multiqc
    transcriptome_gtf = filter_annotate.out.filtered_gtf
    transcriptome_fasta = transcriptome_fasta.out.fasta
    mapping_logs =  minimap2.out.bam_stats.collect()
    versions = ch_versions.collect()
}