include { create_minimap2_index; minimap2_transcriptome } from '../modules/local/minimap2/main'
include { salmon; salmon_tables } from '../modules/local/salmon/main'
include { versions } from '../modules/local/versions/main'

workflow EXPRESSION {
    take:
    full_length_reads
    transcriptome_fasta  
    annotation

    main:
    // Create empty channel for versions
    ch_versions = Channel.empty()
    
     // Map against transcriptome
    create_minimap2_index(transcriptome_fasta,
                            params.minimap_index_extra_opts)
    index_ch = create_minimap2_index.out.index.first()

    ch_versions = ch_versions.mix(create_minimap2_index.out.versions)

    minimap2_transcriptome(full_length_reads,
                            index_ch,
                            params.minimap_extra_opts)
    ch_versions = ch_versions.mix(minimap2_transcriptome.out.versions)

    // Run salmon quant
    salmon(minimap2_transcriptome.out.minimap2_transcriptome_bam,
            transcriptome_fasta.first(),
            params.salmon_extra_opts)
    ch_versions = ch_versions.mix(salmon.out.versions)

     // Write the paths of the salmon_quasi output files to a text file
    quant_paths = salmon.out.quant
        .map { _sample, it -> it.toString() }
        .collectFile(
        name: 'quant_paths.txt',
        newLine: true, sort: true )

    // Salmon output statistics tables
    salmon_tables(quant_paths, annotation, "salmon_tables", params.min_tpm)

    emit:
    salmon_quant = salmon.out.quant
    salmon_multiqc = salmon_tables.out.salmon_summary
    versions = ch_versions.collect()
}     
