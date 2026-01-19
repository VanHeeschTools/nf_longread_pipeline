include { minimap2_index; minimap2_transcriptome } from '../modules/local/minimap2/main'
include { process_alignment_transcriptome } from '../modules/local/process_alignment/main'
include { salmon } from '../modules/local/salmon/main'
include { versions } from '../modules/local/versions/main'

workflow EXPRESSION {
    take:
    full_length_reads
    transcriptome_fasta   

    main:
    // Create empty channel for versions
    ch_versions = Channel.empty()
    
     // Map against transcriptome
    minimap2_index(transcriptome_fasta,
                    params.minimap_index_extra_opts)
    ch_versions = ch_versions.mix(minimap2_index.out.versions)

    minimap2_transcriptome(full_length_reads,
                            minimap2_index.out.index,
                            params.minimap_extra_opts)
    ch_versions = ch_versions.mix(minimap2_transcriptome.out.versions)

    process_alignment_transcriptome(minimap2_transcriptome.out.sam)                    
    ch_versions = ch_versions.mix(process_alignment_transcriptome.out.versions)

    // Run salmon quant
    salmon(process_alignment_transcriptome.out.bam,
            transcriptome_fasta,
            params.salmon_extra_opts)
    ch_versions = ch_versions.mix(salmon.out.versions)
    ch_versions.view()

    // Combine all version information
    //versions (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    emit:
    salmon_quant = salmon.out.quant
    //versions = versions.out
}     
