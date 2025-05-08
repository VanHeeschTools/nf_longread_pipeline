include { MINIMAP2_INDEX; MINIMAP2_TRANSCRIPTOME } from '../modules/local/minimap2/main'
include { PROCESS_ALIGNMENT_TRANSCRIPTOME } from '../modules/local/process_alignment/main'
include { SALMON } from '../modules/local/salmon/main'
include { VERSIONS } from '../modules/local/versions/main'

workflow EXPRESSION {
    take:
    full_length_reads
    transcriptome_fasta   

    main:
    // Create empty channel for versions
    ch_versions = Channel.empty()
    
     // Map against transcriptome
    MINIMAP2_INDEX(transcriptome_fasta,
                    params.minimap_index_extra_opts)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    MINIMAP2_TRANSCRIPTOME(full_length_reads,
                            MINIMAP2_INDEX.out.index,
                            params.minimap_extra_opts)
    ch_versions = ch_versions.mix(MINIMAP2_TRANSCRIPTOME.out.versions)

    PROCESS_ALIGNMENT_TRANSCRIPTOME(MINIMAP2_TRANSCRIPTOME.out.sam)                    
    ch_versions = ch_versions.mix(PROCESS_ALIGNMENT_TRANSCRIPTOME.out.versions)

    // Run Salmon quant
    SALMON(PROCESS_ALIGNMENT_TRANSCRIPTOME.out.bam,
            transcriptome_fasta,
            params.salmon_extra_opts)
    ch_versions = ch_versions.mix(SALMON.out.versions)
    ch_versions.view()

    // Combine all version information
    //VERSIONS (
    //    ch_versions.unique().collectFile(name: 'collated_versions.yml')
    //)

    emit:
    salmon_quant = SALMON.out.quant
    //versions = VERSIONS.out
}     
