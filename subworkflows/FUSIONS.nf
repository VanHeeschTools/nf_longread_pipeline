include { jaffal; ctat_lr_fusion } from '../modules/local/fusion_calling/main'

workflow FUSIONS {
    take:
        full_length_reads  // Path, location of input files
        jaffal_data_dir    // Path, location of directory with required annotation files for Jaffal
        genome_version     // String, genome version, JAFFAL automatically sets hg38 unless this parameter is given
        annotation_version // String, annotation version, JAFFAL automatically sets genCode22 unless this parameter is given
        ctat_lr_data_dir   // Path, location of directory with required annotation files for ctat_lr_fusion

    main:
        // Create empty channel for versions
        ch_versions = Channel.empty()

        // Only keep the fastq files of the full_length_reads tuple
        full_length_fastq_reads = full_length_reads
            .map { _sample_id, file -> file }.collect()

        // Run JAFFAL
        jaffal(full_length_fastq_reads,
            jaffal_data_dir,
            genome_version,
            annotation_version)
        ch_versions = ch_versions.mix(jaffal.out.versions)

        ctat_lr_fusion(full_length_reads, ctat_lr_data_dir)

        // Define JAFFAL output csv and fasta
        jaffal_results_csv = jaffal.out.jaffa_results_csv
        jaffal_results_fasta = jaffal.out.jaffa_results_fasta
        jaffal_mqc = jaffal.out.jaffa_mqc

    emit:
        jaffal_results_csv=jaffal_results_csv
        jaffal_results_fasta=jaffal_results_fasta
        jaffal_mqc=jaffal_mqc
        versions=ch_versions

}