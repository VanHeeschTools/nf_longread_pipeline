include { jaffal } from '../modules/local/jaffal/main'

workflow FUSIONS {
    take:
        full_length_reads  // Path, location of input files
        jaffal_data_dir    // Path, location of directory with required annotation files
        genome_version     // String, genome version, JAFFAL automatically sets hg38 unless genome parameter is given
        annotation_version // String, annotation version, JAFFAL automatically sets genCode22 unless annotation parameter is given

    main:
        // full_length_reads emits: tuple(sample_id, file)
        full_length_fastq_reads = full_length_reads
            .map { _sample_id, file -> file }.collect()

        jaffal(full_length_fastq_reads,
            jaffal_data_dir,
            genome_version,
            annotation_version)
        
        jaffal_results_csv = jaffal.out.jaffa_results_csv
        jaffa_results_fasta = jaffal.out.jaffa_results_fasta

    emit:
        jaffal_results_csv
        jaffa_results_fasta

}
