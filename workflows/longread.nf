include { QC }         from '../subworkflows/QC.nf'
include { ASSEMBLY }   from '../subworkflows/ASSEMBLY.nf'
include { EXPRESSION } from '../subworkflows/EXPRESSION.nf'
include { FUSIONS }    from '../subworkflows/FUSIONS.nf'
include { versions }   from '../modules/local/versions/main'
include { multiqc }    from '../modules/local/multiqc/main'
include { buildSampleFileChannel; copy_samplesheet } from '../modules/local/helperfunctions/main.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow LONGREAD {
    main:

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Sample sheet handling with error check
    if (params.sample_sheet) {
        if (!file(params.sample_sheet).exists()) {
            error "ERROR: Sample sheet file does not exist: ${params.sample_sheet}"
            System.exit(1)
        }
        log.info "Using sample sheet: ${params.sample_sheet}"
        def sample_sheet_ch = Channel.fromPath(params.sample_sheet)
        
         // Read samplesheet and create sampe input channel
        input_data = buildSampleFileChannel(sample_sheet_ch, params.input)
        copy_samplesheet(params.sample_sheet, params.outdir)

    } else {
        log.error("ERROR: params.sample_sheet is null or empty! Please set this parameter.")
        System.exit(1)
    }

    // Load required files
    reference_genome = [file(params.reference_genome, checkIfExists: true), file("${params.reference_genome}.fai", checkIfExists: true)]
    annotation = file(params.reference_gtf, checkIfExists: true)

    // Declare empty channels
    nanoplot_logs = channel.empty()
    pychopper_logs = channel.empty()
    mapping_logs = channel.empty()

    if (params.qc) {
        if (params.direct_rna) {
            // Skip PyChopper, treat input as full_length_reads
            QC(input_data, params.direct_rna)
        } else {
            QC(input_data, params.direct_rna)
        }

        // Collect logs for MultiQC
        nanoplot_logs = QC.out.nanoplot_logs.collect()
        pychopper_logs = QC.out.pychopper_logs.collect()
        nanoplot_html = QC.out.nanoplot_html.collect()
        
        // Collect full_length_reads for downstream steps
        full_length_reads = QC.out.full_length_reads
    } else {

        if (params.direct_rna) {
            // Use provided reads as full_length_reads
            full_length_reads = input_data
        } else {
            // If not running QC and not direct RNA, assume pychopper has been run externally
            // and full_length_reads are in the expected directory
            // Use Channel.fromPath to collect files
            full_length_reads = channel.fromPath("${params.outdir}/pychopper/full_length_reads/*.{fastq,fq,fastq.gz,fq.gz}")
                // Check if channel is empty and provide error message
                .ifEmpty {
                    error "Full length reads not found in ${params.outdir}/pychopper/full_length_reads. Please run QC step, provide full length reads, or set --direct-rna to skip pychopper."
                }
                // Mimic tuple(sample, file) structure from input_ch
                // Replace _full_length_reads and extensions from filename to get sample name
                .map { file ->
                    def sample = file.getBaseName().replaceFirst(/_full_length_reads(\.fastq\.gz|\.fq\.gz|\.fastq|\.fq)?$/, '')
                    tuple(sample, file)
                }
        }
    }

    if (params.assembly) {
        ASSEMBLY(full_length_reads, reference_genome, annotation)
        stringtie_mqc = ASSEMBLY.out.stringtie_mqc
        transcriptome_fasta = ASSEMBLY.out.transcriptome_fasta
        mapping_logs = ASSEMBLY.out.mapping_logs.collect()
    } else {
        log.warn "Assembly step skipped."

        //Assign transcriptome fasta if expression is true
        if (params.expression) {
            if (params.transcriptome_fasta) {
                transcriptome_fasta = file(params.transcriptome_fasta, checkIfExists: true)
                log.warn "Using provided transcriptome FASTA: ${params.transcriptome_fasta}."
            } else {
                error "Without assembly, a reference transcriptome fasta file must be provided as `params.transcriptome_fasta`."
            }
        }
    }

    if (params.expression) {    
        EXPRESSION(full_length_reads, transcriptome_fasta, annotation)
        salmon_multiqc = EXPRESSION.out.salmon_multiqc
    }


    if (params.fusions) {
        FUSIONS(full_length_reads,
            params.jaffal_data_dir,
            params.genome_version,
            params.annotation_version)
        
        jaffal_mqc = FUSIONS.out.jaffal_mqc
    }

    // Collect all tool versions
    ch_versions = Channel.empty()
    if (params.qc)  ch_versions = ch_versions.mix(QC.out.versions)
    if (params.assembly) ch_versions = ch_versions.mix(ASSEMBLY.out.versions)
    if (params.expression) ch_versions = ch_versions.mix(EXPRESSION.out.versions)
    if (params.fusions) ch_versions = ch_versions.mix(FUSIONS.out.versions)

    // Run the VERSIONS process
    versions(ch_versions.collect())
    software_versions_mqc = versions.out.software_versions_mqc

    // Collect all output for MultiQC
    multiqc_files = channel.empty()
    if (params.qc) multiqc_files = multiqc_files.mix(nanoplot_logs)
    if (params.qc) multiqc_files = multiqc_files.mix(pychopper_logs)
    if (params.qc) multiqc_files = multiqc_files.mix(nanoplot_html)
    if (params.assembly) multiqc_files = multiqc_files.mix(stringtie_mqc)
    if (params.assembly) multiqc_files = multiqc_files.mix(mapping_logs)
    if (params.expression) multiqc_files = multiqc_files.mix(salmon_multiqc)
    if (params.fusions) multiqc_files = multiqc_files.mix(jaffal_mqc)

    // Convert to list and check if empty
    multiqc_input = multiqc_files.collect().map { files -> 
        files.isEmpty() ? null : files 
    }

    // Run MultiQC only if there are input files
    MULTIQC_REPORT = channel.empty()
    multiqc_input.branch {
        run: it != null
        skip: it == null
    }.set { multiqc_branch }
    
    multiqc(multiqc_branch.run, file(params.multiqc_config), software_versions_mqc)

    // For the skip branch, emit an empty channel
    multiqc_branch.skip
        .map { [] }
        .set { MULTIQC_REPORT }

    // Merge the MultiQC outputs
    MULTIQC_REPORT = MULTIQC_REPORT.mix(multiqc.out.report)

    emit:
    full_length_reads = params.qc ? QC.out.full_length_reads : channel.empty()
    merged_gtf = params.assembly ? ASSEMBLY.out.transcriptome_gtf : channel.empty()
    expression = params.expression ? EXPRESSION.out.salmon_quant : channel.empty()
    multiqc_report = MULTIQC_REPORT 
}
