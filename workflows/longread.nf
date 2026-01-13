include { QC } from '../subworkflows/QC.nf'
include { ASSEMBLY } from '../subworkflows/ASSEMBLY.nf'
include { EXPRESSION } from '../subworkflows/EXPRESSION.nf'
include { VERSIONS } from '../modules/local/versions/main'
include { MULTIQC } from '../modules/local/multiqc/main'


// Function to validate samplesheet inputs (can be moved to a separate module)
def validateSampleSheet(sample_sheet) {
    return sample_sheet
        .splitCsv(header:true, sep:',')
        .map { row -> 
            if (!row.barcode && !row.sample) {
                error "Invalid sample sheet entry: ${row}. Either 'barcode' or 'sample' columns are required."
            }
            [row.barcode ?: row.sample, row.sample ?: row.barcode]
        }
}


def checkInputVar(input_dir) {
    if (file(input_dir).isDirectory()) {
        log.info "Processing directory: ${params.input}"
        def input_ch2 = channel
                .fromPath("${params.input}/**/*.{fastq,fastq.gz,bam}")
                .ifEmpty { error "No input files found in directory: ${params.input}" }
                .map { file -> 
                    def parent = file.parent.name
                    def barcode = params.barcode ?: (parent == params.input.tokenize('/')[-1] ? 'barcode' : parent)
                    //return tuple(sample, file)
                    [[id:"${barcode}"], [barcode:"${barcode}"], file]
                }
    } else {
        error "Input not specified. Please provide --input parameter."
    }
}


workflow LONGREAD {
    main:
    // Sample sheet handling with error check
    if (params.sample_sheet) {
        if (!file(params.sample_sheet).exists()) {
            error "Sample sheet file does not exist: ${params.sample_sheet}"
        }
        log.info "Using sample sheet: ${params.sample_sheet}"
        sample_sheet_ch = channel.fromPath(params.sample_sheet)
    } else {
        log.warn "No sample sheet provided, continuing without it."
        sample_sheet_ch = channel.empty()
    }


    // Define inputs from params
    // If sample sheet is provided, use it to update sample names
    if (sample_sheet_ch) {
        sample_map = validateSampleSheet(sample_sheet_ch)
 
        //Exit if sample_map is empty
        if (!sample_map) {
            log.error("ERROR: sample_map is null or empty! Check your sample sheet.")
            System.exit(1)
        }
    }

    
    reference = file(params.reference_genome, checkIfExists: true)
    annotation = file(params.reference_gtf, checkIfExists: true)

    ch_test = checkInputVar(params.input)
    sample_map = sample_map.map{ sample, barcode -> tuple( [id:"${sample}"], [sample:"${sample}"], [barcode:"${barcode}"])}
    ch_test_comb = sample_map.combine(ch_test.groupTuple(by: [0,1]), by: [0])
    ch_input2 = ch_test_comb.map{meta, sample, barcode, barcode2, file -> tuple(sample.sample, file)}



    if (params.qc) {
        if (params.direct_rna) {
            // Skip PyChopper, treat input as full_length_reads
            QC(ch_input2, params.direct_rna)
        } else {
            QC(ch_input2, params.direct_rna)
        }
        // Collect logs for MultiQC
        nanoplot_logs = QC.out.nanoplot_logs.collect()
        pychopper_logs = QC.out.pychopper_logs.collect()
        //  Collect full_length_reads for downstream steps
        full_length_reads = QC.out.full_length_reads
    } else {
        // If QC is skipped, set empty channels for logs
        nanoplot_logs = channel.empty()
        pychopper_logs = channel.empty()

        if (params.direct_rna) {
            // Use provided reads as full_length_reads
            full_length_reads = ch_input2
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
        ASSEMBLY(full_length_reads, reference, annotation)
        transcriptome_fasta = ASSEMBLY.out.fasta
        mapping_logs = ASSEMBLY.out.mapping_logs.collect()
        gffcompare_logs = ASSEMBLY.out.gffcompare_logs.collect()
        mapping_logs = channel.empty()
        gffcompare_logs = channel.empty()
    } else {
        mapping_logs = channel.empty()
        gffcompare_logs = channel.empty()
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
        EXPRESSION(full_length_reads, transcriptome_fasta)
    } else {
        log.warn "Expression analysis skipped."
        salmon_logs = channel.empty()
    }

    // TODO: Collect all versions.yml files
    //ch_versions = Channel.empty()
    //ch_versions = ch_versions.mix(MINIMAP2.out.versions)
    //ch_versions = ch_versions.mix(PROCESS_ALIGNMENT.out.versions)

    // Run the VERSIONS process
    //VERSIONS(ch_versions.collect())

    // Collect all output for MultiQC
    multiqc_files = channel.empty()
    if (params.qc) multiqc_files = multiqc_files.mix(nanoplot_logs)
    if (params.qc) multiqc_files = multiqc_files.mix(pychopper_logs)
    if (params.assembly) multiqc_files = multiqc_files.mix(mapping_logs)
    if (params.assembly) multiqc_files = multiqc_files.mix(gffcompare_logs)

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
    
    MULTIQC(multiqc_branch.run, file(params.multiqc_config))

    // For the skip branch, emit an empty channel
    multiqc_branch.skip
        .map { [] }
        .set { MULTIQC_REPORT }

    // Merge the MultiQC outputs
    MULTIQC_REPORT = MULTIQC_REPORT.mix(MULTIQC.out.report)

    emit:
    full_length_reads = params.qc ? QC.out.full_length_reads : channel.empty()
    merged_gtf = params.assembly ? ASSEMBLY.out.transcriptome : channel.empty()
    expression = params.expression ? EXPRESSION.out.salmon_quant : channel.empty()
    multiqc_report = MULTIQC_REPORT 
}