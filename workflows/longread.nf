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

workflow LONGREAD {
    take:
    input_ch
    sample_sheet_ch

    main:
    // Define inputs from params
    // If sample sheet is provided, use it to update sample names
    if (sample_sheet_ch) {

        sample_map = validateSampleSheet(sample_sheet_ch)

        //Exit if sample_map is empty
        if (!sample_map) {
            log.error("ERROR: sample_map is null or empty! Check your sample sheet.")
            System.exit(1)
        }

        input_ch = input_ch.join(sample_map, by: 0)
                    .map { barcode, file, sample -> tuple(sample, file) }
                    .ifEmpty { 
                        log.error("""
                        ERROR: input_ch is empty!
                        - Check the structure of your input directory.
                        - Parent directories must match the sample names in the sample sheet.
                        - Check your sample sheet for misspelled sample names.
                        """.stripIndent())
                        System.exit(1)
                    }
    }

    reference = file(params.reference_genome, checkIfExists: true)
    annotation = file(params.reference_gtf, checkIfExists: true)

    if (params.qc) {
        QC(input_ch)
        nanoplot_logs = QC.out.nanoplot_logs.collect()
        pychopper_logs = QC.out.pychopper_logs.collect()
    } else {
        log.warn "QC step skipped. Ensure full length reads are in the correct location: ${params.outdir}/pychopper/full_length_reads."
        nanoplot_logs = Channel.empty()
        pychopper_logs =  Channel.empty()

        // Check for full length reads in the expected output directory
        if (!file("${params.outdir}/pychopper/full_length_reads").exists()) {
            error "Full length reads directory not found. Please run QC step or provide full length reads reads."
        } else {
            full_length_reads = file("${params.outdir}/pychopper/full_length_reads")
        }
    }

    if (params.assembly) {
        ASSEMBLY(QC.out.full_length_reads, reference, annotation)
        transcriptome_fasta = ASSEMBLY.out.fasta
        mapping_logs = ASSEMBLY.out.mapping_logs.collect()
        gffcompare_logs = ASSEMBLY.out.gffcompare_logs.collect()
    } else {
        mapping_logs = Channel.empty()
        gffcompare_logs = Channel.empty()
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
        EXPRESSION(QC.out.full_length_reads, transcriptome_fasta)
    } else {
        log.warn "Expression analysis skipped."
        salmon_logs = Channel.empty()
    }

    // TODO: Collect all versions.yml files
    //ch_versions = Channel.empty()
    //ch_versions = ch_versions.mix(MINIMAP2.out.versions)
    //ch_versions = ch_versions.mix(PROCESS_ALIGNMENT.out.versions)

    // Run the VERSIONS process
    //VERSIONS(ch_versions.collect())

    // Collect all output for MultiQC
    multiqc_files = Channel.empty()
    if (params.qc) multiqc_files = multiqc_files.mix(nanoplot_logs)
    if (params.qc) multiqc_files = multiqc_files.mix(pychopper_logs)
    if (params.assembly) multiqc_files = multiqc_files.mix(mapping_logs)
    if (params.assembly) multiqc_files = multiqc_files.mix(gffcompare_logs)

    // Convert to list and check if empty
    multiqc_input = multiqc_files.collect().map { files -> 
        files.isEmpty() ? null : files 
    }

    // Run MultiQC only if there are input files
    MULTIQC_REPORT = Channel.empty()
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

    multiqc_files.view()

    emit:
    full_length_reads = params.qc ? QC.out.full_length_reads : Channel.empty()
    merged_gtf = params.assembly ? ASSEMBLY.out.transcriptome : Channel.empty()
    expression = params.expression ? EXPRESSION.out.salmon_quant : Channel.empty()
    multiqc_report = MULTIQC_REPORT
}