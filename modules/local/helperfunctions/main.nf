def printHeader () {
    def logMessage =  """
        .-------------------------------------------------------.
        |                _____                 _      _     _   |
        | _ _ ___ ___   |  |  |___ ___ ___ ___| |_   | |___| |_ |
        || | | .'|   |  |     | -_| -_|_ -|  _|   |  | | .'| . ||
        | \\_/|__,|_|_|  |__|__|___|___|___|___|_|_|  |_|__,|___||
        '-------------------------------------------------------'

        ${workflow.manifest.name} ${workflow.manifest.version}
        ==========================
        """
        log.info logMessage.stripIndent()
}

// Copy an input samplesheet to outdir/samplesheet/
def copy_samplesheet(String input, String outdir) {
    // Check if input samplesheet and output directory paths are given
    if (!input || !outdir)
        return false

    // Check if input samplesheet file exists
    def inputFile = file(input)
    if (!inputFile.exists())
        return false

    // Creates the directory if needed
    def destDir = file("${outdir}/samplesheet")
    destDir.mkdirs()

    // Copy the samplesheet to output directory
    inputFile.copyTo(destDir.resolve(inputFile.name))

    return true
}

// Read samplesheet and emit sample IDs in file order
def readSampleIds(samplesheet) {
    Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample && !row.barcode) {
                error "Invalid samplesheet row: ${row}. 'sample' or 'barcode' is required."
            }
            row.sample ?: row.barcode
        }
}

// Function to validate samplesheet inputs
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

// Function to obtain all samples in given data directory linked to all given barcodes
def readInputDirectory(input_dir) {

    return Channel
        .fromPath("${input_dir}/**/*.{fastq,fastq.gz,bam}")
        .ifEmpty { error "No input files found in directory: ${input_dir}" }
        .map { file ->

            // Barcode is the parent folder name
            def barcode = file.parent.name

            tuple(barcode, file)
        }
}

// Function to create input channel for the workflow using given samplesheet and data directory
def buildSampleFileChannel(sample_sheet_ch, input_dir) {

    def sample_map_ch  = validateSampleSheet(sample_sheet_ch)

    // Create file channel from input directory
    def files_by_barcode = readInputDirectory(input_dir).groupTuple(by: 0)

    return sample_map_ch
        .combine(files_by_barcode, by: [0])
        .map { _barcode, sample, files ->
            tuple(sample, files)
        }
}