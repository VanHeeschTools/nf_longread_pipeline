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

    return channel
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

    /*def sample_map_ch  = validateSampleSheet(sample_sheet_ch)

    // Create file channel from input directory
    def files_by_barcode = readInputDirectory(input_dir).groupTuple(by: 0)
    files_by_barcode = files_by_barcode.map{ file -> [ file.simpleName.split('_')[0], file ] }.groupTuple()

    

    return sample_map_ch
        .join(files_by_barcode) 
        .map { barcode, sample, files -> 
            // .flatten() ensures the list isn't [[f1, f2], [f3, f4]]
            return [ sample, files.flatten() ] 
        }

    */
    def sample_map_ch = validateSampleSheet(sample_sheet_ch)

    // 1. Get the files. readInputDirectory returns [barcode, file]
    // 2. groupTuple(by: 0) collects all files sharing the same barcode into a list
    def files_by_barcode = readInputDirectory(input_dir).groupTuple(by: 0)

    sample_map_ch.join(files_by_barcode).view { it -> "MATCHED: Barcode ${it[0]} -> Sample ${it[1]}" }

    return sample_map_ch
        .join(files_by_barcode)
        .map { barcode, sample, files -> 
            // Use .flatten() to ensure you don't have [[f1], [f2]]
            // but rather [f1, f2]
            return tuple(sample, files.flatten()) 
        }
}

