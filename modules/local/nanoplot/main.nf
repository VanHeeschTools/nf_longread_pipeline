process nanoplot {
    tag "$sample"
    label 'process_low'

    input:
    tuple val(sample), path(reads)

    output:
    path "${sample}_nanoplot"

    script:
    //detect if input is single file or list of files
    def input_files = reads instanceof List ? reads.join(' ') : [reads] 
    //concatenate list of files
    def input_string = "${input_files.join(' ')}"
    //identify format of input files
    def input_format = reads[0].getExtension() == "bam" ? "--bam" : "--fastq_rich"
    """
    NanoPlot $input_format $reads \
        -o ${sample}_nanoplot \
        -p ${sample}_ \
        -t ${task.cpus} \
        --title "$sample" \
        --color darkseagreen \
        --minlength 50 \
        --N50 \
        --store

    echo "NanoPlot completed for sample: $sample"
    """
}