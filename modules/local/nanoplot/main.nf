process nanoplot {
    tag "$sample"
    label 'process_low'

    input:
    tuple val(sample), path(reads)
    val extra_opts

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
        --N50 \
        --raw \
        --tsv_stats \
        ${extra_opts}

    echo "NanoPlot completed for sample: $sample"
    """
}