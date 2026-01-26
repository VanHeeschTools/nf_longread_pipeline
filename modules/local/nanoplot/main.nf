// Create longread statistics, figures and reports using NanoPlot
process nanoplot {
    tag "$sample"
    label 'process_low'

    input:
        tuple val(sample), path(reads)
        val extra_opts

    output:
        path "${sample}_nanoplot", emit: nanoplot_dir
        path "versions.yml", emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        // Identify format of input files
        def input_format = reads[0].getExtension() == "bam" ? "--bam" : "--fastq"
        """
        # Create temp directory for nanoplot
        temp_dir="nanoplot_tmp"
        mkdir -p "\$temp_dir"
        export TMPDIR="\$(realpath "\$temp_dir")"

        NanoPlot $input_format $reads \
            -o ${sample}_nanoplot \
            -p ${sample}_ \
            -t 1 \
            --title "$sample" \
            --N50 \
            --raw \
            --tsv_stats \
            ${extra_opts}

        echo "NanoPlot completed for sample: $sample"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            NanoPlot: \$(NanoPlot -v)
        END_VERSIONS
        """
}