process multiqc {
    label 'process_low'
    
    input:
        path '*'
        path config

    output:
        path "multiqc_report.html", emit: report

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        multiqc . -f -c $config
        """
}

