process multiqc {
    label 'process_low'
    
    input:
        path '*'
        path config

    output:
        path "multiqc_report.html", emit: report

    script:
        """
        multiqc . -f -c $config
        """
}

