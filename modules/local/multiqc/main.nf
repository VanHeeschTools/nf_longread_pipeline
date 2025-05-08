process MULTIQC {
    label 'process_low'
    
    input:
    path '*'
    path config

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    #Add projectDir for custom multiqc scripts
    export PYTHONPATH=$projectDir/bin:$PYTHONPATH

    multiqc . -f -c $config
    """
}

