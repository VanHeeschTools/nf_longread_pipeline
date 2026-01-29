// Run MultiQC using given directories and multiqc yaml file
process multiqc {
    label 'process_low'
    
    input:
        path multiqc_files, stageAs: "?/*" // Path, multiqc input files staged in individual folders
        path config                        // Path, multiqc config file found in /assets
        path software_versions_mqc         // Path, file containing the versions of tools used in the pipeline

    output:
        path "multiqc_report.html", emit: report

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        multiqc . -v -f -c $config
        """
}

