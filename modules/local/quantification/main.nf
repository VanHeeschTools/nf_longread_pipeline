process QUANTIFICATION {
    conda "bioconda::stringtie=3.0.1"
    
    input:
    path transcripts

    output:
    path "gene_abundances.tsv"

    script:
    """
    stringtie $transcripts -e -B -p 8 -G $params.reference_gtf -o gene_abundances.tsv
    """
}