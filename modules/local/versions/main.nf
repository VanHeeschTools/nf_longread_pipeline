process VERSIONS {
    input:
    path versions

    output:
    path "software_versions.yml"
    path "software_versions_mqc.yml"

    script:
    """
    # Combine all version files into a single file
    cat $versions > software_versions.yml

    # Create a simplified version for MultiQC
    echo "---" > software_versions_mqc.yml
    cat software_versions.yml | grep -v "^{" >> software_versions_mqc.yml
    """
}