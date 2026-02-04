process versions {
    label 'process_superlow'

    input:
    path versions, stageAs: "?/*"

    output:
    path "software_versions.yml"
    path "software_mqc_versions.yml", emit: software_versions_mqc

    script:
    """
    # Combine all version files into a single file
    cat $versions > software_versions.yml

    # Simplify for MultiQC:
    awk '
    {
        # Skip lines without a space
        if (index(\$0, " ") == 0) next

        split(\$0, f, " ")
        tool = f[1]
        ver  = f[2]

        gsub(":", "", tool)

        if (tool == "" || ver == "") next
        # Skip already seen tools
        if (seen[tool]++) next

        print tool ": \\"" ver "\\""
    }
    ' software_versions.yml > software_mqc_versions.yml
    """
}

