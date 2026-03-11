// Run fusion detection with JAFFAL
process jaffal {
    label 'process_extreme'

    input:
        path full_length_reads // Path, location of input files
        val jaffal_data_dir    // Path, location of directory with required annotation files
        val genome_version     // String, genome version
        val annotation_version // String, annotation version

    output:
        path "*" // Remove after testing
        path "jaffa_results.csv", emit: jaffa_results_csv
        path "jaffa_results.fasta", emit: jaffa_results_fasta
        path "jaffal_mqc.csv", emit: jaffa_mqc
        path "versions.yml", emit:versions

    when:
        task.ext.when == null || task.ext.when 

    script:
        """
        # Run Jaffal
        bpipe run \
            -p inputPathsAreAbsolute=true \
            -p genome=${genome_version} \
            -p annotation=${annotation_version} \
            -p refBase=${jaffal_data_dir} \
            /JAFFA/JAFFAL.groovy \
            ${full_length_reads.join(' ')}

        # Write output statistics to MultiQC ready file
        awk -F',' 'NR==1 { next }
        {
            sample = \$1
            gsub(/_full_length_reads\\.fastq\$/, "", sample)

            class = \$17

            if (class == "HighConfidence") hc[sample]++
            else if (class == "LowConfidence") lc[sample]++
            else if (class == "PotentialTransSplicing") pts[sample]++
        }
        END {
            print "Sample,HighConfidence,LowConfidence,PotentialTransSplicing"
            for (s in hc) {
                print s "," hc[s]+0 "," lc[s]+0 "," pts[s]+0
            }
        }' jaffa_results.csv > jaffal_mqc.csv


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            JAFFA: 2.4
        END_VERSIONS
        """
}

process ctat_lr_fusion {
    label 'process_extreme'

    input:
        tuple val(sample_id), path(full_length_reads)  // Path, location of input files
        val ctat_lr_data_dir                           // Path, location of directory with required annotation files

    output:
        path "*" // Remove after testing

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        ctat-LR-fusion \
        -T ${full_length_reads} \
        --genome_lib_dir ${ctat_lr_data_dir} \
        --CPU ${task.cpus} \
        --examine_coding_effect \
        --extract_fusion_LR_fasta ${sample_id}_fusion_evidence_reads.fa \
        --vis \
        --output ${sample_id}_ctat_LR_fusion_outdir
        """
}