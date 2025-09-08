# nf_longread_pipeline

A modular and reproducible Nextflow pipeline designed for processing long-read sequencing data from Oxford Nanopore. The pipeline encompasses quality control, mapping, assembly, and quantification.

## Requirements

* Nextflow 23.04.1 or later
* Java 17 or later (up to 24)

### Containerised software

* NanoPlot
* Pychopper
* Minimap2
* Samtools
* Stringtie (3.0 recommended)
* Gffcompare
* Salmon
* R (tested on 4.1.2), including the following packages:
  * tidyverse
  * GenomicRanges
  * rtracklayer
  * data.table
  * dplyr
  * tidyr

Note: The pipeline is developed and tested with one singularity containers for each process.

## Usage

The pipeline can be run with a command as follows:

```bash
nextflow run nf_longread_pipeline/main.nf \
    --input /path/to/fastq/folder/ \
    --sample_sheet /path/to/sample_sheet.csv \
    --outdir /path/to/output/ \
    --reference_genome /path/to/reference_genome.fasta \
    --reference_gtf /path/to/reference_gtf.gtf \
    --reference_transcriptome /path/to/reference_transcriptome.fasta \
    -profile slurm
```

Or by specifying a config file (see [link] as example):

```bash
nextflow run nf_longread_pipeline/mains.nf \
    -c /path/to/local/file.config \
    -profile slurm
```

## Input requirements

### Input/output files


| Parameter      | Description                                                                                     |
| :--------------- | ------------------------------------------------------------------------------------------------- |
| --input        | **(Required)** Path to directory with raw FASTQ files (accepts a subdirectories one level deep) |
| --sample_sheet | **(Required)** CSV file with sample metadata (see format below)                                 |
| --outdir       | **(Optional)** Path to output directory, defaults to `results`                                  |

### Reference files


| Parameter                 | Description                                            |
| :-------------------------- | -------------------------------------------------------- |
| --reference_genome        | **(Required)** FASTA file with reference genome        |
| --reference_gtf           | **(Required)** GTF file with reference annotations     |
| --reference_transcriptome | **(Optional)** FASTA file with reference transcriptome |

### Module Toggles


| Parameter    | Description                                  | Default |
| :------------- | ---------------------------------------------- | --------- |
| --qc         | Enable quality control (NanoPlot, Pychopper) | true    |
| --assembly   | Enable transcriptome assembly                | true    |
| --expression | Enable quantification (Salmon)               | true    |

Note: QC false will assume pychopper has run previously and full length reads are in the expected location `pychopper/full_length_reads`. Provide `--direct_rna` option (see below) the samplesheet already contains the path to full length oriented reads.

### Optional Parameters

#### Nanoplot


| Parameter             | Description                                                                         | Default |
| :---------------------- | ------------------------------------------------------------------------------------- | --------- |
| --nanoplot_extra_opts | Additional opts. See [Nanoplot documentation](https://github.com/wdecoster/NanoPlot) | ""      |

#### Pychopper


| Parameter                 | Description                                                                             | Default |
| :-------------------------- | ----------------------------------------------------------------------------------------- | --------- |
| --direct_rna                | Boolean to treat input reads as oriented & full-length (Pychopper will not run).                         | false  |
| --cdna_kit                | See Pychopper docs for accepted kits                                                    | PCB114  |
| --custom_primers_file     | FASTA file with custom primers                                                          | null    |
| --pychopper_backend       | edlib or phmm                                                                           | edlib   |
| --pychopper_extra_opts    | Additional opts. See [Pychopper documentation](https://github.com/epi2me-labs/pychopper) | ""      |
| --store_full_length_reads | Save full-length reads after trimming                                                   | true    |

#### Minimap2


| Parameter                        | Description                                                                                      | Default |
| :--------------------------------- | -------------------------------------------------------------------------------------------------- | --------- |
| --minimap_extra_opts             | Minimap genome alignment opts (for assembly). See [Minimap docs](https://github.com/lh3/minimap2) | ""      |
| --minimap_index_extra_opts       | Minimap indexing opts.                                                                           | ""      |
| --minimap_transcripts_extra_opts | Minimap transcriptome alignment opts (for quantification)                                        | ""      |

#### StringTie


| Parameter              | Description                  | Default |
| :----------------------- | ------------------------------ | --------- |
| --stringtie_extra_opts | Additional StringTie options | ""      |

#### Filter Annotation


| Parameter        | Description                   | Default |
| :----------------- | ------------------------------- | --------- |
| --min_tpm        | Min TPM to retain transcripts | 0.1     |
| --min_occurrence | Min number of samples present | 1       |

#### Salmon


| Parameter           | Description                           | Default |
| :-------------------- | --------------------------------------- | --------- |
| --salmon_extra_opts | Salmon-specific flags (e.g. --gcBias) | "--ont" |

## Sample sheet format

Sampesheet is a CSV file with the following columns:


| Column  | Description |
| :-------- | ------------- |
| barcode | Barcode id  |
| sample  | Sample name |

```csv
barcode,sample
BC1,sample1
BC2,sample2
BC3,sample3
```

The barcode field can be left empty for singleplex experiments.

```csv
barcode,sample
,sample1
,sample2
,sample3
```

## Output structure

```
results/
├── nanoplot/
├── pychopper/
├── minimap2/
├── stringtie/
├── gffcompare/
├── salmon/
└── multiqc/
```
