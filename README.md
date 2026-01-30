# nf_longread_pipeline

A modular and reproducible Nextflow pipeline designed for processing long-read sequencing data from Oxford Nanopore. The pipeline encompasses quality control, mapping, assembly, expression quantification, and fusion gene detection.

## Installation

### Prerequisites

- **Nextflow** (version 23.04.1 or later) - [Install guide](https://www.nextflow.io/docs/latest/getstarted.html)
- **Java** (version 17 or later, up to 24) - [Download](https://www.oracle.com/java/technologies/downloads/)
- **Singularity** (for containerized execution) - [Install guide](https://docs.sylabs.io/guides/latest/user-guide/index.html)

### Installation Steps

```bash
# Verify Nextflow installation
nextflow -version

# (Optional) Set Nextflow home directory for storing downloaded pipelines
export NXF_HOME=~/.nextflow
mkdir -p $NXF_HOME

# (Optional) Create a Singularity cache directory for faster subsequent runs
export NXF_SINGULARITY_CACHEDIR=~/.singularity/cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Pull the pipeline 
nextflow pull VanHeeschTools/nf_longread_pipeline -r main
```

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

```bash
nextflow run nf_longread_pipeline/main.nf \
    -c /path/to/custom.config \
    -profile slurm
```

See [test/documentation/test.config](test/documentation/test.config) for an example configuration file.

**Resource allocation with SLURM profile:**
- Base process: 2 CPUs, 10 GB, 4 hours
- Low complexity: 8 CPUs, 12 GB, 8 hours
- Medium complexity: 12 CPUs, 36 GB, 16 hours
- High complexity: 16 CPUs, 128 GB, 48 hours
- Extreme complexity: 20 CPUs, 250 GB, 48 hours

Resources automatically scale with task retries (multiplied by attempt number).

## Pipeline Workflow

The pipeline executes the following main steps:

```
Input FASTQ files
    ↓
[QC] Quality Control
    ├─→ NanoPlot: Read length & quality visualization
    └─→ Pychopper: Adapter trimming & orientation correction
    ↓
[ASSEMBLY] Transcriptome Assembly
    ├─→ Minimap2: Alignment to reference genome
    ├─→ StringTie: Per-sample transcript assembly
    └─→ Gffcompare: Merge & compare transcripts across samples
    ↓
[EXPRESSION] Expression Quantification
    └─→ Salmon: Transcript abundance estimation (TPM)
    ↓
[FUSIONS] Fusion Detection (Optional)
    └─→ JAFFAL: Fusion gene prediction
    ↓
Quality Report
    └─→ MultiQC: Aggregate all QC metrics
```

### Conditional Execution

Each major workflow step can be toggled independently:

- **`--qc true/false`** (default: true) - Skip NanoPlot and Pychopper if reads are pre-processed
- **`--assembly true/false`** (default: true) - Skip assembly if using external transcriptome
- **`--expression true/false`** (default: true) - Skip quantification
- **`--fusions true/false`** (default: true) - Skip fusion detection

## Input requirements

### Input/output files


| Parameter      | Description                                                                                     |
| :--------------- | ------------------------------------------------------------------------------------------------- |
| --input        | **(Required)** Path to directory with raw FASTQ files (accepts subdirectories one level deep) |
| --sample_sheet | **(Required)** CSV file with sample metadata (see format below)                                 |
| --outdir       | **(Optional)** Path to output directory, defaults to `results`                                  |

**Input file format notes:**
- FASTQ files can be gzipped or uncompressed (.fastq, .fastq.gz, .fq, .fq.gz)
- Input directory structure: `/path/to/fastq/` containing fastq files or single-level subdirectories with fastq files
- Example: `fastq/sample1.fastq.gz`, `fastq/subdir/sample2.fastq`, `fastq/sample3.fq.gz`

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
| --salmon_extra_opts | Salmon-specific flags (e.g. --seqBias) | "--ont" |

#### Fusions (JAFFAL)


| Parameter                 | Description                                                                             | Default |
| :-------------------------- | --------------------------------------------------------------------------------------- | --------- |
| --fusions                 | Enable fusion gene detection with JAFFAL                                               | true    |
| --genome_version          | Genome assembly version used in JAFFAL database (hg38, hg19, mm39, mm10, etc.)         | hg38    |
| --annotation_version      | Annotation version in JAFFAL database (e.g., genCode48, genCode45, M35)                | genCode48|

## Sample sheet format

Samplesheet is a CSV file with the following columns:


| Column  | Description |
| :-------- | ------------- |
| barcode | Barcode id (leave empty for singleplex) |
| sample  | Sample name (used in output filenames) |

```csv
barcode,sample
BC1,sample1
BC2,sample2
BC3,sample3
```

The barcode field can be left empty for singleplex experiments:

```csv
barcode,sample
,sample1
,sample2
,sample3
```

## Filter & Annotation Parameters

These parameters control transcript filtering and annotation in the merged GTF output:

### Transcript Filtering

Reconstructed transcripts are filtered based on expression levels and sample presence:

| Parameter        | Description                                                           | Default |
| :----------------- | ----------------------------------------------------------------------- | --------- |
| --min_tpm        | Minimum TPM (Transcripts Per Million) to retain transcripts           | 0.1     |
| --min_occurrence | Minimum number of samples a transcript must appear in to be retained  | 1       |

**How filtering works:**
1. Only transcripts with TPM ≥ `--min_tpm` in at least one sample are kept
2. Transcripts must appear in ≥ `--min_occurrence` samples
3. Reference transcripts (from annotation GTF) are always included regardless of TPM, annotated according to wether they have evidence from long-read transcriptome reconstruction: 
- *LR_full_match*:   reference transcripts fully matching long-read models (=)
- *LR_partial_match*:  reference transcripts with partially overlapping long-read models (c,j)
- *LR_other_match*:  reference transcripts with antisense and other long-read overlaps 
- *LR_no_evidence*: reference transcripts without any overlapping long-read model


## Output

### Output directory structure

```
results/
├── nanoplot/
├── pychopper/
├── minimap2/
├── stringtie/
...
├── merged_gtf/
├── salmon/
├── jaffal/
└── multiqc/
```

### Output files

#### nanoplot/

Quality control visualizations for raw sequencing reads per sample. Each sample has a subdirectory containing:

| File | Description |
| :---- | ------------- |
| `*_NanoPlot-report.html` | Interactive HTML report with read quality metrics |
| `*_NanoStats.txt` | Summary statistics (read count, length, quality) |
| `*_NanoPlot-data.tsv.gz` | Underlying data for all plots |
| `*_HistogramReadlength.png`/`.html` | Read length distribution plots |
| `*_LengthvsQualityScatterPlot_*.png`/`.html` | Scatter plots of read length vs quality score |
| `*_Yield_By_Length.png`/`.html` | Cumulative yield by read length |

#### pychopper/

Orientation and quality filtering results from cDNA library prep validation. Outputs per sample include:

| File | Description |
| :---- | ------------- |
| `*_pychopper_report.pdf` | Quality report with adapter detection metrics |
| `*_pychopper_stats.tsv` | Statistics on full-length vs partial reads |
| `full_length_reads/` | Directory containing FASTQ files of full-length oriented reads |

#### minimap2/

Read alignment files mapping sequencing reads to reference genome:

| File | Description |
| :---- | ------------- |
| `*.bam` | Aligned reads in BAM format (if minimap2 ran) |
| `*.bam.bai` | BAM index files |

#### stringtie/

Transcriptome assembly output per sample:

| File | Description |
| :---- | ------------- |
| `*.gff` | Transcript annotations in GFF3 format with expression values (TPM) |
| `*_stringtie.log` | StringTie execution log |

#### merged_gtf/

Consolidated transcript annotation from all samples:

| File | Description |
| :---- | ------------- |
| `output.filtered.extended_reference.gtf` | Extended reference annotation with novel transcripts after filtering |
| `output.filtered.novel_transcripts.gtf` | Novel transcripts not in reference after filtering |
| `output.filtered.tsv` | Novel transcript presence across samples (stringtie TPM) |
| `output.filtered.log` | Logs filtering and annotation of reconstructed transcripts  |
| `output.combined.gtf` | Merged GTF file combining reference and novel transcripts |
| `output.tracking` | Gffcompare tracking file linking transcripts to reference |
| `output.stats` | Gffcompare comparison statistics |
| `output_transcript_presence.tsv` | Binary table of transcript presence/absence across samples |

#### salmon/

Quantification results for transcript expression levels:

| File | Description |
| :---- | ------------- |
| `*/quant.sf` | Per-sample salmon quantification with TPM and count values |

#### salmon_tables/

Aggregated expressiond data for immediate use at transcript and gene level:

| `*_transcript_counts.tsv` | Raw read counts aggregated at transcript level (samples × transcripts) |
| `*_transcript_tpms.tsv` | Transcript-level TPM (Transcripts Per Million) values |
| `*_gene_ID_counts.tsv` | Raw read counts aggregated by gene ID (samples × genes) |
| `*_gene_ID_tpms.tsv` | Gene-level TPM values aggregated by gene ID |
| `*_gene_name_counts.tsv` | Raw read counts aggregated by gene name (samples × genes) |
| `*_gene_name_tpms.tsv` | Gene-level TPM values aggregated by gene name |
| `*_tx2gene.tsv` | Mapping table linking transcript IDs to gene IDs and gene names |
| `*_multiqc_summary_mqc.tsv` | Summary statistics table for MultiQC including read mapping rates, expressed transcript counts, and mean/median TPM values |
| `*_log.txt` | Log file with warnings/errors from salmon_tables generation |

#### jaffal/

Fusion gene detection results:

| File | Description |
| :---- | ------------- |
| `jaffa_results.csv` | Predicted fusion events with supporting evidence |
| `jaffa_results.fasta` | Sequences of predicted fusion junctions |
| `*_full_length_reads.fastq` | Full-length reads supporting fusion predictions |
| `jaffal_mqc.csv` | Summary statistics for MultiQC integration |

#### seqkit_stats/

Basic sequence statistics for alignment files:

| File | Description |
| :---- | ------------- |
| `minimap2_bams_seqkit_stats.tsv` | Read count and length statistics from aligned BAM files |

#### samplesheet/

Input sample metadata:

| File | Description |
| :---- | ------------- |
| `samplesheet.csv` | Copy of input sample sheet with barcode and sample information |

#### output_samplesheet/

Processed sample sheet for downstream analysis:

| File | Description |
| :---- | ------------- |
| `output_samplesheet.csv` | Sample sheet of output files |

#### multiqc/

Aggregated quality control report:

| File | Description |
| :---- | ------------- |
| `multiqc_report.html` | Interactive HTML report combining metrics from all analysis steps |
| `multiqc_data/` | Directory with underlying data and plots for the multiqc report |

## Compute Resource Requirements

### Typical Resource Usage

Resource requirements vary depending on:
- Number of samples
- Read depth (number of reads per sample)
- Genome size
- Selected modules (QC, assembly, expression, fusions)

### Estimated Requirements by Dataset Size

| Dataset Type | Samples | CPUs | Memory | Storage | Time |
| :-------------- | --------- | ------ | -------- | --------- | ------ |
| Small test | 1-3 | 8 | 16 GB | 5 GB | 0.5-1 h |
| Medium | 5-10 | 12 | 36 GB | 50 GB | 2-4 h |
| Large | 20+ | 16-32 | 128 GB | 200+ GB | 8-24 h |

**Storage notes:**
- Input FASTQ files are NOT deleted automatically
- Work directory (`./work/`) can grow large; clean up after successful runs with `nextflow clean`
- Output directory contains final results and is configurable via `--outdir`

## Troubleshooting

### Common Issues & Solutions

#### 1. Sample sheet validation error

**Error message:** `ERROR: Sample sheet file does not exist`

**Solution:**
```bash
# Ensure the sample sheet path is absolute and file exists
ls -la /path/to/sample_sheet.csv

# Check sample sheet format (must have barcode,sample columns)
head /path/to/sample_sheet.csv
```

#### 2. Reference genome file not found

**Error message:** `ERROR: Cannot find file: /path/to/reference_genome.fasta`

**Solution:**
- Verify file path is absolute (not relative)
- Check file permissions: `ls -l /path/to/reference_genome.fasta`
- Ensure file is not corrupted: `zcat /path/to/reference_genome.fasta | head -2`

#### 3. Full-length reads not found (when --qc false)

**Error message:** `Full length reads not found in results/pychopper/full_length_reads`

**Solution:**
- Either re-run with `--qc true` (default) to generate full-length reads
- Or provide `--direct_rna true` if reads are already oriented and full-length
- Or ensure pychopper has been run previously and outputs are in the correct location

#### 4. Singularity container image download fails

**Error message:** `Failed to pull Singularity image`

**Solution:**
```bash
# Set up Singularity cache directory
export NXF_SINGULARITY_CACHEDIR=~/.singularity/cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# Pre-download container images (optional)
# Re-run the pipeline; containers will be cached for future runs
```

#### 5. SLURM job fails with memory error

**Error message:** `slurmstepd: error: Exceeded job memory limit`

**Solution:**
- Check your dataset size and adjust `max_memory` in base.config
- Or submit with higher time limits: `--max_memory 256.GB`
- Reduce parallelization if on a memory-constrained system

#### 6. No output files produced

**Diagnosis steps:**
```bash
# Check execution report
tail -50 .nextflow.log

# Examine work directory for error messages
find work/ -name ".command.err" -exec grep -l "Error" {} \;

# Run with verbose output
nextflow run main.nf ... -resume -v
```

### Getting Help

- **Check logs:** Review `.nextflow.log` and `execution_trace.txt` in output directory
- **Resume interrupted runs:** Use `-resume` flag to restart from last successful process
- **Dry run:** Use `-n` flag to preview workflow without execution
- **Report issues:** Submit issues to the [GitHub repository](https://github.com/VanHeeschTools/nf_longread_pipeline/issues)

## Citation & Acknowledgments

**Authors:** Marina Reixachs Sole, Edwin van der Werf, Rico Hagelaar
