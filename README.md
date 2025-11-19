# Metagenomics KEGG Mapping and Analysis Pipeline

A comprehensive, automated pipeline for taxonomic profiling and functional annotation of metagenomic samples using SingleM and KEGG databases, designed for HPC SLURM clusters.

## Overview

This pipeline processes raw paired-end metagenomic FASTQ files through:

1. **SingleM Taxonomic Profiling** - Community composition analysis
2. **SingleM Treatment Group Summarization** - Combined analyses at all taxonomic levels
3. **Quality Filtering** - Trimmomatic preprocessing
4. **KEGG Functional Annotation** - Diamond BLASTX alignment
5. **Post-processing** - KEGG hit extraction, summarization, and normalization

### Pipeline Workflow

```
Raw FASTQ Files
    ↓
┌─────────────────────────────────────┐
│  STEP 1: SingleM Profiling          │
│  - Taxonomic classification         │
│  - Multi-level summarization        │
│  - Prokaryotic fraction estimation  │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│  STEP 2: Treatment Group Summary    │
│  - Combined summaries by group      │
│  - All taxonomic levels             │
│  (phylum → species)                 │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│  STEP 3: Diamond KEGG Alignment     │
│  - Trimmomatic QC filtering         │
│  - BLASTX against KEGG database     │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│  STEP 4: Extract KEGG Hits          │
│  - Parse Diamond output             │
│  - Map to KEGG orthology numbers    │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│  STEP 5: Sum KEGG Hits              │
│  - Aggregate hits per KO            │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│  STEP 6: Calculate Genome           │
│           Equivalents               │
│  - Mean single-copy gene abundance  │
└─────────────────────────────────────┘
    ↓
┌─────────────────────────────────────┐
│  STEP 7: Normalize Feature Table    │
│  - Genome-equivalent normalization  │
│  - Final feature matrix             │
└─────────────────────────────────────┘
    ↓
Final Results
```

## Features

- **Flexible FASTQ Format Support** - Handles `.fq.gz`, `.fastq.gz`, `_R1/_R2`, `_1/_2` naming conventions
- **Automated Job Dependency** - SLURM job chaining ensures correct execution order
- **Treatment Group Analysis** - Automated summarization across biological replicates
- **Multi-level Taxonomic Profiling** - From phylum to species level
- **Robust Error Handling** - Logging, file locking, and validation at each step
- **Resource Optimization** - Parallel processing with configurable resource allocation
- **Checkpointing** - Skips completed samples automatically

## Requirements

### Software

- **SLURM** workload manager
- **Conda** or Miniconda
- **SingleM** (conda environment)
- **Diamond** (conda environment with Trimmomatic)
- **Python 3** with packages:
  - pandas
  - numpy
  - openpyxl
  - tqdm

### Databases

1. **SingleM Metapackage** (e.g., S5.4.0.GTDB_r226)
2. **KEGG Prokaryotes Database** (.dmnd format for Diamond)
3. **KEGG Annotation File** (prokaryotes.dat)

### Conda Environments

#### SingleM Environment
```bash
conda create -n singlem -c bioconda singlem
```

#### Diamond/Metagenome Assembly Environment
```bash
conda create -n metagenome_assembly -c bioconda diamond trimmomatic pigz
```

#### Python Environment (if not in base)
```bash
conda install pandas numpy openpyxl tqdm
```

## Installation

1. **Clone the repository**
```bash
git clone https://github.com/Naturalist1986/Metagenomic_Kegg_Mapping_and_Anlaysis_Pipeline.git
cd Metagenomic_Kegg_Mapping_and_Anlaysis_Pipeline
```

2. **Run setup script**
```bash
chmod +x setup_pipeline.sh
./setup_pipeline.sh
```

This will:
- Create the `logs/` directory (required for SLURM output)
- Make all scripts executable
- Verify directory structure

3. **Create configuration file**
```bash
cp config.example.sh config.sh
```

4. **Edit configuration file**
```bash
nano config.sh  # or your preferred editor
```

Update the following required fields:
- `FASTQ_DIR` - Path to your raw FASTQ files
- `OUTPUT_BASE_DIR` - Where to write results
- `PIPELINE_DIR` - Path to this pipeline directory
- `SINGLEM_METAPACKAGE_PATH` - Path to SingleM database
- `KEGG_DB_PATH` - Path to KEGG Diamond database
- `KEGG_DAT_FILE` - Path to KEGG annotation file

## Usage

### Basic Usage

```bash
# Submit the entire pipeline
sbatch run_metagenomics_pipeline.sh config.sh
```

### Monitor Progress

```bash
# Check job status
squeue -u $USER

# View logs
tail -f logs/pipeline_*.out

# Check individual step logs
tail -f logs/singlem_*.out
tail -f logs/diamond_*.out
```

### Treatment Groups (Optional)

To analyze samples by treatment groups, edit your `config.sh`:

```bash
# Example: Control vs Treatment
TREATMENT_GROUPS="control:sample1,sample2,sample3;treatment:sample4,sample5,sample6"
```

This creates combined SingleM summaries for each treatment group at all taxonomic levels.

### Running Individual Steps

You can also run individual pipeline steps:

```bash
# Step 1: SingleM only
sbatch --array=0-N scripts/run_singlem_step.sh config.sh

# Step 3: Diamond alignment only (after SingleM)
sbatch --array=0-N scripts/run_diamond_step.sh config.sh
```

## Input Data

### FASTQ Files

Place paired-end FASTQ files in your `FASTQ_DIR`. Supported naming conventions:

- `SAMPLE_1.fq.gz` and `SAMPLE_2.fq.gz`
- `SAMPLE_1.fastq.gz` and `SAMPLE_2.fastq.gz`
- `SAMPLE_R1.fq.gz` and `SAMPLE_R2.fq.gz`
- `SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`

**Example:**
```
/path/to/fastqs/
├── control_rep1_1.fq.gz
├── control_rep1_2.fq.gz
├── control_rep2_1.fq.gz
├── control_rep2_2.fq.gz
├── treatment_rep1_1.fq.gz
├── treatment_rep1_2.fq.gz
└── ...
```

## Output Structure

```
OUTPUT_BASE_DIR/
├── singlem/
│   ├── output/                    # SingleM taxonomic profiles
│   ├── otu-table/                 # OTU tables
│   ├── krona/                     # Krona HTML visualizations
│   ├── summarize_byphylum/        # Phylum-level summaries
│   ├── summarize_byclass/         # Class-level summaries
│   ├── summarize_byorder/         # Order-level summaries
│   ├── summarize_byfamily/        # Family-level summaries
│   ├── summarize_bygenus/         # Genus-level summaries
│   ├── summarize_byspecies/       # Species-level summaries
│   ├── combined_by_treatment/     # Treatment group summaries
│   └── prok_fraction/             # Prokaryotic fraction estimates
├── diamond_kegg/
│   ├── *.tsv.gz                   # Diamond BLASTX results
│   ├── qc_filtered_trimmomatic/   # QC-filtered FASTQ files
│   └── trimmomatic_logs/          # Trimmomatic logs
├── post_diamond_processing/
│   └── *_kegg_hits.tsv           # Extracted KEGG hits
├── kegg_summed/
│   ├── *_kegg_hits_summed.tsv    # Summed KEGG hits per sample
│   └── kegg_stats.xlsx            # Genome equivalent estimates
└── final_results/
    └── normalized_kegg_feature_table.tsv  # Final normalized table
```

## Output Files

### Key Output Files

1. **Normalized KEGG Feature Table** (`final_results/normalized_kegg_feature_table.tsv`)
   - Rows: KEGG orthology (KO) numbers
   - Columns: Sample names
   - Values: Normalized hit counts (hits per genome equivalent)
   - **Primary file for downstream statistical analysis**

2. **Treatment Group Summaries** (`singlem/combined_by_treatment/*.csv`)
   - Combined taxonomic profiles for each treatment group
   - Separate files for each taxonomic level (phylum, class, order, family, genus, species)

3. **Genome Equivalent Estimates** (`kegg_summed/kegg_stats.xlsx`)
   - Sample-level genome equivalent calculations
   - Based on mean single-copy gene abundance

4. **SingleM Taxonomic Profiles** (`singlem/output/*.tsv`)
   - Per-sample taxonomic classification
   - Used for community composition analysis

## Pipeline Steps Detail

### Step 1: SingleM Taxonomic Profiling

- **Input:** Raw FASTQ files
- **Output:** Taxonomic profiles, OTU tables, Krona plots
- **Resources:** 32 cores, 32GB RAM, 6 hours per sample
- **Key Features:**
  - Ribosomal protein gene profiling
  - Taxonomic classification
  - Prokaryotic fraction estimation

### Step 2: Treatment Group Summarization

- **Input:** Individual SingleM profiles
- **Output:** Combined summaries by treatment group
- **Resources:** 8 cores, 16GB RAM, 2 hours
- **Taxonomic Levels:**
  - Phylum
  - Class
  - Order
  - Family
  - Genus
  - Species

### Step 3: Diamond KEGG Alignment

- **Input:** Raw FASTQ files
- **Output:** KEGG database alignments
- **Resources:** 32 cores, 512GB RAM, 7 days per sample
- **Steps:**
  1. Trimmomatic quality filtering
  2. Diamond BLASTX alignment (R1 and R2 separately)
  3. Combination and compression

### Step 4-7: Post-processing

- **Step 4:** Extract KEGG orthology numbers from alignments
- **Step 5:** Sum hits per KO per sample
- **Step 6:** Calculate genome equivalents using single-copy genes
- **Step 7:** Normalize by genome equivalents

## Troubleshooting

### Common Issues

**Issue:** `No log file created after sbatch`
- **Solution:** The `logs/` directory must exist before submitting jobs
- Run `./setup_pipeline.sh` to create it
- Or manually: `mkdir -p logs`

**Issue:** `No FASTQ files found`
- **Solution:** Check that FASTQ files match expected naming patterns
- Verify `FASTQ_DIR` path in config.sh

**Issue:** `KEGG database not found`
- **Solution:** Ensure KEGG database is in Diamond format (.dmnd)
- Run `diamond makedb` if needed

**Issue:** `Conda environment not found`
- **Solution:** Verify environment names in config.sh match your conda environments
- Run `conda env list` to check

**Issue:** `Job failed with memory error`
- **Solution:** Increase memory allocation in SBATCH directives
- Check cluster limits with `sinfo`

**Issue:** `Treatment groups not working`
- **Solution:** Ensure sample names in TREATMENT_GROUPS exactly match FASTQ file prefixes
- Check format: `"group1:sample1,sample2;group2:sample3,sample4"`

### Log Files

All steps generate detailed log files in `logs/`:
- `pipeline_*.out` - Main pipeline log
- `singlem_*.out` - SingleM step logs
- `diamond_*.out` - Diamond step logs
- `after_diamond_*.out` - Post-processing logs

Check these files for detailed error messages.

## Citation

If you use this pipeline, please cite:

- **SingleM:** [Insert SingleM citation]
- **Diamond:** Buchfink B, Xie C, Huson DH. Fast and sensitive protein alignment using DIAMOND. Nature Methods 12, 59-60 (2015).
- **KEGG:** Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## License

[Insert License Information]

## Support

For issues and questions:
- Open an issue on GitHub
- Check existing issues for solutions
- Review log files for detailed error messages

## Authors

[Insert Author Information]

---

**Last Updated:** 2025-11-19
**Version:** 1.0.0
