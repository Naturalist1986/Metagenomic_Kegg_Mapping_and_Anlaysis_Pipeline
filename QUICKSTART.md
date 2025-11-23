# Quick Start Guide

Get your metagenomics pipeline running in 5 minutes!

## Prerequisites

- Access to an HPC cluster with SLURM
- Paired-end metagenomic FASTQ files
- SingleM and Diamond conda environments set up
- KEGG database downloaded

## 5-Step Setup

### 1. Clone and Enter Directory

```bash
git clone https://github.com/Naturalist1986/Metagenomic_Kegg_Mapping_and_Anlaysis_Pipeline.git
cd Metagenomic_Kegg_Mapping_and_Anlaysis_Pipeline
```

### 2. Run Setup Script

```bash
chmod +x setup_pipeline.sh
./setup_pipeline.sh
```

This creates the `logs/` directory and makes all scripts executable.

### 3. Create Configuration File

```bash
cp config.example.sh config.sh
nano config.sh  # or vim, emacs, etc.
```

### 4. Edit REQUIRED Fields in config.sh

**Minimum required changes:**

```bash
# 1. Where are your FASTQ files?
FASTQ_DIR="/path/to/your/raw_fastqs"

# 2. Where should results go?
OUTPUT_BASE_DIR="/path/to/your/output"

# 3. Where is this pipeline installed?
PIPELINE_DIR="/path/to/Metagenomic_Kegg_Mapping_and_Anlaysis_Pipeline"

# 4. Database paths
SINGLEM_METAPACKAGE_PATH="/path/to/singlem_database.smpkg.zb"
KEGG_DB_PATH="/path/to/prokaryotes.kegg.dmnd"
KEGG_DAT_FILE="/path/to/prokaryotes.dat"

# 5. Conda environment names (if different from defaults)
SINGLEM_ENV="singlem"           # your SingleM environment name
DIAMOND_ENV="metagenome_assembly"  # your Diamond environment name
```

**Optional - Treatment Groups:**

```bash
# Example: compare control vs treatment samples
TREATMENT_GROUPS="control:ctrl1,ctrl2,ctrl3;treatment:treat1,treat2,treat3"
```

### 5. Submit Pipeline

```bash
sbatch run_metagenomics_pipeline.sh config.sh
```

That's it! The pipeline will run all steps automatically.

**Pro tips:**

**Auto-Resume** (enabled by default):
- Pipeline automatically detects completed stages
- Just re-run the same command after Diamond finishes
- It will skip completed stages and continue from post-processing

**Manual control:**
```bash
# In config.sh:
AUTO_RESUME=true         # Auto-detect completed stages (default)
SKIP_SINGLEM=true        # Force skip SingleM
SKIP_DIAMOND=true        # Force skip Diamond
MAX_CONCURRENT_JOBS=10   # Limit concurrent jobs (optional)
```

**Resource management:**
- Set `MAX_CONCURRENT_JOBS=10` to run max 10 jobs at once
- Useful for large datasets or cluster policies
- Leave empty for no limit

## Monitor Progress

```bash
# Check job status
squeue -u $USER

# Watch main pipeline log
tail -f logs/pipeline_*.out

# Check specific steps
tail -f logs/singlem_*.out
tail -f logs/diamond_*.out
```

## Expected Timeline

For a typical metagenomic dataset:

| Step | Time per Sample | Total Time (10 samples) |
|------|----------------|------------------------|
| 1. SingleM | 2-6 hours | 2-6 hours (parallel) |
| 2. Summarize | N/A | 30 min - 1 hour |
| 3. Diamond | 12-72 hours | 12-72 hours (parallel) |
| 4-7. Post-processing | 1-3 hours | 1-3 hours (sequential) |

**Total: ~15-80 hours** depending on sample size and complexity

## Quick Troubleshooting

### "No FASTQ files found"
- Check your FASTQ files are named: `SAMPLE_1.fq.gz` and `SAMPLE_2.fq.gz`
- Or: `SAMPLE_R1.fastq.gz` and `SAMPLE_R2.fastq.gz`
- Verify `FASTQ_DIR` path is correct

### "Database not found"
- Verify database paths exist: `ls -lh /path/to/database`
- KEGG database must be in Diamond format (.dmnd)

### "Conda environment not found"
- Check environment names: `conda env list`
- Update `SINGLEM_ENV` and `DIAMOND_ENV` in config.sh

### Job Failed
- Check logs: `cat logs/*_ERROR_*.err`
- Most common issues:
  - Incorrect paths in config.sh
  - Insufficient memory (increase in SBATCH directives)
  - Database format mismatch

## What You Get

After successful completion, find your results in:

```
OUTPUT_BASE_DIR/
├── singlem/                          # Taxonomic profiles
│   └── combined_by_treatment/        # Treatment group summaries ⭐
├── diamond_kegg/                     # Raw alignments
├── kegg_summed/
│   └── kegg_stats.xlsx              # Genome equivalents
└── final_results/
    └── normalized_kegg_feature_table.tsv  # YOUR MAIN RESULT ⭐⭐⭐
```

**Main output:** `normalized_kegg_feature_table.tsv`
- Rows = KEGG orthology numbers
- Columns = Your samples
- Values = Normalized functional gene abundances
- **Ready for statistical analysis!**

## Next Steps

1. **Quality Control**
   - Check SingleM Krona plots: `singlem/krona/*.html`
   - Review genome equivalent estimates: `kegg_summed/kegg_stats.xlsx`

2. **Statistical Analysis**
   - Import `normalized_kegg_feature_table.tsv` into R/Python
   - Perform differential abundance analysis
   - Create visualizations

3. **Biological Interpretation**
   - Map significant KEGGs to pathways
   - Use KEGG Mapper: https://www.genome.jp/kegg/mapper/
   - Interpret in context of your experimental design

## Example Analysis in R

```r
# Load normalized feature table
kegg_table <- read.table("final_results/normalized_kegg_feature_table.tsv",
                         header = TRUE, sep = "\t", row.names = 1)

# Basic exploration
dim(kegg_table)  # Number of KEGGs x samples
summary(kegg_table)

# Simple PCA
pca <- prcomp(t(kegg_table), scale. = TRUE)
plot(pca$x[,1], pca$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "KEGG Functional Profile PCA")
```

## Getting Help

- **Full documentation:** See README.md
- **Issues:** Check logs in `logs/` directory
- **GitHub Issues:** [Open an issue](https://github.com/yourusername/repo/issues)

---

**Happy analyzing! 🧬🔬**
