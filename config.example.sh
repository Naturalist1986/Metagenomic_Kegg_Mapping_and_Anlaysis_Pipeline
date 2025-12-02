#!/bin/bash

################################################################################
# METAGENOMICS PIPELINE CONFIGURATION FILE
#
# Copy this file to config.sh and modify the paths and parameters
# for your specific analysis
#
# Usage: cp config.example.sh config.sh
#        Edit config.sh with your settings
#        sbatch run_metagenomics_pipeline.sh config.sh
################################################################################

# ============================================================================
# DIRECTORY PATHS
# ============================================================================

# Input directory containing raw FASTQ files
# Files should be paired-end with naming: SAMPLE_1.fq.gz and SAMPLE_2.fq.gz
# (or variations: .fastq.gz, _R1/_R2)
FASTQ_DIR="/path/to/your/raw_fastqs"

# Base output directory (all results will be organized under this)
OUTPUT_BASE_DIR="/path/to/your/output"

# Pipeline installation directory (where this config file is located)
PIPELINE_DIR="/path/to/Metagenomic_Kegg_Mapping_and_Anlaysis_Pipeline"

# ============================================================================
# DATABASE PATHS
# ============================================================================

# SingleM metapackage path
SINGLEM_METAPACKAGE_PATH="/path/to/singlem/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb"

# KEGG database for Diamond (should be .dmnd format)
KEGG_DB_PATH="/path/to/kegg/prokaryotes.kegg.dmnd"

# KEGG annotation file (prokaryotes.dat)
KEGG_DAT_FILE="/path/to/kegg/prokaryotes.dat"

# ============================================================================
# CONDA ENVIRONMENTS
# ============================================================================

# Conda base directory (default: ~/miniconda3)
CONDA_BASE="${HOME}/miniconda3"

# SingleM conda environment name
SINGLEM_ENV="singlem"

# Diamond/metagenome assembly conda environment name
DIAMOND_ENV="metagenome_assembly"

# ============================================================================
# PROCESSING PARAMETERS
# ============================================================================

# Skip pipeline stages (set to "true" to skip, "false" or empty to run)
SKIP_SINGLEM=false          # Skip SingleM taxonomic profiling
SKIP_DIAMOND=false          # Skip Diamond KEGG alignment

# Auto-resume from completed stages (checks for existing output files)
AUTO_RESUME=true            # Automatically detect and skip completed stages

# Number of samples to process (auto-detected if not set)
# Leave empty for automatic detection
NUM_SAMPLES=""

# Array offset (usually 0)
ARRAY_OFFSET=0

# Maximum number of concurrent jobs per array stage
# Limits how many jobs run simultaneously (e.g., "10" = max 10 jobs at once)
# Leave empty for no limit (all jobs can run in parallel)
# Example: MAX_CONCURRENT_JOBS=10
MAX_CONCURRENT_JOBS=""

# SingleM threads per sample
SINGLEM_THREADS=32

# Trimmomatic threads
TRIMMOMATIC_THREADS=32

# Diamond threads per sample
DIAMOND_THREADS=30

# Diamond block size (memory management)
# Lower value = less memory usage but slower
# Recommended: 2.0-6.0 for large files, leave empty for default
# Example: DIAMOND_BLOCK_SIZE=4.0
DIAMOND_BLOCK_SIZE=""

# Diamond memory allocation (SLURM --mem parameter)
# Adjust based on your sample sizes and available cluster resources
# Examples: 256G, 512G, 768G, 1T
# Default: 512G
DIAMOND_MEMORY="512G"

# Diamond time limit (SLURM --time parameter)
# Adjust based on sample size and cluster queue availability
# Format: days-hours:minutes:seconds or hours:minutes:seconds
# Examples: 7-0 (7 days), 3-12:00:00 (3.5 days), 48:00:00 (2 days)
# Default: 7-0 (7 days)
DIAMOND_TIME="7-0"

# Keep QC-filtered FASTQ files after Diamond processing
# Set to "true" to keep files, "false" to delete after use
KEEP_QC_FILES=false

# Trimmomatic JAR path (auto-detected if in conda environment)
# Only set this if you have Trimmomatic installed in a custom location
# TRIMMOMATIC_JAR="${CONDA_BASE}/envs/${DIAMOND_ENV}/share/trimmomatic-0.39-2/trimmomatic.jar"

# ============================================================================
# TREATMENT GROUPS (OPTIONAL)
# ============================================================================

# Define treatment groups for SingleM summarization
# Format: "group1:sample1,sample2,sample3;group2:sample4,sample5,sample6"
#
# Example:
# TREATMENT_GROUPS="control:ctrl_1,ctrl_2,ctrl_3;treated:treat_1,treat_2,treat_3"
#
# Leave empty to create a single combined summary for all samples
TREATMENT_GROUPS=""

# ============================================================================
# SLURM RESOURCE ALLOCATION
# ============================================================================

# SLURM account/allocation name (required by some clusters)
# Leave empty if your cluster doesn't require an account specification
# Example: SLURM_ACCOUNT="my_lab_allocation"
SLURM_ACCOUNT=""

# These settings can be customized based on your cluster requirements
# Uncomment and modify as needed

# Memory allocations (already set in individual scripts)
# SINGLEM_MEM="32G"
# DIAMOND_MEM="512G"
# POST_PROCESSING_MEM="128G"

# Time limits (already set in individual scripts)
# SINGLEM_TIME="6:00:00"
# DIAMOND_TIME="7-00:00:00"

# ============================================================================
# VALIDATION
# ============================================================================

# Validate that required paths exist
validate_config() {
    local errors=0

    if [ ! -d "$FASTQ_DIR" ]; then
        echo "ERROR: FASTQ directory does not exist: $FASTQ_DIR"
        errors=$((errors + 1))
    fi

    if [ ! -f "$SINGLEM_METAPACKAGE_PATH" ]; then
        echo "ERROR: SingleM metapackage not found: $SINGLEM_METAPACKAGE_PATH"
        errors=$((errors + 1))
    fi

    if [ ! -f "$KEGG_DB_PATH" ]; then
        echo "ERROR: KEGG database not found: $KEGG_DB_PATH"
        errors=$((errors + 1))
    fi

    if [ ! -f "$KEGG_DAT_FILE" ]; then
        echo "ERROR: KEGG annotation file not found: $KEGG_DAT_FILE"
        errors=$((errors + 1))
    fi

    if [ ! -d "$PIPELINE_DIR" ]; then
        echo "ERROR: Pipeline directory not found: $PIPELINE_DIR"
        errors=$((errors + 1))
    fi

    if [ $errors -gt 0 ]; then
        echo ""
        echo "Configuration validation failed with $errors error(s)"
        echo "Please fix the errors above before running the pipeline"
        return 1
    fi

    echo "Configuration validation successful!"
    return 0
}

# Uncomment to validate configuration when sourced
# validate_config
