#!/bin/bash
#SBATCH --job-name=metagenomics_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=7-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

################################################################################
# METAGENOMICS KEGG MAPPING AND ANALYSIS PIPELINE
#
# This master script orchestrates the complete metagenomics analysis pipeline:
# 1. SingleM taxonomic profiling on raw FASTQ files
# 2. SingleM treatment group summarization at all taxonomic levels
# 3. Diamond KEGG alignment (with Trimmomatic preprocessing)
# 4. Post-processing and normalization
#
# Usage:
#   sbatch run_metagenomics_pipeline.sh config.sh
#
# The config file should define all necessary paths and parameters
################################################################################

set -euo pipefail

# Function for logging with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_message "=========================================="
log_message "METAGENOMICS PIPELINE STARTED"
log_message "=========================================="

# Check if config file is provided
if [ $# -eq 0 ]; then
    log_message "ERROR: No configuration file provided"
    log_message "Usage: sbatch $0 config.sh"
    exit 1
fi

CONFIG_FILE=$1

# Source the configuration file
if [ ! -f "$CONFIG_FILE" ]; then
    log_message "ERROR: Configuration file not found: $CONFIG_FILE"
    exit 1
fi

source "$CONFIG_FILE"
log_message "Configuration loaded from: $CONFIG_FILE"

# Validate required variables from config
# Always required
if [ -z "${FASTQ_DIR:-}" ]; then
    log_message "ERROR: FASTQ_DIR not set in config file"
    exit 1
fi

if [ -z "${OUTPUT_BASE_DIR:-}" ]; then
    log_message "ERROR: OUTPUT_BASE_DIR not set in config file"
    exit 1
fi

if [ -z "${PIPELINE_DIR:-}" ]; then
    log_message "ERROR: PIPELINE_DIR not set in config file"
    exit 1
fi

# Only required if SingleM is not skipped
if [ "${SKIP_SINGLEM:-false}" != "true" ]; then
    if [ -z "${SINGLEM_METAPACKAGE_PATH:-}" ]; then
        log_message "ERROR: SINGLEM_METAPACKAGE_PATH not set (required when SKIP_SINGLEM is not true)"
        exit 1
    fi
fi

# Only required if Diamond is not skipped
if [ "${SKIP_DIAMOND:-false}" != "true" ]; then
    if [ -z "${KEGG_DB_PATH:-}" ]; then
        log_message "ERROR: KEGG_DB_PATH not set (required when SKIP_DIAMOND is not true)"
        exit 1
    fi
    if [ -z "${KEGG_DAT_FILE:-}" ]; then
        log_message "ERROR: KEGG_DAT_FILE not set (required when SKIP_DIAMOND is not true)"
        exit 1
    fi
fi

# Create log directory
mkdir -p logs

# Set default number of samples if not provided
NUM_SAMPLES=${NUM_SAMPLES:-$(ls ${FASTQ_DIR}/*_1.f*q.gz 2>/dev/null | wc -l)}

if [ "$NUM_SAMPLES" -eq 0 ]; then
    log_message "ERROR: No FASTQ files found in $FASTQ_DIR"
    exit 1
fi

log_message "Found $NUM_SAMPLES samples to process"

# Initialize variables for job dependencies
SINGLEM_JOB_ID=""
SUMMARIZE_JOB_ID=""
DIAMOND_JOB_ID=""

################################################################################
# AUTO-RESUME: Check for completed stages
################################################################################

if [ "${AUTO_RESUME:-true}" = "true" ]; then
    log_message "Auto-resume enabled: Checking for completed stages..."

    # Check if SingleM is complete (unless already set to skip)
    if [ "${SKIP_SINGLEM:-false}" != "true" ]; then
        SINGLEM_OUTDIR="${OUTPUT_BASE_DIR}/singlem/output"
        if [ -d "$SINGLEM_OUTDIR" ]; then
            # Count completed SingleM outputs
            COMPLETED_SINGLEM=$(ls -1 ${SINGLEM_OUTDIR}/*.tsv 2>/dev/null | wc -l)
            if [ "$COMPLETED_SINGLEM" -eq "$NUM_SAMPLES" ]; then
                log_message "  ✓ Found $COMPLETED_SINGLEM/$NUM_SAMPLES completed SingleM outputs"
                log_message "  → Auto-skipping SingleM stage"
                SKIP_SINGLEM=true
            else
                log_message "  ○ SingleM incomplete: $COMPLETED_SINGLEM/$NUM_SAMPLES samples"
            fi
        fi
    fi

    # Check if Diamond is complete (unless already set to skip)
    if [ "${SKIP_DIAMOND:-false}" != "true" ]; then
        DIAMOND_OUTDIR="${OUTPUT_BASE_DIR}/diamond_kegg"
        if [ -d "$DIAMOND_OUTDIR" ]; then
            # Count completed Diamond outputs
            COMPLETED_DIAMOND=$(ls -1 ${DIAMOND_OUTDIR}/*.tsv.gz 2>/dev/null | wc -l)
            if [ "$COMPLETED_DIAMOND" -eq "$NUM_SAMPLES" ]; then
                log_message "  ✓ Found $COMPLETED_DIAMOND/$NUM_SAMPLES completed Diamond outputs"
                log_message "  → Auto-skipping Diamond stage"
                SKIP_DIAMOND=true
            else
                log_message "  ○ Diamond incomplete: $COMPLETED_DIAMOND/$NUM_SAMPLES samples"
            fi
        fi
    fi

    log_message "=========================================="
fi

################################################################################
# STEP 1: SINGLEM TAXONOMIC PROFILING
################################################################################

if [ "${SKIP_SINGLEM:-false}" = "true" ]; then
    log_message "=========================================="
    log_message "STEP 1: SingleM - SKIPPED (SKIP_SINGLEM=true)"
    log_message "=========================================="
    log_message "STEP 2: SingleM Summarization - SKIPPED"
    log_message "=========================================="
else
    log_message "=========================================="
    log_message "STEP 1: Running SingleM Taxonomic Profiling"
    log_message "=========================================="

    SINGLEM_JOB_ID=$(sbatch \
        --parsable \
        --array=0-$((NUM_SAMPLES-1)) \
        --output=logs/singlem_%A_%a.out \
        --error=logs/singlem_%A_%a.err \
        ${PIPELINE_DIR}/scripts/run_singlem_step.sh ${CONFIG_FILE})

    log_message "SingleM job submitted: Job ID $SINGLEM_JOB_ID"

    ############################################################################
    # STEP 2: SINGLEM TREATMENT GROUP SUMMARIZATION
    ############################################################################

    log_message "=========================================="
    log_message "STEP 2: SingleM Treatment Group Summarization"
    log_message "=========================================="

    SUMMARIZE_JOB_ID=$(sbatch \
        --parsable \
        --dependency=afterok:${SINGLEM_JOB_ID} \
        --output=logs/singlem_summarize_%j.out \
        --error=logs/singlem_summarize_%j.err \
        ${PIPELINE_DIR}/scripts/singlem_treatment_summarize.sh ${CONFIG_FILE})

    log_message "SingleM summarization job submitted: Job ID $SUMMARIZE_JOB_ID"
fi

################################################################################
# STEP 3: DIAMOND KEGG ALIGNMENT
################################################################################

if [ "${SKIP_DIAMOND:-false}" = "true" ]; then
    log_message "=========================================="
    log_message "STEP 3: Diamond KEGG - SKIPPED (SKIP_DIAMOND=true)"
    log_message "=========================================="
else
    log_message "=========================================="
    log_message "STEP 3: Running Diamond KEGG Alignment"
    log_message "=========================================="

    # Set dependency based on whether SingleM was run
    if [ -n "$SINGLEM_JOB_ID" ]; then
        DEPENDENCY="--dependency=afterok:${SINGLEM_JOB_ID}"
    else
        DEPENDENCY=""
    fi

    DIAMOND_JOB_ID=$(sbatch \
        --parsable \
        --array=0-$((NUM_SAMPLES-1)) \
        $DEPENDENCY \
        --output=logs/diamond_%A_%a.out \
        --error=logs/diamond_%A_%a.err \
        ${PIPELINE_DIR}/scripts/run_diamond_step.sh ${CONFIG_FILE})

    log_message "Diamond job submitted: Job ID $DIAMOND_JOB_ID"
fi

################################################################################
# STEP 4: POST-DIAMOND KEGG PROCESSING
################################################################################

# Only run KEGG processing steps if Diamond was run
if [ -z "$DIAMOND_JOB_ID" ]; then
    log_message "=========================================="
    log_message "STEP 4-7: KEGG Processing - SKIPPED (Diamond not run)"
    log_message "=========================================="
    AFTER_DIAMOND_JOB_ID=""
    SUM_KEGG_JOB_ID=""
    MEAN_SC_JOB_ID=""
    NORMALIZE_JOB_ID=""
else
    log_message "=========================================="
    log_message "STEP 4: Post-Diamond KEGG Hit Processing"
    log_message "=========================================="

    AFTER_DIAMOND_JOB_ID=$(sbatch \
        --parsable \
        --array=0-$((NUM_SAMPLES-1)) \
        --dependency=afterok:${DIAMOND_JOB_ID} \
        --output=logs/after_diamond_%A_%a.out \
        --error=logs/after_diamond_%A_%a.err \
        ${PIPELINE_DIR}/scripts/run_after_diamond_step.sh ${CONFIG_FILE})

    log_message "After-Diamond processing job submitted: Job ID $AFTER_DIAMOND_JOB_ID"

    ############################################################################
    # STEP 5: SUM KEGG HITS
    ############################################################################

    log_message "=========================================="
    log_message "STEP 5: Summing KEGG Hits"
    log_message "=========================================="

    SUM_KEGG_JOB_ID=$(sbatch \
        --parsable \
        --array=0-$((NUM_SAMPLES-1)) \
        --dependency=afterok:${AFTER_DIAMOND_JOB_ID} \
        --output=logs/sum_kegg_%A_%a.out \
        --error=logs/sum_kegg_%A_%a.err \
        ${PIPELINE_DIR}/scripts/run_sum_kegg_step.sh ${CONFIG_FILE})

    log_message "Sum KEGG hits job submitted: Job ID $SUM_KEGG_JOB_ID"

    ############################################################################
    # STEP 6: CALCULATE MEAN SINGLE COPY GENES
    ############################################################################

    log_message "=========================================="
    log_message "STEP 6: Calculating Mean Single-Copy Genes"
    log_message "=========================================="

    MEAN_SC_JOB_ID=$(sbatch \
        --parsable \
        --dependency=afterok:${SUM_KEGG_JOB_ID} \
        --output=logs/mean_single_copy_%j.out \
        --error=logs/mean_single_copy_%j.err \
        ${PIPELINE_DIR}/scripts/run_mean_single_copy_step.sh ${CONFIG_FILE})

    log_message "Mean single-copy calculation job submitted: Job ID $MEAN_SC_JOB_ID"

    ############################################################################
    # STEP 7: CREATE NORMALIZED FEATURE TABLE
    ############################################################################

    log_message "=========================================="
    log_message "STEP 7: Creating Normalized Feature Table"
    log_message "=========================================="

    NORMALIZE_JOB_ID=$(sbatch \
        --parsable \
        --dependency=afterok:${MEAN_SC_JOB_ID} \
        --output=logs/normalize_%j.out \
        --error=logs/normalize_%j.err \
        ${PIPELINE_DIR}/scripts/run_normalize_step.sh ${CONFIG_FILE})

    log_message "Normalization job submitted: Job ID $NORMALIZE_JOB_ID"
fi

################################################################################
# PIPELINE SUMMARY
################################################################################

log_message "=========================================="
log_message "PIPELINE SUBMISSION COMPLETE"
log_message "=========================================="
log_message "Job dependency chain:"

if [ -n "$SINGLEM_JOB_ID" ]; then
    log_message "  1. SingleM:              $SINGLEM_JOB_ID"
    log_message "  2. SingleM Summarize:    $SUMMARIZE_JOB_ID"
else
    log_message "  1. SingleM:              SKIPPED"
    log_message "  2. SingleM Summarize:    SKIPPED"
fi

if [ -n "$DIAMOND_JOB_ID" ]; then
    log_message "  3. Diamond KEGG:         $DIAMOND_JOB_ID"
    log_message "  4. After Diamond:        $AFTER_DIAMOND_JOB_ID"
    log_message "  5. Sum KEGG Hits:        $SUM_KEGG_JOB_ID"
    log_message "  6. Mean Single-Copy:     $MEAN_SC_JOB_ID"
    log_message "  7. Normalize:            $NORMALIZE_JOB_ID"
else
    log_message "  3. Diamond KEGG:         SKIPPED"
    log_message "  4-7. KEGG Processing:    SKIPPED"
fi

log_message ""
log_message "Monitor progress with: squeue -u $USER"
log_message "Check logs in: logs/"
log_message "=========================================="

exit 0
