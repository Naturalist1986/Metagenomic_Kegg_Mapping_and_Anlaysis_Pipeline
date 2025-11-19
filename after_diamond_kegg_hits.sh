#!/bin/bash
#SBATCH --job-name=process_kegg
#SBATCH --output=/sci/backup/ofinkel/moshea/scripts/logs/job_%A_%a.out
#SBATCH --error=/sci/backup/ofinkel/moshea/scripts/logs/job_%A_%a.err
#SBATCH --array=0-33
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00

# Exit on any error, undefined variable, or pipe failure
set -euo pipefail

# Load modules if needed
# module load spack/all
# module load gzip/1.13-x86_64-gcc-12.2.0-oqikzf5

# Define directories and files
INPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_KEGG_output/"
OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_post_processing/"
KEGG_FILE="/sci/backup/aerez/aerez/moshea/kegg_2023/prokaryotes.dat"
SCRIPT="/sci/backup/ofinkel/moshea/scripts/process_single_diamond.py"
LOG_DIR="/sci/backup/ofinkel/moshea/scripts/logs"

# Function for logging with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Verify required directories and files exist
if [ ! -d "$INPUT_DIR" ]; then
    log_message "ERROR: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    log_message "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

if [ ! -d "$LOG_DIR" ]; then
    log_message "Creating log directory: $LOG_DIR"
    mkdir -p "$LOG_DIR"
fi

if [ ! -f "$KEGG_FILE" ]; then
    log_message "ERROR: KEGG file does not exist: $KEGG_FILE"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    log_message "ERROR: Python script does not exist: $SCRIPT"
    exit 1
fi

# Get list of files to process
shopt -s nullglob  # Prevent expansion if no matches
FILES=(${INPUT_DIR}*.tsv.gz)
shopt -u nullglob

# Check if any files were found
if [ ${#FILES[@]} -eq 0 ]; then
    log_message "ERROR: No .tsv.gz files found in $INPUT_DIR"
    exit 1
fi

# Validate array task ID is within bounds
if [ $SLURM_ARRAY_TASK_ID -ge ${#FILES[@]} ]; then
    log_message "ERROR: Array task ID $SLURM_ARRAY_TASK_ID exceeds number of files (${#FILES[@]})"
    exit 1
fi

# Get the file to process
FILE_TO_PROCESS=${FILES[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$FILE_TO_PROCESS" .tsv.gz)
OUTPUT_FILE="${OUTPUT_DIR}${BASENAME}_kegg_hits.tsv"
TEMP_OUTPUT="${OUTPUT_FILE}.tmp"
LOCK_FILE="${OUTPUT_FILE}.lock"

log_message "Starting job for task $SLURM_ARRAY_TASK_ID"
log_message "Processing file: $FILE_TO_PROCESS"
log_message "Output will be: $OUTPUT_FILE"

# Check if output file already exists and is complete
if [ -f "$OUTPUT_FILE" ]; then
    # Check if file has content (not just created but empty)
    if [ -s "$OUTPUT_FILE" ]; then
        log_message "Output file already exists and has content. Skipping."
        exit 0
    else
        log_message "Output file exists but is empty. Removing and reprocessing."
        rm -f "$OUTPUT_FILE"
    fi
fi

# Use file locking to prevent concurrent processing of same file
(
    flock -n 200 || {
        log_message "Another job is already processing this file. Exiting."
        exit 0
    }
    
    # Clean up any previous temporary files
    rm -f "$TEMP_OUTPUT"
    
    # Log system resources before processing
    log_message "Memory available: $(free -h | grep Mem | awk '{print $7}')"
    log_message "Starting processing with $SLURM_CPUS_PER_TASK CPUs"
    
    # Run the processing script with error handling
    if python "$SCRIPT" --input "$FILE_TO_PROCESS" --output "$TEMP_OUTPUT" --kegg "$KEGG_FILE"; then
        # Move temporary file to final location only if successful
        mv "$TEMP_OUTPUT" "$OUTPUT_FILE"
        log_message "Successfully processed $BASENAME"
        
        # Optional: Validate output file
        if [ ! -s "$OUTPUT_FILE" ]; then
            log_message "WARNING: Output file is empty after processing"
            exit 1
        fi
        
        # Log output file size for verification
        OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        log_message "Output file size: $OUTPUT_SIZE"
    else
        EXIT_CODE=$?
        log_message "ERROR: Processing failed with exit code $EXIT_CODE"
        rm -f "$TEMP_OUTPUT"  # Clean up temporary file
        exit $EXIT_CODE
    fi
    
) 200>"$LOCK_FILE"

# Clean up lock file
rm -f "$LOCK_FILE"

log_message "Job completed successfully for task $SLURM_ARRAY_TASK_ID"