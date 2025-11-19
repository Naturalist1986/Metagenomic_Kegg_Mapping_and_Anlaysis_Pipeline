#!/bin/bash
#SBATCH --job-name=mean_single_copy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00

################################################################################
# Mean Single-Copy Genes Calculation Step
#
# This script calculates the mean of single-copy KEGG genes across all samples
# to estimate genome equivalents
################################################################################

set -euo pipefail

# Function for logging
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Load configuration
if [ $# -eq 0 ]; then
    log_message "ERROR: No configuration file provided"
    exit 1
fi

CONFIG_FILE=$1
source "$CONFIG_FILE"

# Set up directories
SUM_KEGG_DIR="${OUTPUT_BASE_DIR}/kegg_summed"

log_message "=========================================="
log_message "Calculating Mean Single-Copy Genes"
log_message "=========================================="

# Check if input directory exists and has files
if [ ! -d "$SUM_KEGG_DIR" ]; then
    log_message "ERROR: Input directory not found: $SUM_KEGG_DIR"
    exit 1
fi

tsv_files=( ${SUM_KEGG_DIR}/*_kegg_hits_summed.tsv )
if [ ! -e "${tsv_files[0]}" ]; then
    log_message "ERROR: No summed KEGG files found in $SUM_KEGG_DIR"
    exit 1
fi

log_message "Found ${#tsv_files[@]} files to process"

OUTPUT_FILE="${SUM_KEGG_DIR}/kegg_stats.xlsx"

# Check if output already exists
if [ -f "$OUTPUT_FILE" ]; then
    log_message "Output file already exists: $OUTPUT_FILE"
    log_message "Removing and regenerating..."
    rm -f "$OUTPUT_FILE"
fi

# Update the mean_single_copy.py script with current paths
MEAN_SC_SCRIPT="${PIPELINE_DIR}/mean_single_copy.py"

if [ ! -f "$MEAN_SC_SCRIPT" ]; then
    log_message "ERROR: Python script not found: $MEAN_SC_SCRIPT"
    exit 1
fi

# Create a temporary modified version of the script with updated paths
TEMP_SCRIPT="/tmp/mean_single_copy_temp_${SLURM_JOB_ID:-$$}.py"

# Update the input/output directories in the script
sed -e "s|^input_dir = .*|input_dir = \"${SUM_KEGG_DIR}/\"|" \
    -e "s|^output_file = .*|output_file = \"${OUTPUT_FILE}\"|" \
    "$MEAN_SC_SCRIPT" > "$TEMP_SCRIPT"

log_message "Running mean single-copy gene calculation..."

python "$TEMP_SCRIPT"

if [ $? -ne 0 ]; then
    log_message "ERROR: Mean single-copy calculation failed"
    rm -f "$TEMP_SCRIPT"
    exit 1
fi

# Clean up temp script
rm -f "$TEMP_SCRIPT"

# Verify output exists
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    OUTPUT_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
    log_message "Success! Output file created: $OUTPUT_FILE"
    log_message "File size: $OUTPUT_SIZE"
else
    log_message "ERROR: Output file was not created or is empty"
    exit 1
fi

log_message "=========================================="
log_message "Mean Single-Copy Calculation Complete"
log_message "=========================================="

exit 0
