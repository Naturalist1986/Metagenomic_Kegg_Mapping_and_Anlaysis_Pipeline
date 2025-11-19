#!/bin/bash
#SBATCH --job-name=singlem_summarize
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 2:00:00
#SBATCH --mem=16G

################################################################################
# SingleM Treatment Group Summarization
#
# This script creates combined summaries for all samples in treatment groups
# at each taxonomic level (phylum, class, order, family, genus, species)
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

# Initialize conda
if [ -f "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh"
else
    log_message "ERROR: Cannot find conda initialization script"
    exit 1
fi

# Activate SingleM environment
conda activate ${SINGLEM_ENV:-singlem}
log_message "Activated conda environment: ${SINGLEM_ENV:-singlem}"

export SINGLEM_METAPACKAGE_PATH="${SINGLEM_METAPACKAGE_PATH}"

SINGLEM_OUTDIR="${OUTPUT_BASE_DIR}/singlem"
COMBINED_DIR="${SINGLEM_OUTDIR}/combined_by_treatment"

# Create output directory
mkdir -p "$COMBINED_DIR"

# Taxonomic levels to process
TAXONOMIC_LEVELS=("phylum" "class" "order" "family" "genus" "species")

log_message "=========================================="
log_message "Starting Treatment Group Summarization"
log_message "=========================================="

# Check if TREATMENT_GROUPS is defined in config
if [ -z "${TREATMENT_GROUPS:-}" ]; then
    log_message "WARNING: TREATMENT_GROUPS not defined in config"
    log_message "Creating a single combined summary for all samples"

    # Process each taxonomic level
    for level in "${TAXONOMIC_LEVELS[@]}"; do
        log_message "Processing taxonomic level: $level"

        input_dir="${SINGLEM_OUTDIR}/summarize_by${level}"
        output_file="${COMBINED_DIR}/all_samples_${level}.csv"

        # Check if input directory exists and has files
        if [ ! -d "$input_dir" ]; then
            log_message "WARNING: Input directory not found: $input_dir"
            continue
        fi

        tsv_files=( ${input_dir}/*.tsv )
        if [ ! -e "${tsv_files[0]}" ]; then
            log_message "WARNING: No TSV files found in $input_dir"
            continue
        fi

        log_message "Found ${#tsv_files[@]} samples for level: $level"
        log_message "Creating combined summary: $output_file"

        # Run SingleM summarise for all samples at this level
        singlem summarise \
            --input-taxonomic-profile ${input_dir}/*.tsv \
            --output-species-by-site-relative-abundance "$output_file" \
            --output-species-by-site-level "$level"

        if [ $? -eq 0 ]; then
            log_message "Successfully created combined summary for $level"
        else
            log_message "ERROR: Failed to create combined summary for $level"
            exit 1
        fi
    done
else
    # Process by treatment groups
    log_message "Processing samples by treatment groups"

    # Parse treatment groups (format: "group1:sample1,sample2;group2:sample3,sample4")
    IFS=';' read -ra GROUPS <<< "$TREATMENT_GROUPS"

    for group_entry in "${GROUPS[@]}"; do
        IFS=':' read -r group_name sample_list <<< "$group_entry"

        log_message "----------------------------------------"
        log_message "Processing treatment group: $group_name"
        log_message "----------------------------------------"

        # Convert comma-separated samples to array
        IFS=',' read -ra SAMPLES <<< "$sample_list"

        # Process each taxonomic level for this treatment group
        for level in "${TAXONOMIC_LEVELS[@]}"; do
            log_message "  Level: $level"

            input_dir="${SINGLEM_OUTDIR}/summarize_by${level}"
            output_file="${COMBINED_DIR}/${group_name}_${level}.csv"

            # Build list of input files for this group
            input_files=()
            for sample in "${SAMPLES[@]}"; do
                # Trim whitespace
                sample=$(echo "$sample" | xargs)
                tsv_file="${input_dir}/${sample}.tsv"

                if [ -e "$tsv_file" ]; then
                    input_files+=("$tsv_file")
                else
                    log_message "    WARNING: File not found for sample $sample: $tsv_file"
                fi
            done

            # Check if we have any files for this group
            if [ ${#input_files[@]} -eq 0 ]; then
                log_message "    WARNING: No input files found for group $group_name at level $level"
                continue
            fi

            log_message "    Found ${#input_files[@]} samples for group: $group_name"
            log_message "    Creating: $output_file"

            # Run SingleM summarise for this treatment group
            singlem summarise \
                --input-taxonomic-profile "${input_files[@]}" \
                --output-species-by-site-relative-abundance "$output_file" \
                --output-species-by-site-level "$level"

            if [ $? -eq 0 ]; then
                log_message "    ✓ Successfully created summary"
            else
                log_message "    ERROR: Failed to create summary for $group_name at level $level"
                exit 1
            fi
        done
    done
fi

log_message "=========================================="
log_message "Treatment Group Summarization Complete"
log_message "Output directory: $COMBINED_DIR"
log_message "=========================================="

exit 0
