#!/bin/bash
#SBATCH --job-name=singlem
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -t 6:00:00
#SBATCH --mem=32G

################################################################################
# SingleM Taxonomic Profiling Step
#
# This script runs SingleM on individual FASTQ samples
# Handles various FASTQ file formats: .fq.gz, .fastq.gz, .fq, .fastq
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

# Initialize conda for this script
if [ -f "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh"
else
    log_message "ERROR: Cannot find conda initialization script"
    exit 1
fi

# Activate SingleM environment
conda activate ${SINGLEM_ENV:-singlem}
log_message "Activated conda environment: ${SINGLEM_ENV:-singlem}"

# Set SingleM metapackage path
export SINGLEM_METAPACKAGE_PATH="${SINGLEM_METAPACKAGE_PATH}"

# Calculate array index offset
OFFSET=${ARRAY_OFFSET:-0}
LINE_NUM=$((SLURM_ARRAY_TASK_ID + OFFSET))

# Find all FASTQ files (handles various formats)
# Priority: _1.fq.gz, _1.fastq.gz, _R1.fq.gz, _R1.fastq.gz
FASTQ_PATTERNS=(
    "${FASTQ_DIR}/*_1.fq.gz"
    "${FASTQ_DIR}/*_1.fastq.gz"
    "${FASTQ_DIR}/*_R1.fq.gz"
    "${FASTQ_DIR}/*_R1.fastq.gz"
)

# Try each pattern until we find files
Samples=()
for pattern in "${FASTQ_PATTERNS[@]}"; do
    files=( $pattern )
    if [ -e "${files[0]}" ]; then
        # Extract sample names based on pattern
        if [[ "$pattern" == *"_1.fq.gz" ]]; then
            Samples=( $(ls $pattern 2>/dev/null | xargs -n 1 basename | sed 's/_1\.fq\.gz$//') )
        elif [[ "$pattern" == *"_1.fastq.gz" ]]; then
            Samples=( $(ls $pattern 2>/dev/null | xargs -n 1 basename | sed 's/_1\.fastq\.gz$//') )
        elif [[ "$pattern" == *"_R1.fq.gz" ]]; then
            Samples=( $(ls $pattern 2>/dev/null | xargs -n 1 basename | sed 's/_R1\.fq\.gz$//') )
        elif [[ "$pattern" == *"_R1.fastq.gz" ]]; then
            Samples=( $(ls $pattern 2>/dev/null | xargs -n 1 basename | sed 's/_R1\.fastq\.gz$//') )
        fi
        break
    fi
done

# Check if samples were found
if [ ${#Samples[@]} -eq 0 ]; then
    log_message "ERROR: No paired-end FASTQ files found in $FASTQ_DIR"
    log_message "Tried patterns: ${FASTQ_PATTERNS[@]}"
    exit 1
fi

# Validate array index
if [ $LINE_NUM -ge ${#Samples[@]} ]; then
    log_message "ERROR: Array index $LINE_NUM exceeds number of samples (${#Samples[@]})"
    exit 1
fi

# Get sample name
SAMPLE=${Samples[$LINE_NUM]}
log_message "Processing sample: $SAMPLE (index: $LINE_NUM)"

# Auto-detect file extensions
F1=""
F2=""
for ext1 in "_1.fq.gz" "_1.fastq.gz" "_R1.fq.gz" "_R1.fastq.gz"; do
    ext2="${ext1/_1/_2}"
    ext2="${ext2/_R1/_R2}"
    test_f1="${FASTQ_DIR}/${SAMPLE}${ext1}"
    test_f2="${FASTQ_DIR}/${SAMPLE}${ext2}"

    if [ -e "$test_f1" ] && [ -e "$test_f2" ]; then
        F1="$test_f1"
        F2="$test_f2"
        log_message "Found FASTQ pair: $ext1 / $ext2"
        break
    fi
done

# Verify files exist
if [ ! -e "$F1" ] || [ ! -e "$F2" ]; then
    log_message "ERROR: FASTQ files not found for sample $SAMPLE"
    log_message "  Expected R1: $F1"
    log_message "  Expected R2: $F2"
    exit 1
fi

# Create output directories
SINGLEM_OUTDIR="${OUTPUT_BASE_DIR}/singlem"
mkdir -p $SINGLEM_OUTDIR/{output,otu-table,krona,archive-otu-table,summarize,summarize_byphylum,summarize_byclass,summarize_byorder,summarize_byfamily,summarize_bygenus,summarize_byspecies,prok_fraction/per-taxon-read-fraction}

# Define output files
OUTFILE="${SINGLEM_OUTDIR}/output/${SAMPLE}.tsv"

# Taxonomic level output files
SUMMARIZE_FILES=(
    "${SINGLEM_OUTDIR}/summarize/${SAMPLE}.tsv"
    "${SINGLEM_OUTDIR}/summarize_byphylum/${SAMPLE}.tsv:phylum"
    "${SINGLEM_OUTDIR}/summarize_byclass/${SAMPLE}.tsv:class"
    "${SINGLEM_OUTDIR}/summarize_byorder/${SAMPLE}.tsv:order"
    "${SINGLEM_OUTDIR}/summarize_byfamily/${SAMPLE}.tsv:family"
    "${SINGLEM_OUTDIR}/summarize_bygenus/${SAMPLE}.tsv:genus"
    "${SINGLEM_OUTDIR}/summarize_byspecies/${SAMPLE}.tsv:species"
)

# Run SingleM pipe
if [ -e "$OUTFILE" ]; then
    log_message "Output file exists, skipping: $OUTFILE"
else
    log_message "Running SingleM version: $(singlem --version)"
    singlem pipe \
        -1 "$F1" \
        -2 "$F2" \
        -p "$OUTFILE" \
        --otu-table "${SINGLEM_OUTDIR}/otu-table/${SAMPLE}-otu_table.tsv" \
        --threads ${SINGLEM_THREADS:-32} \
        --taxonomic-profile-krona "${SINGLEM_OUTDIR}/krona/${SAMPLE}.html" \
        --archive-otu-table "${SINGLEM_OUTDIR}/archive-otu-table/${SAMPLE}-archive-otu-table.tsv"

    if [ $? -ne 0 ]; then
        log_message "ERROR: SingleM pipe failed"
        exit 1
    fi
    log_message "SingleM pipe completed successfully"
fi

# Run summarization at each taxonomic level
for entry in "${SUMMARIZE_FILES[@]}"; do
    IFS=':' read -r outfile level <<< "$entry"

    if [ -e "$outfile" ]; then
        log_message "Summarize file exists, skipping: $outfile"
        continue
    fi

    if [ -z "$level" ]; then
        # Full taxonomy
        log_message "Running SingleM summarize (full taxonomy)"
        singlem summarise \
            --input-taxonomic-profile "$OUTFILE" \
            --output-taxonomic-profile-with-extras "$outfile"
    else
        # Specific taxonomic level
        log_message "Running SingleM summarize at level: $level"
        singlem summarise \
            --input-taxonomic-profile "$OUTFILE" \
            --output-species-by-site-relative-abundance "$outfile" \
            --output-species-by-site-level "$level"
    fi

    if [ $? -ne 0 ]; then
        log_message "ERROR: SingleM summarize failed for level: ${level:-full}"
        exit 1
    fi
done

# Calculate prokaryotic fraction
PROK_OUTFILE="${SINGLEM_OUTDIR}/prok_fraction/${SAMPLE}.tsv"

if [ -e "$PROK_OUTFILE" ]; then
    log_message "Prokaryotic fraction file exists, skipping: $PROK_OUTFILE"
else
    log_message "Calculating prokaryotic fraction"
    singlem microbial_fraction \
        -1 "$F1" \
        -2 "$F2" \
        -p "$OUTFILE" \
        --output-tsv "$PROK_OUTFILE" \
        --output-per-taxon-read-fractions "${SINGLEM_OUTDIR}/prok_fraction/per-taxon-read-fraction/${SAMPLE}-per-taxon-read-fracs.tsv"

    if [ $? -ne 0 ]; then
        log_message "ERROR: SingleM microbial_fraction failed"
        exit 1
    fi
fi

log_message "=========================================="
log_message "SingleM processing completed for: $SAMPLE"
log_message "=========================================="

exit 0
