#!/bin/bash
#SBATCH --job-name=diamond_kegg
#SBATCH -N 1
#SBATCH -c 32
# Memory and time are set dynamically via sbatch --mem and --time in main pipeline script

################################################################################
# Diamond KEGG Alignment Step (with Trimmomatic preprocessing)
#
# This script:
# 1. Runs Trimmomatic for quality filtering
# 2. Performs Diamond BLASTX alignment against KEGG database
# 3. Combines and compresses outputs
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

# Activate metagenome assembly environment
conda activate ${DIAMOND_ENV:-metagenome_assembly}
log_message "Activated conda environment: ${DIAMOND_ENV:-metagenome_assembly}"

# Verify Diamond is available
if ! command -v diamond &> /dev/null; then
    log_message "ERROR: Diamond not found in PATH"
    exit 1
fi

log_message "Diamond version: $(diamond --version)"

# Calculate array index
OFFSET=${ARRAY_OFFSET:-0}
LINE_NUM=$((SLURM_ARRAY_TASK_ID + OFFSET))

# Find all FASTQ files
FASTQ_PATTERNS=(
    "${FASTQ_DIR}/*_1.fq.gz"
    "${FASTQ_DIR}/*_1.fastq.gz"
    "${FASTQ_DIR}/*_R1.fq.gz"
    "${FASTQ_DIR}/*_R1.fastq.gz"
)

Samples=()
for pattern in "${FASTQ_PATTERNS[@]}"; do
    files=( $pattern )
    if [ -e "${files[0]}" ]; then
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

if [ ${#Samples[@]} -eq 0 ]; then
    log_message "ERROR: No FASTQ files found"
    exit 1
fi

if [ $LINE_NUM -ge ${#Samples[@]} ]; then
    log_message "ERROR: Array index $LINE_NUM exceeds number of samples"
    exit 1
fi

SAMPLE=${Samples[$LINE_NUM]}
log_message "=========================================="
log_message "Processing sample: $SAMPLE"
log_message "=========================================="

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
        break
    fi
done

if [ ! -e "$F1" ] || [ ! -e "$F2" ]; then
    log_message "ERROR: FASTQ files not found for sample $SAMPLE"
    exit 1
fi

# Set up directories
OUTDIR="${OUTPUT_BASE_DIR}/diamond_kegg"
QC_DIR="${OUTDIR}/qc_filtered_trimmomatic"
LOG_DIR="${OUTDIR}/trimmomatic_logs"
DIAMOND_TMP="${OUTDIR}/diamond_tmp"
TMP_DIR=${SLURM_TMPDIR:-/tmp/diamond_${SLURM_JOB_ID:-$$}}

mkdir -p "$OUTDIR" "$QC_DIR" "$LOG_DIR" "$DIAMOND_TMP" "$TMP_DIR"

# Define output files
QC1="${QC_DIR}/${SAMPLE}_R1_paired.fq.gz"
QC2="${QC_DIR}/${SAMPLE}_R2_paired.fq.gz"
QC1_UNPAIRED="${QC_DIR}/${SAMPLE}_R1_unpaired.fq.gz"
QC2_UNPAIRED="${QC_DIR}/${SAMPLE}_R2_unpaired.fq.gz"
OUT1="${OUTDIR}/${SAMPLE}_R1.tsv"
OUT2="${OUTDIR}/${SAMPLE}_R2.tsv"
OUTBOTH="${OUTDIR}/${SAMPLE}.tsv.gz"
TRIM_LOG="${LOG_DIR}/${SAMPLE}_trimmomatic.log"

# Check if final output exists
if [ -f "$OUTBOTH" ] && [ -s "$OUTBOTH" ]; then
    log_message "Final output already exists: $OUTBOTH"
    log_message "Skipping sample $SAMPLE"
    exit 0
fi

# ============================================================================
# STEP 1: TRIMMOMATIC QUALITY FILTERING
# ============================================================================

log_message "----------------------------------------"
log_message "Step 1: Quality filtering with Trimmomatic"
log_message "----------------------------------------"

QC_FILES_CREATED=false
if [ -f "$QC1" ] && [ -f "$QC2" ] && [ -s "$QC1" ] && [ -s "$QC2" ]; then
    log_message "QC files already exist. Skipping Trimmomatic."
else
    log_message "Running Trimmomatic..."
    QC_FILES_CREATED=true

    # Set Trimmomatic JAR path
    TRIMMOMATIC_JAR=${TRIMMOMATIC_JAR:-${CONDA_PREFIX}/share/trimmomatic*/trimmomatic.jar}

    # Find the actual JAR file if wildcard was used
    if [[ "$TRIMMOMATIC_JAR" == *"*"* ]]; then
        TRIMMOMATIC_JAR=$(ls ${CONDA_PREFIX}/share/trimmomatic*/trimmomatic.jar 2>/dev/null | head -1)
    fi

    if [ ! -f "$TRIMMOMATIC_JAR" ]; then
        log_message "ERROR: Trimmomatic JAR not found at $TRIMMOMATIC_JAR"
        exit 1
    fi

    log_message "Using Trimmomatic JAR: $TRIMMOMATIC_JAR"

    # Set adapter file
    ADAPTER_DIR=$(dirname $TRIMMOMATIC_JAR)/adapters
    if [ -d "$ADAPTER_DIR" ]; then
        ADAPTER_FILE="$ADAPTER_DIR/TruSeq3-PE-2.fa"
    else
        ADAPTER_FILE="${TMP_DIR}/adapters.fa"
        cat > "$ADAPTER_FILE" << 'EOF'
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
EOF
    fi

    # Run Trimmomatic
    java -jar "$TRIMMOMATIC_JAR" PE \
        -threads ${TRIMMOMATIC_THREADS:-32} \
        -phred33 \
        "$F1" "$F2" \
        "$QC1" "$QC1_UNPAIRED" \
        "$QC2" "$QC2_UNPAIRED" \
        ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:20 \
        MINLEN:50 \
        2> "$TRIM_LOG"

    if [ $? -ne 0 ]; then
        log_message "ERROR: Trimmomatic failed"
        exit 1
    fi

    # Clean up unpaired reads
    rm -f "$QC1_UNPAIRED" "$QC2_UNPAIRED"
    log_message "Trimmomatic completed successfully"
fi

# ============================================================================
# STEP 2: DIAMOND KEGG ALIGNMENT
# ============================================================================

log_message "----------------------------------------"
log_message "Step 2: Diamond KEGG alignment"
log_message "----------------------------------------"

# Function to run Diamond on a single file
run_diamond() {
    local input_file=$1
    local output_file=$2
    local thread_count=${DIAMOND_THREADS:-30}

    if [ -f "$output_file" ]; then
        rm -f "$output_file"
    fi

    log_message "Running Diamond on $(basename $input_file)..."

    # Build Diamond command with optional block size
    DIAMOND_CMD="diamond blastx \
        --query $input_file \
        --db $KEGG_DB_PATH \
        --max-target-seqs 1 \
        --evalue 1e-5 \
        --verbose \
        --tmpdir $DIAMOND_TMP \
        --out $output_file \
        --threads $thread_count \
        --outfmt 6 \
        --sensitive \
        --ignore-warnings"

    # Add block size if configured (helps with memory management)
    if [ -n "${DIAMOND_BLOCK_SIZE:-}" ]; then
        DIAMOND_CMD="$DIAMOND_CMD --block-size ${DIAMOND_BLOCK_SIZE}"
        log_message "  Using block size: ${DIAMOND_BLOCK_SIZE}"
    fi

    eval "$DIAMOND_CMD"

    if [ $? -ne 0 ]; then
        log_message "ERROR: Diamond failed for $(basename $input_file)"
        return 1
    fi

    # Add header
    sed -i "1iqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" "$output_file"

    log_message "Completed Diamond for $(basename $input_file)"
    return 0
}

# Process R1
log_message "Processing R1..."
run_diamond "$QC1" "$OUT1"
if [ $? -ne 0 ]; then
    log_message "ERROR: Diamond failed for R1"
    exit 1
fi

# Process R2
log_message "Processing R2..."
run_diamond "$QC2" "$OUT2"
if [ $? -ne 0 ]; then
    log_message "ERROR: Diamond failed for R2"
    exit 1
fi

# ============================================================================
# STEP 3: COMBINE AND COMPRESS OUTPUTS
# ============================================================================

log_message "----------------------------------------"
log_message "Step 3: Combining and compressing outputs"
log_message "----------------------------------------"

if command -v pigz &> /dev/null; then
    COMPRESS_CMD="pigz -p 32 -c"
else
    COMPRESS_CMD="gzip -c"
fi

log_message "Compressing outputs..."
$COMPRESS_CMD "$OUT1" > "$OUTBOTH"
if [ $? -ne 0 ]; then
    log_message "ERROR: Failed to compress R1"
    exit 1
fi

tail -n +2 "$OUT2" | $COMPRESS_CMD >> "$OUTBOTH"
if [ $? -ne 0 ]; then
    log_message "ERROR: Failed to append R2"
    exit 1
fi

# Verify output
if [ -s "$OUTBOTH" ]; then
    log_message "Success: $OUTBOTH created"
    rm -f "$OUT1" "$OUT2"

    # Optionally clean up QC files
    if [ "$QC_FILES_CREATED" = true ] && [ "${KEEP_QC_FILES:-false}" != "true" ]; then
        log_message "Cleaning up temporary QC files..."
        rm -f "$QC1" "$QC2"
    fi
else
    log_message "ERROR: Output file is empty"
    exit 1
fi

# Clean up temp directory
rm -rf "$TMP_DIR"

log_message "=========================================="
log_message "Diamond processing completed for: $SAMPLE"
log_message "=========================================="

exit 0
