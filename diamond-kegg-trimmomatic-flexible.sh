#!/bin/bash
#SBATCH -o diamond-kegg-trimmomatic-%j.out
#SBATCH -N 1
#SBATCH -t 7-0
#SBATCH -c 32
#SBATCH --mem=512G
##SBATCH --mail-type=end
##SBATCH --mail-user=moshea@mail.huji.ac.il

# ============================================================================
# FLEXIBLE DIAMOND-KEGG PIPELINE
# Usage: 
#   sbatch diamond-kegg-trimmomatic-flexible.sh SAMPLE1 [SAMPLE2 SAMPLE3 ...]
#   OR
#   bash diamond-kegg-trimmomatic-flexible.sh SAMPLE1 [SAMPLE2 SAMPLE3 ...]
# 
# Where SAMPLE is the prefix of your files (without _1.fq.gz or _2.fq.gz)
# Example: 
#   sbatch diamond-kegg-trimmomatic-flexible.sh Sample_A Sample_B
#   This will process Sample_A_1.fq.gz, Sample_A_2.fq.gz, Sample_B_1.fq.gz, Sample_B_2.fq.gz
# ============================================================================

# Check if at least one sample name was provided
if [ $# -eq 0 ]; then
    echo "Error: No sample names provided!"
    echo "Usage: $0 SAMPLE1 [SAMPLE2 SAMPLE3 ...]"
    echo "Example: $0 Sample_A Sample_B"
    echo "Note: Provide sample prefixes without _1.fq.gz or _2.fq.gz extensions"
    exit 1
fi

# === 1. Initialize conda ===
# adjust this path if your miniconda is elsewhere
source ~/miniconda3/etc/profile.d/conda.sh

# === 2. Activate your env ===
conda activate metagenome_assembly

# (optional) sanity checks
echo "=============================================="
echo "Job started at: $(date)"
echo "Hostname: $(hostname)"
echo "Which diamond: $(which diamond)"
diamond --version
echo "Number of samples to process: $#"
echo "Sample list: $@"
echo "=============================================="

# Set Trimmomatic JAR path (adjust as needed for your system)
TRIMMOMATIC_JAR=${TRIMMOMATIC_JAR:-/sci/home/moshea/miniconda3/envs/metagenome_assembly/share/trimmomatic-0.39-2/trimmomatic.jar}

# Define directories
INDIR=/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/raw_fastqs
OUTDIR=/sci/backup/ofinkel/moshea/Efrat_Metagenomes_KEGG_output
QC_DIR=${OUTDIR}/qc_filtered_trimmomatic
LOG_DIR=${OUTDIR}/trimmomatic_logs
DIAMOND_TMP=/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/KEGG_tmp

# CRITICAL: Use local node scratch for temp files to prevent network I/O crashes
TMP_DIR=${SLURM_TMPDIR:-/tmp/moshea_diamond_${SLURM_JOB_ID:-$$}}

# Create necessary directories
mkdir -p $OUTDIR
mkdir -p $QC_DIR
mkdir -p $TMP_DIR
mkdir -p $LOG_DIR
mkdir -p $DIAMOND_TMP

echo "Input directory: $INDIR"
echo "Output directory: $OUTDIR"
echo "QC directory: $QC_DIR"
echo "Temp directory: $TMP_DIR"
echo "Log directory: $LOG_DIR"
echo "Trimmomatic JAR: $TRIMMOMATIC_JAR"

# Function to process a single sample
process_sample() {
    local sample=$1
    echo ""
    echo "=============================================="
    echo "Processing sample: $sample"
    echo "Start time: $(date)"
    echo "=============================================="
    
    # Define input and output files
    IN1=${INDIR}/${sample}_1.fq.gz
    IN2=${INDIR}/${sample}_2.fq.gz
    QC1=${QC_DIR}/${sample}_R1_paired.fq.gz
    QC2=${QC_DIR}/${sample}_R2_paired.fq.gz
    QC1_UNPAIRED=${QC_DIR}/${sample}_R1_unpaired.fq.gz
    QC2_UNPAIRED=${QC_DIR}/${sample}_R2_unpaired.fq.gz
    OUT1=${OUTDIR}/${sample}_R1.tsv
    OUT2=${OUTDIR}/${sample}_R2.tsv
    OUTBOTH=${OUTDIR}/${sample}.tsv.gz
    TRIM_LOG=${LOG_DIR}/${sample}_trimmomatic.log
    
    # Verify input files exist
    if [ ! -f "$IN1" ] || [ ! -f "$IN2" ]; then
        echo "ERROR: Input files for $sample do not exist!"
        echo "  Expected: $IN1"
        echo "  Expected: $IN2"
        return 1
    fi
    
    # Check if final output already exists
    if [ -f "$OUTBOTH" ] && [ -s "$OUTBOTH" ]; then
        echo "Final compressed output already exists: $OUTBOTH"
        echo "Skipping sample $sample"
        return 0
    fi
    
    echo "--------------------------------------------------"
    echo "Step 1: Quality filtering with Trimmomatic"
    echo "--------------------------------------------------"
    
    QC_FILES_CREATED=false
    if [ -f "$QC1" ] && [ -f "$QC2" ] && [ -s "$QC1" ] && [ -s "$QC2" ]; then
        echo "QC files already exist. Skipping Trimmomatic."
    else
        echo "QC files not found. Running Trimmomatic..."
        QC_FILES_CREATED=true
        
        ADAPTER_DIR=$(dirname $TRIMMOMATIC_JAR)/adapters
        if [ -d "$ADAPTER_DIR" ]; then
            ADAPTER_FILE="$ADAPTER_DIR/TruSeq3-PE-2.fa"
        else
            ADAPTER_FILE="${TMP_DIR}/adapters_${sample}.fa"
            cat > $ADAPTER_FILE << 'EOF'
>PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
EOF
        fi
        
        java -jar $TRIMMOMATIC_JAR PE \
            -threads 32 \
            -phred33 \
            $IN1 $IN2 \
            $QC1 $QC1_UNPAIRED \
            $QC2 $QC2_UNPAIRED \
            ILLUMINACLIP:${ADAPTER_FILE}:2:30:10:2:keepBothReads \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:20 \
            MINLEN:50 \
            2> $TRIM_LOG
        
        if [ $? -ne 0 ]; then
            echo "ERROR: Trimmomatic failed for sample $sample"
            return 1
        fi
        rm -f $QC1_UNPAIRED $QC2_UNPAIRED
    fi
    
    echo "--------------------------------------------------"
    echo "Step 2: DIAMOND mapping to KEGG database"
    echo "--------------------------------------------------"
    
    # Function to run DIAMOND on a single file
    run_diamond() {
        local input_file=$1
        local output_file=$2
        local thread_count=30
        local sample_name=$(basename $input_file .fq)
        
        if [ -f "$output_file" ]; then
            rm -f "$output_file"
        fi
        
        echo "Starting DIAMOND on $input_file with $thread_count threads..."
        
        diamond blastx \
            --query $input_file \
            --db /sci/backup/aerez/aerez/moshea/kegg_2023/prokaryotes.kegg.dmnd \
            --max-target-seqs 1 \
            --evalue 1e-5 \
            --verbose \
            --tmpdir $DIAMOND_TMP \
            --out $output_file \
            --threads $thread_count \
            --outfmt 6 \
            --sensitive \
            --ignore-warnings
        
        local status=$?
        
        if [ $status -ne 0 ]; then
            echo "Error processing $input_file (exit code: $status)"
            return $status
        fi
        
        # Add header
        sed -i "1iqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" $output_file
        echo "Completed DIAMOND for $input_file"
        return 0
    }
    
    # Process R1
    echo ">>> Processing R1 ($QC1)"
    run_diamond $QC1 $OUT1
    if [ $? -ne 0 ]; then
        echo "ERROR: DIAMOND failed for R1. Aborting sample $sample."
        return 1
    fi
    
    # Process R2
    echo ">>> Processing R2 ($QC2)"
    run_diamond $QC2 $OUT2
    if [ $? -ne 0 ]; then
        echo "ERROR: DIAMOND failed for R2. Aborting sample $sample."
        return 1
    fi
    
    echo "--------------------------------------------------"
    echo "Step 3: Compressing and combining outputs"
    echo "--------------------------------------------------"
    
    if command -v pigz &> /dev/null; then
        COMPRESS_CMD="pigz -p 32 -c"
    else
        COMPRESS_CMD="gzip -c"
    fi
    
    echo "Compressing outputs..."
    $COMPRESS_CMD $OUT1 > $OUTBOTH || { echo "Error compressing R1"; return 1; }
    tail -n +2 $OUT2 | $COMPRESS_CMD >> $OUTBOTH || { echo "Error appending R2"; return 1; }
    
    if [ -s "$OUTBOTH" ]; then
        echo "Success: $OUTBOTH created."
        rm -f $OUT1 $OUT2
        if [ "$QC_FILES_CREATED" = true ]; then
            echo "Cleaning up temp QC files..."
            rm -f $QC1 $QC2
        fi
    else
        echo "ERROR: Output file is empty."
        return 1
    fi
    
    echo "Sample $sample completed successfully at $(date)"
    return 0
}

# Main processing loop
FAILED_SAMPLES=()
SUCCESSFUL_SAMPLES=()

for sample in "$@"; do
    if process_sample "$sample"; then
        SUCCESSFUL_SAMPLES+=("$sample")
    else
        FAILED_SAMPLES+=("$sample")
        echo "WARNING: Sample $sample failed. Continuing with next sample..."
    fi
done

# Clean up temp directory
rm -rf $TMP_DIR

# Final summary
echo ""
echo "=============================================="
echo "FINAL SUMMARY"
echo "=============================================="
echo "Job completed at: $(date)"
echo "Total samples processed: $#"
echo "Successful: ${#SUCCESSFUL_SAMPLES[@]}"
echo "Failed: ${#FAILED_SAMPLES[@]}"

if [ ${#SUCCESSFUL_SAMPLES[@]} -gt 0 ]; then
    echo ""
    echo "Successfully processed samples:"
    for s in "${SUCCESSFUL_SAMPLES[@]}"; do
        echo "  - $s"
    done
fi

if [ ${#FAILED_SAMPLES[@]} -gt 0 ]; then
    echo ""
    echo "Failed samples:"
    for s in "${FAILED_SAMPLES[@]}"; do
        echo "  - $s"
    done
    exit 1
fi

echo ""
echo "All samples processed successfully!"
exit 0