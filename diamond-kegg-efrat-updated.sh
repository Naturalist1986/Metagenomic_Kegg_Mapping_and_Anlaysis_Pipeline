#!/bin/bash
#SBATCH -o diamond-kegg-trimmomatic-%A_%a.out
#SBATCH --array=0-33
#SBATCH -N 1
#SBATCH -t 7-0
#SBATCH -c 32
#SBATCH --mem=512G
##SBATCH --mail-type=end
##SBATCH --mail-user=moshea@mail.huji.ac.il

# Module and environment setup
# If Trimmomatic is not available as a module, you can download it:
# wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
# unzip Trimmomatic-0.39.zip
# TRIMMOMATIC_JAR="$(pwd)/Trimmomatic-0.39/trimmomatic-0.39.jar"
# === 1. Initialize conda ===
# adjust this path if your miniconda is elsewhere
source ~/miniconda3/etc/profile.d/conda.sh

# === 2. Activate your env ===
conda activate metagenome_assembly

# (optional) sanity checks
echo "Hostname: $(hostname)"
echo "Which diamond: $(which diamond)"
diamond --version

# Set Trimmomatic JAR path (adjust as needed for your system)
TRIMMOMATIC_JAR=${TRIMMOMATIC_JAR:-/sci/home/moshea/miniconda3/envs/metagenome_assembly/share/trimmomatic-0.39-2/trimmomatic.jar}

# Configuration
OFFSET=0
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)

# Define directories
INDIR=/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/raw_fastqs
OUTDIR=/sci/backup/ofinkel/moshea/Efrat_Metagenomes_KEGG_output
QC_DIR=${OUTDIR}/qc_filtered_trimmomatic
LOG_DIR=${OUTDIR}/trimmomatic_logs

# CRITICAL CHANGE: Use local node scratch for temp files to prevent network I/O crashes
# If SLURM_TMPDIR is not defined by your cluster, fallback to /tmp
TMP_DIR=${SLURM_TMPDIR:-/tmp/moshea_diamond_${SLURM_JOB_ID}}

mkdir -p $OUTDIR
mkdir -p $QC_DIR
mkdir -p $TMP_DIR
mkdir -p $LOG_DIR

echo "Input directory: $INDIR"
echo "Output directory: $OUTDIR"
echo "QC directory: $QC_DIR"
echo "Temp directory: $TMP_DIR"
echo "Log directory: $LOG_DIR"
echo "Trimmomatic JAR: $TRIMMOMATIC_JAR"

# Auto-detect all sample files
cd $INDIR
SAMPLE_FILES=($(ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u))

if [ ${#SAMPLE_FILES[@]} -eq 0 ]; then
    echo "ERROR: No *_1.fq.gz files found in $INDIR"
    exit 1
fi

echo "Total samples detected: ${#SAMPLE_FILES[@]}"
echo "Sample list: ${SAMPLE_FILES[@]}"

if [ $LINE_NUM -ge "${#SAMPLE_FILES[@]}" ]; then
    echo "$LINE_NUM exceeds array size (${#SAMPLE_FILES[@]}). Quitting"
    exit 0
fi

sample=${SAMPLE_FILES[$LINE_NUM]}
echo "Processing sample: $sample (index: $LINE_NUM)"

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
    exit 1
fi

# Check if final output already exists
if [ -f "$OUTBOTH" ] && [ -s "$OUTBOTH" ]; then
    echo "Final compressed output already exists: $OUTBOTH"
    exit 0
fi

echo "--------------------------------------------------"
echo $(date)
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
        ADAPTER_FILE="${TMP_DIR}/adapters.fa"
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
        exit 1
    fi
    rm -f $QC1_UNPAIRED $QC2_UNPAIRED
fi

echo "--------------------------------------------------"
echo $(date)
echo "Step 2: DIAMOND mapping to KEGG database (Sequential)"
echo "--------------------------------------------------"

function run_diamond {
    local input_file=$1
    local output_file=$2
    # INCREASED THREADS: Using 30 threads since we are running sequentially
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
        --tmpdir /sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/KEGG_tmp/ \
        --out $output_file \
        --threads $thread_count \
        --outfmt 6 \
        --sensitive \
        --ignore-warnings
    
    local status=$?
    rm -rf "$unique_tmp_dir"
    
    if [ $status -ne 0 ]; then
        echo "Error processing $input_file (exit code: $status)"
        return $status
    fi
    
    sed -i "1iqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" $output_file
    echo "Completed DIAMOND for $input_file"
    return 0
}

# --- SEQUENTIAL EXECUTION BLOCK ---

# Process R1
echo ">>> Processing R1 ($QC1)"
run_diamond $QC1 $OUT1
if [ $? -ne 0 ]; then
    echo "ERROR: DIAMOND failed for R1. Aborting."
    exit 1
fi

# Process R2 (only runs if R1 succeeds)
echo ">>> Processing R2 ($QC2)"
run_diamond $QC2 $OUT2
if [ $? -ne 0 ]; then
    echo "ERROR: DIAMOND failed for R2. Aborting."
    exit 1
fi

echo "--------------------------------------------------"
echo $(date)
echo "Step 3: Compressing and combining outputs"
echo "--------------------------------------------------"

if command -v pigz &> /dev/null; then
    COMPRESS_CMD="pigz -p 32 -c"
else
    COMPRESS_CMD="gzip -c"
fi

echo "Compressing outputs..."
$COMPRESS_CMD $OUT1 > $OUTBOTH || { echo "Error compressing R1"; exit 1; }
tail -n +2 $OUT2 | $COMPRESS_CMD >> $OUTBOTH || { echo "Error appending R2"; exit 1; }

if [ -s "$OUTBOTH" ]; then
    echo "Success: $OUTBOTH created."
    rm -f $OUT1 $OUT2
    if [ "$QC_FILES_CREATED" = true ]; then
        echo "Cleaning up temp QC files..."
        rm -f $QC1 $QC2
    fi
else
    echo "ERROR: Output file is empty."
    exit 1
fi

echo "Job completed successfully at $(date)"