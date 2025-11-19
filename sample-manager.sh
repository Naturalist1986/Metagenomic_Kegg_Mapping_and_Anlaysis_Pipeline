#!/bin/bash

# ============================================================================
# SAMPLE LIST MANAGEMENT UTILITY
# Helps create and manage sample lists for batch processing
# ============================================================================

INDIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/raw_fastqs"
OUTDIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_KEGG_output"

# Function to create sample lists by different criteria
create_sample_list() {
    echo "Sample List Creator"
    echo "==================="
    echo ""
    echo "1. Create list of ALL samples"
    echo "2. Create list of UNPROCESSED samples only"
    echo "3. Create list by pattern matching"
    echo "4. Create list by selecting specific samples"
    echo "5. Create random subset of samples"
    echo "6. Split samples into multiple batch files"
    echo ""
    read -p "Select option (1-6): " option
    
    case $option in
        1)
            # All samples
            OUTPUT_FILE="all_samples.txt"
            cd $INDIR
            ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u > $OUTPUT_FILE
            echo "Created $OUTPUT_FILE with $(wc -l < $OUTPUT_FILE) samples"
            ;;
            
        2)
            # Unprocessed samples only
            OUTPUT_FILE="unprocessed_samples.txt"
            cd $INDIR
            > $OUTPUT_FILE
            for file in *_1.fq.gz; do
                sample=$(echo $file | sed 's/_1.fq.gz//g')
                if [ ! -f "${OUTDIR}/${sample}.tsv.gz" ]; then
                    echo $sample >> $OUTPUT_FILE
                fi
            done
            echo "Created $OUTPUT_FILE with $(wc -l < $OUTPUT_FILE) unprocessed samples"
            ;;
            
        3)
            # Pattern matching
            read -p "Enter pattern (e.g., Treatment_*, Sample_[0-9]*): " pattern
            OUTPUT_FILE="samples_${pattern//\*/ALL}.txt"
            cd $INDIR
            ls ${pattern}_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u > $OUTPUT_FILE
            if [ -s $OUTPUT_FILE ]; then
                echo "Created $OUTPUT_FILE with $(wc -l < $OUTPUT_FILE) samples matching '$pattern'"
            else
                echo "No samples found matching pattern '$pattern'"
                rm -f $OUTPUT_FILE
            fi
            ;;
            
        4)
            # Interactive selection
            OUTPUT_FILE="selected_samples.txt"
            cd $INDIR
            SAMPLES=($(ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u))
            > $OUTPUT_FILE
            
            echo "Select samples (y/n for each, or 'all' to select all remaining):"
            SELECT_ALL=false
            for i in "${!SAMPLES[@]}"; do
                sample=${SAMPLES[$i]}
                if [ "$SELECT_ALL" = true ]; then
                    echo $sample >> $OUTPUT_FILE
                else
                    read -p "[$((i+1))/${#SAMPLES[@]}] Include $sample? (y/n/all): " choice
                    case $choice in
                        y|Y) echo $sample >> $OUTPUT_FILE ;;
                        all|ALL) 
                            SELECT_ALL=true
                            echo $sample >> $OUTPUT_FILE
                            ;;
                    esac
                fi
            done
            echo "Created $OUTPUT_FILE with $(wc -l < $OUTPUT_FILE) selected samples"
            ;;
            
        5)
            # Random subset
            read -p "Enter number of random samples to select: " num
            OUTPUT_FILE="random_${num}_samples.txt"
            cd $INDIR
            ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u | shuf -n $num > $OUTPUT_FILE
            echo "Created $OUTPUT_FILE with $(wc -l < $OUTPUT_FILE) random samples"
            ;;
            
        6)
            # Split into batches
            read -p "Enter number of samples per batch: " batch_size
            PREFIX="batch"
            cd $INDIR
            ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u | \
                split -l $batch_size -d --additional-suffix=.txt - ${PREFIX}_
            
            # Count batch files
            batch_count=$(ls ${PREFIX}_*.txt 2>/dev/null | wc -l)
            echo "Created $batch_count batch files with up to $batch_size samples each:"
            for file in ${PREFIX}_*.txt; do
                echo "  - $file: $(wc -l < $file) samples"
            done
            ;;
            
        *)
            echo "Invalid option"
            return 1
            ;;
    esac
}

# Function to validate sample lists
validate_sample_list() {
    if [ -z "$1" ]; then
        read -p "Enter sample list file to validate: " FILE
    else
        FILE=$1
    fi
    
    if [ ! -f "$FILE" ]; then
        echo "Error: File $FILE not found"
        return 1
    fi
    
    echo "Validating $FILE..."
    echo ""
    
    VALID=0
    INVALID=0
    PROCESSED=0
    UNPROCESSED=0
    
    while read -r sample; do
        # Skip empty lines and comments
        [[ -z "$sample" || "$sample" =~ ^# ]] && continue
        
        # Check if input files exist
        if [ -f "${INDIR}/${sample}_1.fq.gz" ] && [ -f "${INDIR}/${sample}_2.fq.gz" ]; then
            ((VALID++))
            
            # Check if already processed
            if [ -f "${OUTDIR}/${sample}.tsv.gz" ]; then
                echo "[✓] $sample - Valid, already processed"
                ((PROCESSED++))
            else
                echo "[○] $sample - Valid, not yet processed"
                ((UNPROCESSED++))
            fi
        else
            echo "[✗] $sample - INVALID (input files not found)"
            ((INVALID++))
        fi
    done < "$FILE"
    
    echo ""
    echo "Summary:"
    echo "  Total entries: $((VALID + INVALID))"
    echo "  Valid samples: $VALID"
    echo "    - Processed: $PROCESSED"
    echo "    - Unprocessed: $UNPROCESSED"
    echo "  Invalid samples: $INVALID"
}

# Function to generate SLURM array script from sample list
generate_array_script() {
    if [ -z "$1" ]; then
        read -p "Enter sample list file: " FILE
    else
        FILE=$1
    fi
    
    if [ ! -f "$FILE" ]; then
        echo "Error: File $FILE not found"
        return 1
    fi
    
    SAMPLE_COUNT=$(grep -v '^#' "$FILE" | grep -v '^$' | wc -l)
    SCRIPT_NAME="array_script_$(basename $FILE .txt).sh"
    
    cat > $SCRIPT_NAME << 'EOF'
#!/bin/bash
#SBATCH -o diamond-kegg-array-%A_%a.out
EOF
    
    echo "#SBATCH --array=0-$((SAMPLE_COUNT-1))" >> $SCRIPT_NAME
    
    cat >> $SCRIPT_NAME << 'EOF'
#SBATCH -N 1
#SBATCH -t 7-0
#SBATCH -c 32
#SBATCH --mem=512G

# Read sample from list file
EOF
    echo "SAMPLE_LIST=\"$FILE\"" >> $SCRIPT_NAME
    cat >> $SCRIPT_NAME << 'EOF'
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $SAMPLE_LIST | grep -v '^#')

# Run the flexible script with this sample
bash diamond-kegg-trimmomatic-flexible.sh "$SAMPLE"
EOF
    
    echo "Generated SLURM array script: $SCRIPT_NAME"
    echo "  - Processes $SAMPLE_COUNT samples from $FILE"
    echo "  - Submit with: sbatch $SCRIPT_NAME"
}

# Main menu
show_menu() {
    echo ""
    echo "Diamond-KEGG Sample Management Utility"
    echo "======================================"
    echo ""
    echo "1. Create sample list"
    echo "2. Validate sample list"
    echo "3. Generate SLURM array script from list"
    echo "4. Show processing status"
    echo "5. Exit"
    echo ""
}

# Main loop
while true; do
    show_menu
    read -p "Select option (1-5): " choice
    
    case $choice in
        1)
            create_sample_list
            ;;
        2)
            validate_sample_list
            ;;
        3)
            generate_array_script
            ;;
        4)
            # Show processing status
            cd $INDIR
            ALL_SAMPLES=($(ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u))
            PROCESSED=0
            for sample in "${ALL_SAMPLES[@]}"; do
                [ -f "${OUTDIR}/${sample}.tsv.gz" ] && ((PROCESSED++))
            done
            echo ""
            echo "Processing Status:"
            echo "  Total samples: ${#ALL_SAMPLES[@]}"
            echo "  Processed: $PROCESSED"
            echo "  Remaining: $((${#ALL_SAMPLES[@]} - PROCESSED))"
            echo "  Progress: $((PROCESSED * 100 / ${#ALL_SAMPLES[@]}))%"
            ;;
        5)
            echo "Exiting..."
            exit 0
            ;;
        *)
            echo "Invalid option"
            ;;
    esac
done