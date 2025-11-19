#!/bin/bash

# ============================================================================
# HELPER SCRIPT FOR DIAMOND-KEGG PIPELINE
# Provides utilities for running the flexible pipeline in different modes
# ============================================================================

PIPELINE_SCRIPT="diamond-kegg-trimmomatic-flexible.sh"
INDIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/raw_fastqs"

# Function to display usage
show_usage() {
    echo "Diamond-KEGG Pipeline Helper"
    echo ""
    echo "Usage: $0 COMMAND [OPTIONS]"
    echo ""
    echo "Commands:"
    echo "  list              List all available samples"
    echo "  run-one SAMPLE    Run pipeline for a single sample"
    echo "  run-list FILE     Run pipeline for samples listed in a file (one per line)"
    echo "  run-range N M     Run pipeline for samples N through M (based on sorted list)"
    echo "  run-pattern PAT   Run pipeline for samples matching pattern (e.g., 'Sample_A*')"
    echo "  run-all           Run pipeline for all available samples"
    echo "  check             Check status of all samples (what's been processed)"
    echo ""
    echo "Examples:"
    echo "  $0 list"
    echo "  $0 run-one Sample_A"
    echo "  $0 run-list my_samples.txt"
    echo "  $0 run-range 0 5"
    echo "  $0 run-pattern 'Treatment_*'"
    echo "  $0 check"
}

# Function to get all available samples
get_all_samples() {
    cd $INDIR
    ls *_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u
}

# Function to check processing status
check_status() {
    echo "Checking processing status..."
    echo ""
    
    local OUTDIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_KEGG_output"
    local all_samples=($(get_all_samples))
    local processed=0
    local unprocessed=0
    
    echo "Sample Status Report:"
    echo "====================="
    
    for sample in "${all_samples[@]}"; do
        if [ -f "${OUTDIR}/${sample}.tsv.gz" ] && [ -s "${OUTDIR}/${sample}.tsv.gz" ]; then
            echo "[✓] $sample - Processed"
            ((processed++))
        else
            echo "[ ] $sample - Not processed"
            ((unprocessed++))
        fi
    done
    
    echo ""
    echo "Summary:"
    echo "  Total samples: ${#all_samples[@]}"
    echo "  Processed: $processed"
    echo "  Unprocessed: $unprocessed"
}

# Main script logic
case "$1" in
    list)
        echo "Available samples in $INDIR:"
        echo "=============================="
        get_all_samples | nl -v 0
        ;;
        
    run-one)
        if [ -z "$2" ]; then
            echo "Error: Please specify a sample name"
            exit 1
        fi
        echo "Submitting job for sample: $2"
        sbatch $PIPELINE_SCRIPT "$2"
        ;;
        
    run-list)
        if [ -z "$2" ] || [ ! -f "$2" ]; then
            echo "Error: Please provide a valid file containing sample names"
            exit 1
        fi
        SAMPLES=$(cat "$2" | grep -v '^#' | grep -v '^$' | tr '\n' ' ')
        echo "Submitting job for $(echo $SAMPLES | wc -w) samples from $2"
        sbatch $PIPELINE_SCRIPT $SAMPLES
        ;;
        
    run-range)
        if [ -z "$2" ] || [ -z "$3" ]; then
            echo "Error: Please specify start and end indices"
            exit 1
        fi
        ALL_SAMPLES=($(get_all_samples))
        SELECTED_SAMPLES=""
        for i in $(seq $2 $3); do
            if [ $i -lt ${#ALL_SAMPLES[@]} ]; then
                SELECTED_SAMPLES="$SELECTED_SAMPLES ${ALL_SAMPLES[$i]}"
            fi
        done
        echo "Submitting job for samples $2 through $3"
        sbatch $PIPELINE_SCRIPT $SELECTED_SAMPLES
        ;;
        
    run-pattern)
        if [ -z "$2" ]; then
            echo "Error: Please specify a pattern"
            exit 1
        fi
        cd $INDIR
        SAMPLES=$(ls $2_1.fq.gz 2>/dev/null | sed 's/_1.fq.gz//g' | sort -u | tr '\n' ' ')
        if [ -z "$SAMPLES" ]; then
            echo "Error: No samples matching pattern '$2'"
            exit 1
        fi
        echo "Submitting job for samples matching pattern '$2'"
        sbatch $PIPELINE_SCRIPT $SAMPLES
        ;;
        
    run-all)
        SAMPLES=$(get_all_samples | tr '\n' ' ')
        echo "Submitting job for ALL samples ($(echo $SAMPLES | wc -w) total)"
        read -p "Are you sure? (y/N): " confirm
        if [ "$confirm" = "y" ] || [ "$confirm" = "Y" ]; then
            sbatch $PIPELINE_SCRIPT $SAMPLES
        else
            echo "Cancelled"
        fi
        ;;
        
    check)
        check_status
        ;;
        
    *)
        show_usage
        exit 1
        ;;
esac