#!/bin/bash

#SBATCH --account=aerez
#SBATCH --job-name=kegg_sum
#SBATCH --output=/sci/backup/ofinkel/moshea/scripts/logs/job_%A_%a.out
#SBATCH --error=/sci/backup/ofinkel/moshea/scripts/logs/job_%A_%a.err
#SBATCH --array=0-33
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00

# module load spack/all

INPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_post_processing/"
OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_kegg_coverage_summed/"
SCRIPT="/sci/backup/ofinkel/moshea/scripts/sum_kegg_hits.py"

# Create the output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Get the list of files to process
FILES=($(ls ${INPUT_DIR}*.tsv))

# Get the specific file for this array job
FILE_TO_PROCESS=${FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the SRR accession from the file name
SRR_ACCESSION=$(basename "$FILE_TO_PROCESS" | sed 's/-diamond_kegg_hits.tsv//')

# Define the output file path
OUTPUT_FILE=${OUTPUT_DIR}${SRR_ACCESSION}-kegg_hits_summed.tsv

# Skip if the output file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Output file $OUTPUT_FILE already exists. Skipping."
    exit 0
fi

# Run the Python script with a single file
python $SCRIPT --input "$FILE_TO_PROCESS" --output "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "Error processing file $FILE_TO_PROCESS"
    exit 1
fi
