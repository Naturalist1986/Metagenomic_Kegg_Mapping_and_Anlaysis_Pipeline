#!/bin/bash
#SBATCH --array=1-34
#SBATCH -o logs/singlem_%A_%a.out
#SBATCH -n 1 # node count
#SBATCH -c 32
#SBATCH -t 6:00:00
#SBATCH --mem=32G

#conda activate singlem  # Published version 0.19
export SINGLEM_METAPACKAGE_PATH='/sci/backup/aerez/aerez/moshea/singlem_db/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb'

OFFSET=-1
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)

# Updated input directory
indir=/sci/backup/ofinkel/moshea/Noa_Metagenomes_Fastqs/

# Generate array of sample base names from *_1.fq.gz files
# This extracts the prefix before _1.fq.gz for all paired-end samples
Samples=( $(ls $indir/*_1.fq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_1\.fq\.gz$//') )

# Check if array is populated
if [ ${#Samples[@]} -eq 0 ]; then
    echo "No samples found in $indir matching *_1.fq.gz pattern. Quitting"
    exit 1
fi

# For testing with single sample, uncomment:
# Samples=( ces_3_DKDN250001561-1A_22KVLLLT4_L6 )

outdir=/sci/backup/ofinkel/moshea/Noa_Results/Singlem
LINE_NUM=$(echo "$SLURM_ARRAY_TASK_ID + $OFFSET" | bc)

# Check if LINE_NUM exceeds array size
[ $LINE_NUM -ge "${#Samples[@]}" ] && echo "$LINE_NUM exceeds array size. Quitting" && exit 1

s=${Samples[$LINE_NUM]}

# Define input files based on new naming convention
f1=$indir/${s}_1.fq.gz
f2=$indir/${s}_2.fq.gz

# Check if files exist
[ ! -e "$f1" ] && echo "File $f1 does not exist. Quitting" && exit 1
[ ! -e "$f2" ] && echo "File $f2 does not exist. Quitting" && exit 1

# Create output directories
mkdir -p $outdir
mkdir -p $outdir/output
mkdir -p $outdir/otu-table
mkdir -p $outdir/krona
mkdir -p $outdir/archive-otu-table
mkdir -p $outdir/summarize
mkdir -p $outdir/summarize_byfamily
mkdir -p $outdir/summarize_bygenus
mkdir -p $outdir/summarize_byspecies

outfile=$outdir/output/$s.tsv
summarizefile=$outdir/summarize/$s.tsv
summarizefamilyfile=$outdir/summarize_byfamily/$s.tsv
summarizegenusfile=$outdir/summarize_bygenus/$s.tsv
summarizespeciesfile=$outdir/summarize_byspecies/$s.tsv

# Run singlem pipe
if [ -e "$outfile" ]; then
    echo "$outfile exists. Skipping step."
else
    echo "singlem version $(singlem --version)"
    cmd="singlem pipe -1 $f1 -2 $f2 -p $outfile --otu-table $outdir/otu-table/$s-otu_table.tsv --threads 4 --taxonomic-profile-krona $outdir/krona/$s.html --archive-otu-table $outdir/archive-otu-table/$s-archive-otu-table.tsv"
    echo "Running command: $cmd"
    if ! eval "$cmd"; then
        echo "Error in singlem pipe. Exiting."
        exit 1
    fi
fi

# Summarize - full taxonomy
if [ -e "$summarizefile" ]; then
    echo "$summarizefile exists. Skipping step."
else
    cmd="singlem summarise --input-taxonomic-profile $outfile --output-taxonomic-profile-with-extras $summarizefile"
    echo "Running command: $cmd"
    if ! eval "$cmd"; then
        echo "Error in singlem summarise. Exiting."
        exit 1
    fi
fi

# Summarize by family
if [ -e "$summarizefamilyfile" ]; then
    echo "$summarizefamilyfile exists. Skipping step."
else
    cmd="singlem summarise --input-taxonomic-profile $outfile --output-species-by-site-relative-abundance $summarizefamilyfile --output-species-by-site-level family"
    echo "Running command: $cmd"
    if ! eval "$cmd"; then
        echo "Error in singlem summarise (family). Exiting."
        exit 1
    fi
fi

# Summarize by genus
if [ -e "$summarizegenusfile" ]; then
    echo "$summarizegenusfile exists. Skipping step."
else
    cmd="singlem summarise --input-taxonomic-profile $outfile --output-species-by-site-relative-abundance $summarizegenusfile --output-species-by-site-level genus"
    echo "Running command: $cmd"
    if ! eval "$cmd"; then
        echo "Error in singlem summarise (genus). Exiting."
        exit 1
    fi
fi

# Summarize by species
if [ -e "$summarizespeciesfile" ]; then
    echo "$summarizespeciesfile exists. Skipping step."
else
    cmd="singlem summarise --input-taxonomic-profile $outfile --output-species-by-site-relative-abundance $summarizespeciesfile --output-species-by-site-level species"
    echo "Running command: $cmd"
    if ! eval "$cmd"; then
        echo "Error in singlem summarise (species). Exiting."
        exit 1
    fi
fi

echo "---------------------------------------------------"
echo "Finished singlem Successfully."
echo "---------------------------------------------------"
echo "Calculating prok fraction"

prevout=$outdir
outdir="$outdir/prok_fraction"
mkdir -p $outdir
mkdir -p $outdir/per-taxon-read-fraction/
outfile=$outdir/$s.tsv

if [ -e "$outfile" ]; then
    echo "$outfile exists. Skipping step."
else
    cmd="singlem microbial_fraction -1 $f1 -2 $f2 -p $prevout/output/$s.tsv --output-tsv $outfile --output-per-taxon-read-fractions $outdir/per-taxon-read-fraction/${s}-per-taxon-read-fracs.tsv"
    echo "Running command: $cmd"
    if ! eval "$cmd"; then
        echo "Error in singlem microbial_fraction. Exiting."
        exit 1
    fi
fi

echo "---------------------------------------------------"
echo "Finished prok fraction Successfully."
echo "---------------------------------------------------"
