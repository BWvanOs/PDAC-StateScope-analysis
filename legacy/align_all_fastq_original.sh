#!/bin/bash
set -euo pipefail

# Only activate if not already in the correct environment
if [[ "${CONDA_DEFAULT_ENV:-}" != "dnaseq_env" ]]; then
    if [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
        conda activate dnaseq_env
    else
        echo "Could not find conda.sh. Exiting."
        exit 1
    fi
fi

# Define tool paths explicitly to avoid path issues
SAMTOOLS="/LOCATION/OF/SAMTOOLS"
PICARD="picard" 

# Paths
INPUT_DIR="/LOCATION/OF/INPUTFILES"
OUTPUT_DIR="/LOCATION/OF/OUTPUT/"
REF_DIR="$HOME/hg19_reference" ##location of reference
REF_FASTA="$REF_DIR/hg19.fa"

# Step 21: Loop over all R1 files
shopt -s nullglob
R1_FILES=("$INPUT_DIR"/*.R1.fastq.gz)

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "No R1 FASTQ files found in $INPUT_DIR"
    exit 1
fi

for R1 in "${R1_FILES[@]}"; do
    BASENAME=$(basename "$R1" .R1.fastq.gz)
    ##Derive anme of R2 form R1
    R2="$INPUT_DIR/$BASENAME.R2.fastq.gz"
    FINAL_BAM="$OUTPUT_DIR/$BASENAME.final.bam"

    if [ ! -f "$R2" ]; then
        echo "Skipping $BASENAME: R2 file not found"
        continue
    fi

    if [ -f "$FINAL_BAM" ]; then
        echo "Skipping $BASENAME: final BAM already exists"
        continue
    fi

    echo "Processing sample: $BASENAME"

    SORTED_BAM="$OUTPUT_DIR/$BASENAME.sorted.bam"
    RG_BAM="$OUTPUT_DIR/$BASENAME.rg.bam"
    METRICS="$OUTPUT_DIR/$BASENAME.dup_metrics.txt"

    # Align and sort directly (so no SAM file)
    echo "Aligning and sorting..."
    bwa mem -t 16 -M "$REF_FASTA" "$R1" "$R2" | \
    $SAMTOOLS sort -@ 8 -o "$SORTED_BAM"

    # Add read groups
    echo "Adding read groups..."
    $PICARD AddOrReplaceReadGroups \
        --INPUT "$SORTED_BAM" \
        --OUTPUT "$RG_BAM" \
        --RGID 1 \
        --RGLB lib1 \
        --RGPL illumina \
        --RGPU unit1 \
        --RGSM "$BASENAME"
    rm -f "$SORTED_BAM"

    # Mark the duplicates with picard
    echo "Marking duplicates..."
    $PICARD MarkDuplicates \
        --INPUT "$RG_BAM" \
        --OUTPUT "$FINAL_BAM" \
        --METRICS_FILE "$METRICS" \
        --CREATE_INDEX true
    rm -f "$RG_BAM"

    echo "Finished sample: $BASENAME"
done

echo "All samples processed."
