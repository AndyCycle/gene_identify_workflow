#!/bin/bash
# This script processes a single sample for Step 3 of the VEP pipeline.
# It is called by the main vep_project_pipeline.sh script in parallel.

# --- Inherit/Define necessary variables ---
# These variables should be passed as arguments or defined based on standard locations
# For simplicity, we assume they are passed as arguments here.
# Sample ID is the first argument.
SAMPLE_ID="$1"
GENOMES1KG_DIR="$2"
OUTPUT_BASE_DIR="$3"
TEMP_SORT_DIR="$4" # Pass if needed by bcftools sort in inner loops (less likely here)

# --- Exit on error ---
set -e
set -o pipefail

echo "  Processing sample: $SAMPLE_ID - Start time: $(date)"

# --- Configuration specific to this sample ---
SAMPLE_RAW_VCFS_DIR="$OUTPUT_BASE_DIR/200_sample_raw_vcfs" # Needs to be consistent with main script
MERGED_SAMPLE_VCF="$SAMPLE_RAW_VCFS_DIR/${SAMPLE_ID}.1KG.raw.vcf.gz"

# Define the list of chromosomes
CHRS=( {1..22} X )

# Check if the final output for this sample is already done
if [ ! -f "$MERGED_SAMPLE_VCF" ]; then

    echo "    Processing sample: $SAMPLE_ID"
    TEMP_CHR_VCFS=() # Array to store temp VCFs for this sample

    # Extract variants for this sample from each chromosome
    for CHR in "${CHRS[@]}"; do
        # --- Construct the correct filename based on chromosome (duplicate logic from main script) ---
        GENOMES1KG_CHR_VCF=""
        case "$CHR" in
            [1-9]|1[0-9]|2[0-2])
                GENOMES1KG_CHR_VCF="${GENOMES1KG_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
                ;;
            X)
                GENOMES1KG_CHR_VCF="${GENOMES1KG_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
                ;;
            *)
                echo "    WARNING: Unrecognized chromosome format: ${CHR}. Skipping for sample $SAMPLE_ID."
                continue # Skip to next CHR in loop
                ;;
        esac
        # --- End filename construction ---

        if [ -f "$GENOMES1KG_CHR_VCF" ]; then
            TEMP_VCF="$SAMPLE_RAW_VCFS_DIR/${SAMPLE_ID}.chr${CHR}.temp.vcf.gz"
            # Check if temp file exists from a previous partial run and remove (optional)
            # if [ -f "$TEMP_VCF" ]; then rm "$TEMP_VCF" "${TEMP_VCF}.tbi" 2>/dev/null; fi

            echo "      Extracting variants for Chr${CHR}..."
            # -s $SAMPLE_ID : Extract data for this sample
            # -q 0.0001:minor : Keep sites where sample is NOT homozygous reference
            bcftools view -s "$SAMPLE_ID" -q 0.0001:minor "$GENOMES1KG_CHR_VCF" -Oz -o "$TEMP_VCF"
            bcftools index --tbi "$TEMP_VCF" # Index temp file for concat
            TEMP_CHR_VCFS+=("$TEMP_VCF") # Add temp VCF to array
        else
             echo "      WARNING: 1000 Genomes Chr${CHR} VCF not found at expected path: $GENOMES1KG_CHR_VCF. Skipping for sample $SAMPLE_ID."
        fi
    done # End chromosome loop

    # Merge per-chromosome VCFs for this sample
    if [ ${#TEMP_CHR_VCFS[@]} -eq 0 ]; then
        echo "    ERROR: No chromosome VCFs found or processed for sample $SAMPLE_ID. Skipping merge."
        # Optional: Exit with an error code to signal failure for this task
        # exit 1
    else
        echo "    Merging chromosome VCFs for $SAMPLE_ID..."
        bcftools concat "${TEMP_CHR_VCFS[@]}" -a -D -Oz -o "$MERGED_SAMPLE_VCF"

        echo "    Indexing merged VCF for $SAMPLE_ID..."
        bcftools index --tbi "$MERGED_SAMPLE_VCF"

        echo "    Cleaning up temporary chromosome VCFs for $SAMPLE_ID..."
        rm "${TEMP_CHR_VCFS[@]}"
        rm "${TEMP_CHR_VCFS[@]/%.vcf.gz/.vcf.gz.tbi}" 2>/dev/null # Remove indices, suppress errors if index didn't exist

        echo "  Finished processing sample: $SAMPLE_ID. Output: $MERGED_SAMPLE_VCF"
    fi

else
     echo "  Sample $SAMPLE_ID output found: $MERGED_SAMPLE_VCF. Skipping."
fi

echo "  Processing sample: $SAMPLE_ID - End time: $(date)"

# Exit successfully even if some chromosome files were skipped, unless concat failed.
# If concat failed (exit 1 above), the outer set -e will catch it if not commented out.
exit 0