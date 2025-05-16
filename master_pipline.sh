#!/bin/bash
#SBATCH --job-name=vep_project_pipeline
#SBATCH --output=%j_pipeline.log
#SBATCH --partition=x86_64_8cpu  # Adjust based on your cluster
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8       # Max CPUs for any single step, VEP uses --fork
#SBATCH --time=72:00:00         # Adjust based on expected total runtime

export PATH=/share/home/ggwsxysz1lsy/miniconda3/bin:$PATH

# 激活指定环境
source activate cyc_vep

# --- Exit on error ---
set -e
set -o pipefail

echo "Pipeline Start time: $(date)"

# --- Configuration ---
OMIM_FEATURES="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/omim/omim_features.txt"
OMIM_FEATURES_TEST="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/omim/omim_features_T.txt"
CLINVAR_VCF="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/ClinVar/vcf_GRCh37/clinvar_20250409.vcf.gz"
GENOMES1KG_DIR="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/1000Genomes/hgdownload.soe.ucsc.edu/gbdb/hg19/1000Genomes/phase3"
SAMPLE_IDS_FILE="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/1000Genomes/random_200_samples.txt"

# Define the directory containing per-chromosome 1000G VCFs (adjust if needed)
GENOMES1KG_CHR_VCF_PATTERN="${GENOMES1KG_DIR}/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

OUTPUT_BASE_DIR="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/annotated_datasets"
mkdir -p "$OUTPUT_BASE_DIR"

TEMP_SORT_DIR="$OUTPUT_BASE_DIR/temp_sort"
mkdir -p "$TEMP_SORT_DIR" # Ensure this directory exists and is writable


# --- VEP Cache and Reference FASTA ---
VEP_CACHE_DIR="/share/home/ggwsxysz1lsy/.vep"
REF_FASTA="/share/home/ggwsxysz1lsy/.vep/homo_sapiens/105_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

# --- Check dependencies ---
if ! command -v bcftools &> /dev/null; then echo "ERROR: bcftools not found"; exit 1; fi
if ! command -v vep &> /dev/null; then echo "ERROR: vep not found"; exit 1; fi
if ! command -v bgzip &> /dev/null; then echo "ERROR: bgzip not found"; exit 1; fi
if ! command -v tabix &> /dev/null; then echo "ERROR: tabix not found"; exit 1; fi # VEP output might need tabix indexing

# --- Check input files ---
if [ ! -f "$OMIM_FEATURES" ]; then echo "ERROR: OMIM features file not found: $OMIM_FEATURES"; exit 1; fi
if [ ! -f "$OMIM_FEATURES_TEST" ]; then echo "ERROR: OMIM test features file not found: $OMIM_FEATURES_TEST"; exit 1; fi
if [ ! -f "$CLINVAR_VCF" ]; then echo "ERROR: ClinVar VCF not found: $CLINVAR_VCF"; exit 1; fi
if [ ! -d "$GENOMES1KG_DIR" ]; then echo "ERROR: 1000 Genomes directory not found: $GENOMES1KG_DIR"; exit 1; fi
if [ ! -f "$SAMPLE_IDS_FILE" ]; then echo "ERROR: Sample IDs file not found: $SAMPLE_IDS_FILE"; exit 1; fi
if [ ! -d "$VEP_CACHE_DIR" ]; then echo "ERROR: VEP Cache directory not found: $VEP_CACHE_DIR"; exit 1; fi
if [ ! -f "$REF_FASTA" ]; then echo "ERROR: Reference FASTA not found: $REF_FASTA"; exit 1; fi
# Check FASTA index, create if missing
if [ ! -f "${REF_FASTA}.fai" ]; then echo "WARNING: Reference FASTA index (.fai) not found for $REF_FASTA. Attempting to create..."; samtools faidx "$REF_FASTA" || { echo "ERROR: Failed to index FASTA"; exit 1; } fi


# --- Define Standard VEP Annotation Parameters (Now Comprehensive) ---
VEP_ANNOTATION_PARAMS=(
    --format vcf --vcf --compress_output bgzip --force_overwrite --cache --dir_cache "$VEP_CACHE_DIR"
    --assembly GRCh37 --fasta "$REF_FASTA" --offline --fork ${SLURM_CPUS_PER_TASK:-4}
    --everything # !!! Enable everything
    --verbose # Keep verbose for debugging

    # --- Add plugins for comprehensive predictive and other data ---
    # Ensure these plugins are INSTALLED using vep_install
    # --plugin dbNSFP
    # --plugin SpliceAI
    # --plugin LoFtool # If needed
    # Add other plugins here as needed
)


# --- STEP 1: Build Set_P_Train ---
echo "--- Starting Step 1: Building Set_P_Train ---"
SET_P_TRAIN_CLINVAR_FILTERED="$OUTPUT_BASE_DIR/Set_P_Train_ClinVar_filtered.vcf.gz"
SET_P_TRAIN_OMIM_MATCHED="$OUTPUT_BASE_DIR/Set_P_Train_OMIM_matched.vcf.gz" # Output of custom script
SET_P_TRAIN_VEP_VCF="$OUTPUT_BASE_DIR/Set_P_Train.vep.vcf.gz"
SET_P_TRAIN_VEP_STATS="$OUTPUT_BASE_DIR/Set_P_Train.vep_stats.html"

if [ ! -f "$SET_P_TRAIN_VEP_VCF" ]; then # Check if final output exists

    echo "Step 1.1: Filtering ClinVar for Pathogenic/Likely_pathogenic variants with OMIM association..."
    if [ ! -f "$SET_P_TRAIN_CLINVAR_FILTERED" ]; then
        bcftools view "$CLINVAR_VCF" \
            -i '(INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic") && INFO/CLNDISDB~"OMIM:"' \
            -Oz -o "$SET_P_TRAIN_CLINVAR_FILTERED"
        bcftools index --tbi "$SET_P_TRAIN_CLINVAR_FILTERED"
        echo "Step 1.1 Complete: $SET_P_TRAIN_CLINVAR_FILTERED"
    else
        echo "Step 1.1 output found: $SET_P_TRAIN_CLINVAR_FILTERED. Skipping."
    fi

    echo "Step 1.2: Matching ClinVar variants to OMIM IDs in omim_features.txt (Requires Custom Script)"
    if [ ! -f "$SET_P_TRAIN_OMIM_MATCHED" ]; then
        # --- !!! CUSTOM SCRIPT CALL !!! ---
        YOUR_CUSTOM_SCRIPT_PATH="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/pipline" # <--- REPLACE WITH ACTUAL PATH
        PYTHON_SCRIPT="$YOUR_CUSTOM_SCRIPT_PATH/filter_clinvar_by_omim.py"
        if [ ! -f "$PYTHON_SCRIPT" ]; then echo "ERROR: Custom script not found: $PYTHON_SCRIPT"; exit 1; fi

        echo "Running custom script to match OMIM IDs..."
        python "$PYTHON_SCRIPT" \
            --vcf "$SET_P_TRAIN_CLINVAR_FILTERED" \
            --omim "$OMIM_FEATURES" \
            --output "$SET_P_TRAIN_OMIM_MATCHED" # Python script now outputs UNSORTED VCF
        # --- !!! END CUSTOM SCRIPT CALL !!! ---

        # Check if the unsorted VCF was created
        if [ ! -f "$SET_P_TRAIN_OMIM_MATCHED" ]; then
            echo "ERROR: Custom script for Step 1.2 failed or did not create output VCF: $SET_P_TRAIN_OMIM_MATCHED"
            exit 1
        fi
        echo "Step 1.2 Complete: $SET_P_TRAIN_OMIM_MATCHED (Unsorted)"

        echo "Step 1.2.1: Sorting the filtered VCF..."
        SET_P_TRAIN_OMIM_MATCHED_SORTED="${SET_P_TRAIN_OMIM_MATCHED}.sorted"
        if [ ! -f "$SET_P_TRAIN_OMIM_MATCHED_SORTED" ]; then
            # Use bcftools sort to sort by coordinate
            bcftools sort "$SET_P_TRAIN_OMIM_MATCHED" -Oz -o "$SET_P_TRAIN_OMIM_MATCHED_SORTED" --temp-dir "$TEMP_SORT_DIR" # Use a temp dir with enough space
            echo "Step 1.2.1 Complete: $SET_P_TRAIN_OMIM_MATCHED_SORTED"
        else
            echo "Step 1.2.1 output found: $SET_P_TRAIN_OMIM_MATCHED_SORTED. Skipping sort."
        fi

        echo "Step 1.2.2: Indexing the sorted VCF..."
        if [ ! -f "$SET_P_TRAIN_OMIM_MATCHED_SORTED.tbi" ]; then
            bcftools index --tbi "$SET_P_TRAIN_OMIM_MATCHED_SORTED"
            echo "Step 1.2.2 Complete: $SET_P_TRAIN_OMIM_MATCHED_SORTED.tbi"
        else
            echo "Step 1.2.2 output found: $SET_P_TRAIN_OMIM_MATCHED_SORTED.tbi. Skipping index."
        fi

        # Replace the original unsorted file with the sorted one for next steps
        echo "Replacing unsorted VCF with sorted VCF..."
        mv "$SET_P_TRAIN_OMIM_MATCHED_SORTED" "$SET_P_TRAIN_OMIM_MATCHED"
        mv "$SET_P_TRAIN_OMIM_MATCHED_SORTED.tbi" "$SET_P_TRAIN_OMIM_MATCHED.tbi"
        echo "Unsorted VCF replaced by sorted VCF and its index."

        # Now $SET_P_TRAIN_OMIM_MATCHED is the sorted and indexed VCF for Step 1.3

    else # Skip Step 1.2 if final output already exists (now implies sorted and indexed)
        echo "Step 1.2 final output found: $SET_P_TRAIN_OMIM_MATCHED. Skipping Steps 1.2, 1.2.1, 1.2.2."
        # Ensure index exists if skipping
        if [ ! -f "$SET_P_TRAIN_OMIM_MATCHED.tbi" ]; then
            echo "WARNING: Final output VCF found but index missing. Attempting to index $SET_P_TRAIN_OMIM_MATCHED..."
            bcftools index --tbi "$SET_P_TRAIN_OMIM_MATCHED" || { echo "ERROR: Failed to index existing file."; exit 1; }
        fi
    fi

    echo "Step 1.3: VEP Annotation of Set_P_Train variants..."
    # Use the output from the custom script as input for VEP
    if [ ! -f "$SET_P_TRAIN_OMIM_MATCHED" ]; then
        echo "ERROR: Input for VEP annotation (Step 1.3) not found: $SET_P_TRAIN_OMIM_MATCHED"
        exit 1
    fi
    vep --input_file "$SET_P_TRAIN_OMIM_MATCHED" \
        --output_file "$SET_P_TRAIN_VEP_VCF" \
        --stats_file "$SET_P_TRAIN_VEP_STATS" \
        "${VEP_ANNOTATION_PARAMS[@]}" # Expand the predefined parameters

    bcftools index --tbi "$SET_P_TRAIN_VEP_VCF"
    echo "Step 1.3 Complete: $SET_P_TRAIN_VEP_VCF"

    # Optional: Clean up intermediate file from step 1.1/1.2 if successful
    # rm "$SET_P_TRAIN_CLINVAR_FILTERED" "${SET_P_TRAIN_CLINVAR_FILTERED}.tbi"
    # rm "$SET_P_TRAIN_OMIM_MATCHED" "${SET_P_TRAIN_OMIM_MATCHED}.tbi"


else
    echo "Step 1 final output found: $SET_P_TRAIN_VEP_VCF. Skipping Step 1."
fi

echo "--- Finished Step 1 ---"


# --- STEP 2: Build Set_NP_Train ---
echo "--- Starting Step 2: Building Set_NP_Train ---"
CLINVAR_BENIGN_FILTERED="$OUTPUT_BASE_DIR/ClinVar_Benign_LikelyBenign.vcf.gz"
# We need a VCF of all sites from 1000G Phase 3 without sample info.
# This is resource intensive for full genome. Let's do it per chromosome and concatenate.
# WARNING: Creating a single VCF of ALL sites from 1000G Phase 3 is very large and time consuming.
# The intersection with ClinVar Benign might also be slow.
# Consider alternative Set_NP definition (e.g., high frequency 1000G sites not in ClinVar P/LP).
GENOMES1KG_ALL_SITES_VCF="$OUTPUT_BASE_DIR/1000Genomes_Phase3_AllSites.vcf.gz"
SET_NP_TRAIN_RAW_VCF="$OUTPUT_BASE_DIR/Set_NP_Train_raw.vcf.gz" # Intersection output
SET_NP_TRAIN_VEP_VCF="$OUTPUT_BASE_DIR/Set_NP_Train.vep.vcf.gz"
SET_NP_TRAIN_VEP_STATS="$OUTPUT_BASE_DIR/Set_NP_Train.vep_stats.html"

if [ ! -f "$SET_NP_TRAIN_VEP_VCF" ]; then # Check if final output exists

    echo "Step 2.1: Filtering ClinVar for Benign/Likely_benign variants..."
    if [ ! -f "$CLINVAR_BENIGN_FILTERED" ]; then
        bcftools view "$CLINVAR_VCF" \
            -i 'INFO/CLNSIG~"Benign" || INFO/CLNSIG~"Likely_benign"' \
            -Oz -o "$CLINVAR_BENIGN_FILTERED"
        bcftools index --tbi "$CLINVAR_BENIGN_FILTERED"
        echo "Step 2.1 Complete: $CLINVAR_BENIGN_FILTERED"
    else
        echo "Step 2.1 output found: $CLINVAR_BENIGN_FILTERED. Skipping."
    fi

    echo "Step 2.2: Creating a VCF of all sites from 1000 Genomes Phase 3 (without sample info)..."
    if [ ! -f "$GENOMES1KG_ALL_SITES_VCF" ]; then
        echo "  Processing chromosomes..."
        TEMP_ALL_SITES_VCFS=()
        # Define chromosomes in correct order, including X and Y
        CHRS=( {1..22} X Y )

        for CHR in "${CHRS[@]}"; do
            # --- Construct the correct filename based on chromosome ---
            GENOMES1KG_CHR_VCF=""
            case "$CHR" in
                [1-9]|1[0-9]|2[0-2])
                    GENOMES1KG_CHR_VCF="${GENOMES1KG_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
                    ;;
                X)
                    GENOMES1KG_CHR_VCF="${GENOMES1KG_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
                    ;;
                Y)
                    GENOMES1KG_CHR_VCF="${GENOMES1KG_DIR}/ALL.chr${CHR}.phase3_integrated_v1b.20130502.genotypes.vcf.gz"
                    ;;
                *)
                    echo "    WARNING: Unrecognized chromosome format: ${CHR}. Skipping."
                    continue # Skip to next CHR in loop
                    ;;
            esac
            # --- End filename construction ---

            # Define the temporary site VCF name for this chromosome
            temp_vcf="$OUTPUT_BASE_DIR/ALL.chr${CHR}.phase3.sites.vcf.gz"

            if [ -f "$GENOMES1KG_CHR_VCF" ]; then
                # --- Check if the temporary site VCF already exists ---
                if [ ! -f "$temp_vcf" ] || [ ! -f "$temp_vcf.tbi" ]; then
                    echo "    Processing $GENOMES1KG_CHR_VCF -> $temp_vcf"
                    # -G: Drop genotypes, keep only site information
                    bcftools view -G "$GENOMES1KG_CHR_VCF" -Oz -o "$temp_vcf"
                    # Add indexing for concatenation
                    bcftools index --tbi "$temp_vcf"
                    echo "    Processed Chr${CHR}. Output: $temp_vcf"
                else
                    echo "    Temporary site VCF found for Chr${CHR}: $temp_vcf. Skipping processing."
                fi
                TEMP_ALL_SITES_VCFS+=("$temp_vcf") # Add temp VCF to array regardless of whether it was just created or pre-existing
            else
                 echo "    WARNING: 1000 Genomes Chr${CHR} VCF not found at expected path: $GENOMES1KG_CHR_VCF. Skipping."
            fi
        done # End chromosome loop

        if [ ${#TEMP_ALL_SITES_VCFS[@]} -eq 0 ]; then
             echo "ERROR: No 1000 Genomes chromosome site VCFs found or processed to create all sites VCF."
             exit 1
        fi

        echo "  Concatenating chromosome site VCFs in order..."
        # bcftools concat takes multiple input files directly.
        # Remove --temp-dir here as concat doesn't support it.
        # Ensure you use the array directly: "${TEMP_ALL_SITES_VCFS[@]}"
        bcftools concat "${TEMP_ALL_SITES_VCFS[@]}" -a -D -Oz -o "$GENOMES1KG_ALL_SITES_VCF"

        echo "  Indexing concatenated VCF..."
        bcftools index --tbi "$GENOMES1KG_ALL_SITES_VCF"

        echo "  Cleaning up temporary site VCFs and indices..."
        rm "${TEMP_ALL_SITES_VCFS[@]}"
        rm "${TEMP_ALL_SITES_VCFS[@]/%.vcf.gz/.vcf.gz.tbi}"

        echo "Step 2.2 Complete: $GENOMES1KG_ALL_SITES_VCF"
    else
        echo "Step 2.2 output found: $GENOMES1KG_ALL_SITES_VCF. Skipping."
    fi


    echo "Step 2.3: Finding intersection of 1000G sites and ClinVar Benign/LikelyBenign..."
    if [ ! -f "$SET_NP_TRAIN_RAW_VCF" ] || [ ! -f "$SET_NP_TRAIN_RAW_VCF.tbi" ]; then
        if [ ! -f "$GENOMES1KG_ALL_SITES_VCF" ] || [ ! -f "$GENOMES1KG_ALL_SITES_VCF.tbi" ]; then
             echo "ERROR: Input 1000G All Sites VCF for Step 2.3 not found or indexed: $GENOMES1KG_ALL_SITES_VCF. Ensure Step 2.2 completed."
             exit 1
        fi
        if [ ! -f "$CLINVAR_BENIGN_FILTERED" ] || [ ! -f "$CLINVAR_BENIGN_FILTERED.tbi" ]; then
             echo "ERROR: Input ClinVar Benign VCF for Step 2.3 not found or indexed: $CLINVAR_BENIGN_FILTERED. Ensure Step 2.1 completed."
             exit 1
        fi

        # bcftools isec outputs multiple files based on sets. We want variants present in BOTH (-n+=2).
        # Use --temp-dir for isec!
        ISECT_NP_TEMP_DIR="$OUTPUT_BASE_DIR/isec_np_temp_$$" # Use PID to make temp dir unique
        mkdir -p "$ISECT_NP_TEMP_DIR" # Ensure temp dir is created before use

        echo "  Running bcftools isec..."
        # Use -n +2 for "present in at least 2 files". For 2 input files, this is the intersection.
        bcftools isec -n +2 "$GENOMES1KG_ALL_SITES_VCF" "$CLINVAR_BENIGN_FILTERED" \
             -p "$ISECT_NP_TEMP_DIR" -Oz

        # Check if the intersection file was created. It should be named "0000.vcf.gz" in the output dir.
        INTERSECTION_FILE="$ISECT_NP_TEMP_DIR/0000.vcf.gz"
        if [ ! -f "$INTERSECTION_FILE" ]; then
            echo "ERROR: bcftools isec failed or did not create the intersection file: $INTERSECTION_FILE"
            rm -r "$ISECT_NP_TEMP_DIR" # Clean up temp dir on failure
            exit 1
        fi

        echo "  Moving intersection file to final location..."
        mv "$INTERSECTION_FILE" "$SET_NP_TRAIN_RAW_VCF"
        # Check if isec also created the index, otherwise create it.
        if [ -f "$INTERSECTION_FILE.tbi" ]; then
             mv "$INTERSECTION_FILE.tbi" "$SET_NP_TRAIN_RAW_VCF.tbi"
        else
             echo "  Index not created by isec, creating index for $SET_NP_TRAIN_RAW_VCF..."
             bcftools index --tbi "$SET_NP_TRAIN_RAW_VCF"
        fi

        echo "  Cleaning up isec temporary directory..."
        rm -r "$ISECT_NP_TEMP_DIR"

        echo "Step 2.3 Complete: $SET_NP_TRAIN_RAW_VCF"
    else
        echo "Step 2.3 output found: $SET_NP_TRAIN_RAW_VCF. Skipping."
    fi


    echo "Step 2.4: VEP Annotation of Set_NP_Train variants..."
     if [ ! -f "$SET_NP_TRAIN_RAW_VCF" ]; then
        echo "ERROR: Input for VEP annotation (Step 2.4) not found: $SET_NP_TRAIN_RAW_VCF"
        exit 1
    fi
    vep --input_file "$SET_NP_TRAIN_RAW_VCF" \
        --output_file "$SET_NP_TRAIN_VEP_VCF" \
        --stats_file "$SET_NP_TRAIN_VEP_STATS" \
        "${VEP_ANNOTATION_PARAMS[@]}" # Expand the predefined parameters

    bcftools index --tbi "$SET_NP_TRAIN_VEP_VCF"
    echo "Step 2.4 Complete: $SET_NP_TRAIN_VEP_VCF"

    # Optional: Clean up intermediate files from step 2.1, 2.2, 2.3 if successful
    # rm "$CLINVAR_BENIGN_FILTERED" "${CLINVAR_BENIGN_FILTERED}.tbi"
    # rm "$GENOMES1KG_ALL_SITES_VCF" "${GENOMES1KG_ALL_SITES_VCF}.tbi"
    # Intermediate files from isec are removed above

else
    echo "Step 2 final output found: $SET_NP_TRAIN_VEP_VCF. Skipping Step 2."
fi
echo "--- Finished Step 2 ---"


# --- STEP 3: Extract and Merge 200 Sample Raw Variant VCFs (Parallelized) ---
echo "--- Starting Step 3: Extracting and Merging 200 Sample Raw Variant VCFs ---"
SAMPLE_RAW_VCFS_DIR="$OUTPUT_BASE_DIR/200_sample_raw_vcfs"
mkdir -p "$SAMPLE_RAW_VCFS_DIR"

# Define the path to the new script
PROCESS_SINGLE_SAMPLE_SCRIPT="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/pipline/process_single_sample_step3.sh" # <--- Set this path correctly!

# Check if the single sample processing script exists and is executable
if [ ! -f "$PROCESS_SINGLE_SAMPLE_SCRIPT" ]; then echo "ERROR: Single sample processing script not found: $PROCESS_SINGLE_SAMPLE_SCRIPT"; exit 1; fi
if [ ! -x "$PROCESS_SINGLE_SAMPLE_SCRIPT" ]; then echo "ERROR: Single sample processing script is not executable: $PROCESS_SINGLE_SAMPLE_SCRIPT"; exit 1; fi


# Determine the number of parallel jobs to run. Use allocated CPUs.
# You might want to use fewer than total CPUs if processes are memory intensive or I/O bound.
# Let's use the allocated number of CPUs as the default.
NUM_PARALLEL_JOBS=${SLURM_CPUS_PER_TASK:-4} # Default to 4 if not run under Slurm

echo "Step 3: Processing samples from $SAMPLE_IDS_FILE in parallel ($NUM_PARALLEL_JOBS jobs)"

# Use xargs to run the single sample script in parallel
# cat "$SAMPLE_IDS_FILE" reads the sample IDs from the file
# xargs -I {} : replace {} with each line read from stdin
# -P $NUM_PARALLEL_JOBS : run up to NUM_PARALLEL_JOBS processes in parallel
# Pass necessary variables from the main script as arguments to the single sample script
# The arguments passed are: $SAMPLE_ID, $GENOMES1KG_DIR, $OUTPUT_BASE_DIR, $TEMP_SORT_DIR (if needed by single script)
cat "$SAMPLE_IDS_FILE" | xargs -I {} -P "$NUM_PARALLEL_JOBS" \
    "$PROCESS_SINGLE_SAMPLE_SCRIPT" {} "$GENOMES1KG_DIR" "$OUTPUT_BASE_DIR" "$TEMP_SORT_DIR"

# Note: xargs exits with status 0 if all commands it invokes exit with 0.
# If any invocation of process_single_sample_step3.sh exits non-zero, xargs will usually exit non-zero.
# set -e should catch this.

echo "Step 3 Complete: Finished extracting and merging raw VCFs for 200 samples (parallelized)."

echo "--- Finished Step 3 ---"



# --- STEP 4: Build Set_P_Test ---
echo "--- Starting Step 4: Building Set_P_Test ---"
SET_P_TEST_CLINVAR_FILTERED="$OUTPUT_BASE_DIR/Set_P_Test_ClinVar_filtered.vcf.gz"
# SET_P_TEST_OMIM_MATCHED="$OUTPUT_BASE_DIR/Set_P_Test_OMIM_matched.vcf.gz" # This variable seems unused if python script outputs _temp
SET_P_TEST_RAW_VCF="$OUTPUT_BASE_DIR/Set_P_Test_raw.vcf.gz" # Final output with flag

if [ ! -f "$SET_P_TEST_RAW_VCF" ]; then # Check if final output of Step 4 exists

    echo "Step 4.1: Filtering ClinVar for Pathogenic/Likely_pathogenic variants with OMIM association (Test Set)..."
    if [ ! -f "$SET_P_TEST_CLINVAR_FILTERED" ]; then
        bcftools view "$CLINVAR_VCF" \
            -i '(INFO/CLNSIG~"Pathogenic" || INFO/CLNSIG~"Likely_pathogenic") && INFO/CLNDISDB~"OMIM:"' \
            -Oz -o "$SET_P_TEST_CLINVAR_FILTERED"
        bcftools index --tbi "$SET_P_TEST_CLINVAR_FILTERED"
        echo "Step 4.1 Complete: $SET_P_TEST_CLINVAR_FILTERED"
    else
        echo "Step 4.1 output found: $SET_P_TEST_CLINVAR_FILTERED. Skipping."
    fi # Closes if for Step 4.1

    echo "Step 4.2.1: Filtering ClinVar test variants by OMIM IDs (Requires Custom Script)"
    SET_P_TEST_OMIM_MATCHED_TEMP="$OUTPUT_BASE_DIR/Set_P_Test_OMIM_matched_temp.vcf.gz" # Temp output before flag

    if [ ! -f "$SET_P_TEST_OMIM_MATCHED_TEMP" ]; then
        YOUR_CUSTOM_SCRIPT_PATH="/share/home/ggwsxysz1lsy/USER/chenyanchao/rare_gene_identification/datasets/pipline" # <--- REPLACE WITH ACTUAL PATH
        PYTHON_SCRIPT="$YOUR_CUSTOM_SCRIPT_PATH/filter_clinvar_by_omim.py"
        if [ ! -f "$PYTHON_SCRIPT" ]; then echo "ERROR: Custom script not found: $PYTHON_SCRIPT"; exit 1; fi

        echo "Running custom script to match OMIM test IDs..."
        python "$PYTHON_SCRIPT" \
            --vcf "$SET_P_TEST_CLINVAR_FILTERED" \
            --omim "$OMIM_FEATURES_TEST" \
            --output "$SET_P_TEST_OMIM_MATCHED_TEMP" # Python script outputs UNSORTED TEMP VCF
        echo "Step 4.2.1 Complete: $SET_P_TEST_OMIM_MATCHED_TEMP (Unsorted)"
        if [ ! -f "$SET_P_TEST_OMIM_MATCHED_TEMP" ]; then
            echo "ERROR: Custom script for Step 4.2.1 failed or did not create output VCF: $SET_P_TEST_OMIM_MATCHED_TEMP"
            exit 1
        fi
    else
        echo "Step 4.2.1 output found: $SET_P_TEST_OMIM_MATCHED_TEMP. Skipping filtering."
    fi # Closes if for Step 4.2.1


    echo "Step 4.2.2: Sorting the filtered test VCF..."
    SET_P_TEST_OMIM_MATCHED_SORTED="${SET_P_TEST_OMIM_MATCHED_TEMP}.sorted"
    if [ ! -f "$SET_P_TEST_OMIM_MATCHED_SORTED" ]; then
        bcftools sort "$SET_P_TEST_OMIM_MATCHED_TEMP" -Oz -o "$SET_P_TEST_OMIM_MATCHED_SORTED" --temp-dir  "$TEMP_SORT_DIR" # Use a temp dir
        echo "Step 4.2.2 Complete: $SET_P_TEST_OMIM_MATCHED_SORTED"
    else
        echo "Step 4.2.2 output found: $SET_P_TEST_OMIM_MATCHED_SORTED. Skipping sort."
    fi # Closes if for Step 4.2.2
    # Clean up unsorted temp file
    # rm "$SET_P_TEST_OMIM_MATCHED_TEMP" "${SET_P_TEST_OMIM_MATCHED_TEMP}.tbi" 2>/dev/null


    echo "Step 4.3: Adding INSERTED_TEST_VAR flag to sorted Set_P_Test variants (using bcftools 1.21 compatible method)..."

    # Define the new INFO field header line string
    INSERTED_FLAG_HEADER_LINE='##INFO=<ID=INSERTED_TEST_VAR,Number=0,Type=Flag,Description="Variant inserted for testing (from Set_P_Test)">'

    # Create a temporary file containing the sites (CHROM and POS) from the sorted VCF
    TEMP_SITES_FILE_UNSORTED="$OUTPUT_BASE_DIR/Set_P_Test_sites_temp_unsorted_$$" # Unsorted temp file
    TEMP_SITES_FILE_SORTED="$OUTPUT_BASE_DIR/Set_P_Test_sites_temp_sorted_$$"     # Sorted temp file
    echo "  Extracting sites to tag to temporary file: $TEMP_SITES_FILE_UNSORTED"

    # Extract CHROM and POS from the sorted VCF (excluding header lines)
    bcftools view "$SET_P_TEST_OMIM_MATCHED_SORTED" | grep -v "^#" | awk '{print $1"\t"$2}' > "$TEMP_SITES_FILE_UNSORTED"

    # Check if the temporary sites file was created and has content
    if [ ! -s "$TEMP_SITES_FILE_UNSORTED" ]; then
        echo "ERROR: Failed to create temporary sites file for tagging, or file is empty: $TEMP_SITES_FILE_UNSORTED"
        # rm "$TEMP_SITES_FILE_UNSORTED" # Optional cleanup
        exit 1
    fi
    echo "  Temporary unsorted sites file created successfully."

    # --- Add Sorting and Indexing steps for the sites file ---
    echo "  Sorting temporary sites file..."
    # Sort the file numerically by chromosome (column 1) and position (column 2)
    # Use -k1,1V -k2,2n for variant call format (VCF) sorting order if chromosomes are mixed strings/numbers.
    # For simple CHROM POS, basic sort -k1,1 -k2,2n might suffice, but VCF sort is safer.
    # bcftools sort is for VCF/BCF, so we'll use a standard sort command.
    sort -k1,1V -k2,2n "$TEMP_SITES_FILE_UNSORTED" > "$TEMP_SITES_FILE_SORTED"

    # Check if sorting was successful
    if [ ! -s "$TEMP_SITES_FILE_SORTED" ]; then
        echo "ERROR: Sorting of temporary sites file failed or resulted in empty file: $TEMP_SITES_FILE_SORTED"
        rm "$TEMP_SITES_FILE_UNSORTED" # Clean up unsorted temp file
        exit 1
    fi
    echo "  Temporary sorted sites file created successfully."

    echo "  Indexing temporary sorted sites file..."
    # tabix requires the file to be bgzip compressed and sorted.
    # bgzip the sorted sites file.
    bgzip "$TEMP_SITES_FILE_SORTED"
    # The bgzip command replaces the original file with a .gz version. Update the variable.
    TEMP_SITES_FILE_SORTED_GZ="$TEMP_SITES_FILE_SORTED.gz"
    if [ ! -f "$TEMP_SITES_FILE_SORTED_GZ" ]; then
        echo "ERROR: bgzip compression of temporary sites file failed: $TEMP_SITES_FILE_SORTED_GZ"
        rm "$TEMP_SITES_FILE_UNSORTED" "$TEMP_SITES_FILE_SORTED" # Clean up temp files
        exit 1
    fi
    echo "  Temporary sorted sites file compressed with bgzip."

    # Index the bgzip compressed file using tabix.
    # -p vcf tells tabix to use VCF format indexing rules (even though it's not a full VCF, it has CHROM/POS)
    # Depending on tabix version and specific format, -p bed might also work if the file has 3 columns (CHROM START END).
    # For CHROM POS, -p vcf is often the correct type for coordinate sorting.
    tabix -p vcf "$TEMP_SITES_FILE_SORTED_GZ"

    # Check if indexing was successful
    if [ ! -f "$TEMP_SITES_FILE_SORTED_GZ.tbi" ]; then
        echo "ERROR: tabix indexing of temporary sites file failed: $TEMP_SITES_FILE_SORTED_GZ.tbi"
        rm "$TEMP_SITES_FILE_UNSORTED" "$TEMP_SITES_FILE_SORTED" "$TEMP_SITES_FILE_SORTED_GZ" # Clean up temp files
        exit 1
    fi
    echo "  Temporary sorted sites file indexed successfully."
    # --- End Sorting and Indexing steps ---


    echo "  Adding flag using bcftools annotate -m..."
    # Use bcftools annotate -a sites_file -c CHROM,POS -m +TAG -H header_line input.vcf -o output.vcf
    # Now use the bgzip compressed and indexed file as the annotation file: "$TEMP_SITES_FILE_SORTED_GZ"
    bcftools annotate -a "$TEMP_SITES_FILE_SORTED_GZ" -c CHROM,POS \
        -m +INSERTED_TEST_VAR \
        -H "$INSERTED_FLAG_HEADER_LINE" \
        "$SET_P_TEST_OMIM_MATCHED_SORTED" \
        -Oz -o "$SET_P_TEST_RAW_VCF"
        # --temp-dir is not supported in bcftools 1.21 annotate

    if [ ! -f "$SET_P_TEST_RAW_VCF" ]; then
        echo "ERROR: bcftools annotate failed to create output VCF: $SET_P_TEST_RAW_VCF"
        rm "$TEMP_SITES_FILE_UNSORTED" "$TEMP_SITES_FILE_SORTED" "$TEMP_SITES_FILE_SORTED_GZ" "$TEMP_SITES_FILE_SORTED_GZ.tbi" 2>/dev/null # Clean up temp files
        exit 1
    fi

    echo "  Cleaning up temporary sites files and indices..."
    rm "$TEMP_SITES_FILE_UNSORTED" "$TEMP_SITES_FILE_SORTED" "$TEMP_SITES_FILE_SORTED_GZ" "$TEMP_SITES_FILE_SORTED_GZ.tbi" 2>/dev/null

    echo "  Indexing the final output VCF..."
    bcftools index --tbi "$SET_P_TEST_RAW_VCF"
    echo "Step 4.3 Complete: $SET_P_TEST_RAW_VCF (Sorted and flagged)"

    # Optional: Clean up intermediate sorted VCF if Step 4.3 was successful
    # rm "$SET_P_TEST_OMIM_MATCHED_SORTED" "${SET_P_TEST_OMIM_MATCHED_SORTED}.tbi" 2>/dev/null

    # Optional: Clean up intermediate files from step 4.1, 4.2.1, 4.2.2 if successful
    # rm "$SET_P_TEST_CLINVAR_FILTERED" "${SET_P_TEST_CLINVAR_FILTERED}.tbi" 2>/dev/null
    # rm "$SET_P_TEST_OMIM_MATCHED_TEMP" "${SET_P_TEST_OMIM_MATCHED_TEMP}.tbi" 2>/dev/null
    # rm "$SET_P_TEST_OMIM_MATCHED_SORTED" "${SET_P_TEST_OMIM_MATCHED_SORTED}.tbi" 2>/dev/null

else # Corresponds to the outer if [ ! -f "$SET_P_TEST_RAW_VCF" ]
    echo "Step 4 final output found: $SET_P_TEST_RAW_VCF. Skipping Step 4."
    # Ensure index exists if skipping
    if [ ! -f "$SET_P_TEST_RAW_VCF.tbi" ]; then
        echo "WARNING: Final output VCF found but index missing. Attempting to index $SET_P_TEST_RAW_VCF..."
        bcftools index --tbi "$SET_P_TEST_RAW_VCF" || { echo "ERROR: Failed to index existing file."; exit 1; }
    fi
fi # Closes the outer if for Step 4
echo "--- Finished Step 4 ---"


# --- STEP 5: Build OKG-T Dataset (Merge Set_P_Train and Set_NP_Train) ---
echo "--- Starting Step 5: Building OKG-T Dataset ---"
OKG_T_VCF="$OUTPUT_BASE_DIR/OKG_T.vep.vcf.gz"

if [ ! -f "$OKG_T_VCF" ]; then # Check if final output exists
    echo "Step 5: Merging Set_P_Train and Set_NP_Train to build OKG-T..."
     if [ ! -f "$SET_P_TRAIN_VEP_VCF" ] || [ ! -f "$SET_NP_TRAIN_VEP_VCF" ]; then
        echo "ERROR: Input VEP annotated VCFs for Step 5 not found. Ensure Step 1 and 2 completed."
        exit 1
    fi
    # Use bcftools merge. Set_P_Train and Set_NP_Train have no sample columns.
    # -m none: keeps all alleles for multiallelic sites if they overlap
    # --force-samples: allows merge even if sample lists don't match (they should both be empty)
    bcftools merge "$SET_P_TRAIN_VEP_VCF" "$SET_NP_TRAIN_VEP_VCF" \
        -m none -Oz -o "$OKG_T_VCF" \
        --force-samples

    bcftools index --tbi "$OKG_T_VCF"
    echo "Step 5 Complete: $OKG_T_VCF"
    echo "Number of variants in OKG-T:"
    bcftools index -n "$OKG_T_VCF"
else
    echo "Step 5 final output found: $OKG_T_VCF. Skipping Step 5."
fi
echo "--- Finished Step 5 ---"


# --- STEP 6: Build OKG-I Dataset (Insert Test Variants, VEP Annotate, Filter) ---
echo "--- Starting Step 6: Building OKG-I Dataset ---"
OKG_I_DIR="$OUTPUT_BASE_DIR/OKG_I_sample_vcfs"
mkdir -p "$OKG_I_DIR"

# Use the raw sample VCFs generated in Step 3 and the Set_P_Test VCF from Step 4
SAMPLE_RAW_VCFS_DIR="$OUTPUT_BASE_DIR/200_sample_raw_vcfs"
SET_P_TEST_RAW_VCF="$OUTPUT_BASE_DIR/Set_P_Test_raw.vcf.gz" # Contains INSERTED_TEST_VAR flag

# Check if input files for Step 6 are ready
if [ ! -d "$SAMPLE_RAW_VCFS_DIR" ] || [ ! -f "$SET_P_TEST_RAW_VCF" ]; then
     echo "ERROR: Input directories/files for Step 6 not found. Ensure Step 3 and 4 completed."
     exit 1
fi
# Check if at least one sample raw VCF exists
if [ -z "$(find "$SAMPLE_RAW_VCFS_DIR" -name "*.1KG.raw.vcf.gz" -print -quit)" ]; then
    echo "ERROR: No raw sample VCFs found in $SAMPLE_RAW_VCFS_DIR. Ensure Step 3 completed successfully."
    exit 1
fi


# --- Loop through each sample raw VCF ---
# Iterate based on the sample IDs file, not directory contents, for consistency
while read SAMPLE_ID; do
    SAMPLE_RAW_VCF="$SAMPLE_RAW_VCFS_DIR/${SAMPLE_ID}.1KG.raw.vcf.gz"
    MERGED_VCF="$OKG_I_DIR/${SAMPLE_ID}.merged_test_vars.vcf.gz"
    ANNOTATED_VCF="$OKG_I_DIR/${SAMPLE_ID}.vep_annotated.vcf.gz"
    SPLIT_VEP_VCF="$OKG_I_DIR/${SAMPLE_ID}.vep_split.vcf.gz" # Intermediate for easier filtering
    FILTERED_VCF="$OKG_I_DIR/${SAMPLE_ID}.filtered_candidate.vcf.gz"
    VEP_STATS="$OKG_I_DIR/${SAMPLE_ID}.vep_stats.html"

    # Check if the final output for this sample is already done
    if [ ! -f "$FILTERED_VCF" ]; then

        echo "Processing sample $SAMPLE_ID for OKG-I..."

        # 6.1 Merge sample raw VCF with Set_P_Test raw VCF
        echo "  Step 6.1: Merging raw sample variants with Set_P_Test variants..."
        if [ ! -f "$MERGED_VCF" ]; then
            if [ ! -f "$SAMPLE_RAW_VCF" ]; then
                 echo "    ERROR: Raw sample VCF not found for $SAMPLE_ID: $SAMPLE_RAW_VCF. Skipping sample."
                 continue # Skip to next sample ID
            fi
            bcftools merge "$SAMPLE_RAW_VCF" "$SET_P_TEST_RAW_VCF" -Oz -o "$MERGED_VCF" --force-samples -m none
            bcftools index --tbi "$MERGED_VCF"
            echo "  Step 6.1 Complete: $MERGED_VCF"
        else
             echo "  Step 6.1 output found: $MERGED_VCF. Skipping."
        fi

        # 6.2 Run VEP Annotation on the merged VCF
        echo "  Step 6.2: Running VEP annotation on merged VCF..."
         if [ ! -f "$ANNOTATED_VCF" ]; then
            if [ ! -f "$MERGED_VCF" ]; then
                 echo "    ERROR: Input for VEP annotation (Step 6.2) not found: $MERGED_VCF. Skipping sample."
                 continue
            fi
            # VEP fork parameter here is per VEP run. Adjust as needed.
            vep --input_file "$MERGED_VCF" \
                --output_file "$ANNOTATED_VCF" \
                --stats_file "$VEP_STATS" \
                "${VEP_ANNOTATION_PARAMS[@]}"

            bcftools index --tbi "$ANNOTATED_VCF"
            echo "  Step 6.2 Complete: $ANNOTATED_VCF"
        else
             echo "  Step 6.2 output found: $ANNOTATED_VCF. Skipping."
        fi

        # 6.3 Split VEP CSQ field for easier filtering
        echo "  Step 6.3: Running bcftools +split-vep..."
        if [ ! -f "$SPLIT_VEP_VCF" ]; then
            if [ ! -f "$ANNOTATED_VCF" ]; then
                 echo "    ERROR: Input for split-vep (Step 6.3) not found: $ANNOTATED_VCF. Skipping sample."
                 continue
            fi
            # Define the fields you expect in the CSQ field from your VEP annotation settings.
            # CHECK A SAMPLE ANNOTATED VCF HEADER to get the exact list after running Step 6.2 successfully!
            CSQ_FIELDS='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|HGVS_OFFSET|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|gnomAD_AF' # Removed LoFtool, add other fields if VEP adds them without explicit plugins
            bcftools +split-vep -i "$ANNOTATED_VCF" -o "$SPLIT_VEP_VCF" -O z -f "$CSQ_FIELDS" -c CSQ -s worst
            bcftools index --tbi "$SPLIT_VEP_VCF"
            echo "  Step 6.3 Complete: $SPLIT_VEP_VCF"
        else
            echo "  Step 6.3 output found: $SPLIT_VEP_VCF. Skipping."
        fi

        # 6.4 Filter down to ~400 candidate variants
        echo "  Step 6.4: Filtering for ~400 candidate variants..."
        if [ ! -f "$FILTERED_VCF" ]; then
             if [ ! -f "$SPLIT_VEP_VCF" ]; then
                 echo "    ERROR: Input for filtering (Step 6.4) not found: $SPLIT_VEP_VCF. Skipping sample."
                 continue
            fi
            # --- DEFINE YOUR FILTER_EXPRESSION HERE ---
            # This expression must rely ONLY on fields available in the split VEP output (CSQ fields + INSERTED_TEST_VAR).
            # It needs careful testing and adjustment of thresholds to achieve ~400 variants per sample on average.
            # Using standard VEP fields: VEP_IMPACT, VEP_CLIN_SIG, VEP_SIFT, VEP_PolyPhen, VEP_gnomAD_AF.
            # VEP_SIFT and VEP_PolyPhen will only have 'b' (basic) scores if plugins aren't used.
            # The logic is: Keep inserted variants OR keep rare (low gnomAD_AF) variants that are potentially impactful based on standard VEP annotation.
            FILTER_EXPRESSION='( INFO/INSERTED_TEST_VAR = 1 ) || \
                               ( VEP_gnomAD_AF < 0.001 && \
                                 ( VEP_IMPACT = "HIGH" || VEP_IMPACT = "MODERATE" || \
                                   VEP_CLIN_SIG ~ "Pathogenic" || VEP_CLIN_SIG ~ "Likely_pathogenic" || VEP_CLIN_SIG ~ "VUS" || \
                                   VEP_SIFT ~ "deleterious" || VEP_PolyPhen ~ "damaging" \
                                 ) \
                               )'
            # Adjust VEP_gnomAD_AF threshold (e.g., 0.001, 0.0005, 0.0001) and potentially include more VEP fields
            # based on what standard VEP outputs without plugins to tune the ~400 variants count.
            # You might need to run Step 6.1-6.3 for a few samples, check the split output format and counts,
            # then adjust this FILTER_EXPRESSION and re-run Step 6 for all samples.

            bcftools filter -i "$FILTER_EXPRESSION" -o "$FILTERED_VCF" -O z "$SPLIT_VEP_VCF"

            bcftools index --tbi "$FILTERED_VCF"

            # Report variant count for checking
            VARIANT_COUNT=$(bcftools index -n "$FILTERED_VCF" 2>/dev/null || echo 0) # Handle case where index fails on empty file
            echo "  Sample $SAMPLE_ID filtered variants: $VARIANT_COUNT"

            echo "  Step 6.4 Complete: $FILTERED_VCF"
        else
            echo "  Step 6.4 output found: $FILTERED_VCF. Skipping."
        fi


        # Optional: Clean up intermediate files for this sample if successful
        # rm "$MERGED_VCF" "${MERGED_VCF}.tbi"
        # rm "$ANNOTATED_VCF" "${ANNOTATED_VCF}.tbi"
        # rm "$SPLIT_VEP_VCF" "${SPLIT_VEP_VCF}.tbi"

        echo "Finished processing sample $SAMPLE_ID for OKG-I."

    else
        echo "Sample $SAMPLE_ID final output found: $FILTERED_VCF. Skipping processing for this sample."
    fi

done < "$SAMPLE_IDS_FILE" # Read sample IDs from file

echo "--- Finished Step 6 ---"

# --- POST-PROCESSING STEPS (Require Custom Scripts or external tools like InterVar) ---
echo "--- Starting Post-Processing Steps ---"

# --- Post-processing 1: Add Disease Descriptions to VEP Output VCFs ---
# This step applies to OKG-T_VCF and potentially the final filtered OKG-I VCFs
# You need a custom script to add disease info from omim_features.txt/omim_features_T.txt
# This involves parsing VEP output (gene symbols, OMIM IDs from Existing_variation etc.)
# and matching it to your OMIM feature files.
echo "Post-processing 1: Adding disease descriptions (Requires Custom Script)"
OKG_T_WITH_DISEASE="$OUTPUT_BASE_DIR/OKG_T.vep.disease.vcf.gz"
if [ ! -f "$OKG_T_WITH_DISEASE" ]; then
    echo "Running custom script to add disease info to OKG-T. Command placeholder:"
    echo "YOUR_CUSTOM_SCRIPT_PATH/add_disease_description.py --vcf \"$OKG_T_VCF\" --omim \"$OMIM_FEATURES\" -o \"$OKG_T_WITH_DISEASE\""
    # Execute your custom script here:
    # YOUR_CUSTOM_SCRIPT_PATH/add_disease_description.py --vcf "$OKG_T_VCF" --omim "$OMIM_FEATURES" -o "$OKG_T_WITH_DISEASE"
    # Example: Simulate script completion (REMOVE THIS IN REAL SCRIPT)
    touch "$OKG_T_WITH_DISEASE"
     if [ ! -f "$OKG_T_WITH_DISEASE" ]; then echo "ERROR: Custom script for adding disease info failed for OKG-T."; exit 1; fi
    echo "Post-processing 1 for OKG-T Complete: $OKG_T_WITH_DISEASE"
else
    echo "Post-processing 1 for OKG-T output found: $OKG_T_WITH_DISEASE. Skipping."
fi

# You might also want to add disease descriptions to the OKG-I sample VCFs.
# This would involve looping through the filtered OKG-I VCFs ($FILTERED_VCF for each sample)
# And running your custom script on each of them.
echo "NOTE: Adding disease descriptions to OKG-I sample VCFs is also needed and requires looping through each filtered sample VCF."


# --- Post-processing 2: Prepare Input for InterVar (ACMG Classification) ---
# This step applies to the final annotated VCFs (OKG-T_VCF, or OKG-T_WITH_DISEASE if you created it,
# and the final filtered OKG-I VCFs after adding disease info).
# InterVar requires specific input formats. You need a custom script or use InterVar's pre-processing tools.
# Due to lack of many predictive scores without plugins, InterVar results might be limited.
echo "Post-processing 2: Preparing input for InterVar (Requires Custom Script/Tool)"
OKG_T_INTERVAR_INPUT="$OUTPUT_BASE_DIR/OKG_T.intervar_input.txt"
if [ ! -f "$OKG_T_INTERVAR_INPUT" ]; then
    echo "Running custom script/tool to prepare InterVar input for OKG-T. Command placeholder:"
    # Replace OKG_T_VCF with OKG_T_WITH_DISEASE if you added disease info
    echo "YOUR_CUSTOM_SCRIPT_PATH/prepare_intervar_input.py --vcf \"$OKG_T_VCF\" -o \"$OKG_T_INTERVAR_INPUT\""
    # Execute your custom script/tool here:
    # YOUR_CUSTOM_SCRIPT_PATH/prepare_intervar_input.py --vcf "$OKG_T_VCF" -o "$OKG_T_INTERVAR_INPUT"
     # Example: Simulate script completion (REMOVE THIS IN REAL SCRIPT)
    touch "$OKG_T_INTERVAR_INPUT"
    if [ ! -f "$OKG_T_INTERVAR_INPUT" ]; then echo "ERROR: Custom script/tool for InterVar input failed for OKG-T."; exit 1; fi
    echo "Post-processing 2 for OKG-T Complete: $OKG_T_INTERVAR_INPUT"
else
    echo "Post-processing 2 for OKG-T output found: $OKG_T_INTERVAR_INPUT. Skipping."
fi

# You also need to prepare InterVar input for each filtered OKG-I sample VCF.
# This requires looping through each sample's final VCF and running your preparation script/tool.
echo "NOTE: Preparing InterVar input for OKG-I sample VCFs is also needed and requires looping through each filtered sample VCF."

# --- Running InterVar ---
# After preparing input files, you will need to run InterVar itself.
# Consult InterVar documentation for command usage.
echo "NOTE: Running InterVar for ACMG classification is the next step after preparing input files."


echo "--- Finished Post-Processing Steps ---"

echo "Pipeline End time: $(date)"