#!/usr/bin/env python3

import pysam
import argparse
import os
import sys

def main():
    """
    Filters ClinVar VCF variants based on whether their associated OMIM IDs
    in the CLNDISDB field are present in a provided OMIM features file.
    """
    parser = argparse.ArgumentParser(
        description='Filter ClinVar VCF variants based on OMIM IDs in a features file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Show default values
    )
    parser.add_argument('--vcf', required=True,
                        help='Input ClinVar VCF file (gzipped and indexed). Assumed to have CLNDISDB INFO field.')
    parser.add_argument('--omim', required=True,
                        help='OMIM features text file (tab-separated, OMIM ID in col 1).')
    parser.add_argument('--output', required=True,
                        help='Output filtered VCF file (will be gzipped and indexed).')
    args = parser.parse_args()

    print(f"[{os.path.basename(__file__)}] Filtering VCF '{args.vcf}' using OMIM features from '{args.omim}'")

    # --- Step 1: Load allowed OMIM IDs from the features file ---
    allowed_omim_ids = set()
    try:
        with open(args.omim, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: # Skip empty lines
                    continue
                # Split only on the first tab, as the second column might contain spaces or tabs
                parts = line.split('\t', 1)
                if len(parts) > 0:
                    omim_id = parts[0].strip()
                    if omim_id: # Ensure the extracted ID is not empty
                        allowed_omim_ids.add(omim_id)
    except FileNotFoundError:
        print(f"[{os.path.basename(__file__)}] ERROR: OMIM features file not found: {args.omim}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[{os.path.basename(__file__)}] ERROR reading OMIM features file {args.omim}: {e}", file=sys.stderr)
        sys.exit(1)

    if not allowed_omim_ids:
        print(f"[{os.path.basename(__file__)}] WARNING: No OMIM IDs loaded from {args.omim}. Output VCF will be empty.", file=sys.stderr)
    else:
        print(f"[{os.path.basename(__file__)}] Loaded {len(allowed_omim_ids)} unique OMIM IDs from {args.omim}")

    # --- Step 2: Process the VCF ---
    try:
        # Open input VCF (pysam automatically handles .gz and indexing)
        # Need to make sure .tbi or .csi index file exists alongside the .vcf.gz
        try:
            vcf_in = pysam.VariantFile(args.vcf, 'r')
        except pysam.utils.TbifileError as e:
             print(f"[{os.path.basename(__file__)}] ERROR: Could not open VCF index for {args.vcf}. Please ensure .tbi or .csi index exists.", file=sys.stderr)
             print(f"[{os.path.basename(__file__)}] Error details: {e}", file=sys.stderr)
             sys.exit(1)
        except FileNotFoundError:
             print(f"[{os.path.basename(__file__)}] ERROR: Input VCF file not found: {args.vcf}", file=sys.stderr)
             sys.exit(1)


        # Open output VCF with input header. 'w' mode overwrites existing file.
        vcf_out = pysam.VariantFile(args.output, 'w', header=vcf_in.header)

        filtered_count = 0
        total_count = 0

        # Iterate through records
        print(f"[{os.path.basename(__file__)}] Processing VCF records...")
        for record in vcf_in:
            total_count += 1
            found_match = False

            # Check if CLNDISDB INFO field exists
            if 'CLNDISDB' in record.info:
                clndisdb_value = record.info['CLNDISDB'] # This can be a string or list depending on VCF header

                # Ensure we iterate over the value(s)
                if not isinstance(clndisdb_value, (list, tuple)):
                     clndisdb_values = [clndisdb_value] # Wrap in list if it's a single string
                else:
                     clndisdb_values = clndisdb_value # It's already a list/tuple


                for clndisdb_entry_str in clndisdb_values:
                    if not isinstance(clndisdb_entry_str, str):
                         # Should not happen with standard VCFs, but safety check
                         print(f"[{os.path.basename(__file__)}] WARNING: Unexpected type in CLNDISDB for record {record.chrom}:{record.pos}: {type(clndisdb_entry_str)}. Skipping entry.", file=sys.stderr)
                         continue

                    # CLNDISDB is a comma-separated string of Source:ID pairs
                    source_id_pairs = clndisdb_entry_str.split(',')
                    for pair in source_id_pairs:
                         pair = pair.strip() # Remove leading/trailing whitespace
                         if pair.lower().startswith('omim:'): # Case-insensitive check for "OMIM:"
                             try:
                                 # Split only on the first colon to get the ID
                                 omim_id = pair.split(':', 1)[1].strip()
                                 # Check if the extracted OMIM ID is in our allowed set
                                 if omim_id and omim_id in allowed_omim_ids:
                                     found_match = True
                                     # No need to check other CLNDISDB entries for this record if one matches
                                     break
                                 elif not omim_id:
                                     # Handle cases like "OMIM:" with no ID
                                     print(f"[{os.path.basename(__file__)}] WARNING: Empty OMIM ID found in CLNDISDB for {record.chrom}:{record.pos} - '{pair}'. Skipping.", file=sys.sys.stderr)

                             except IndexError:
                                 # Handle malformed OMIM: entry that doesn't have an ID part
                                 print(f"[{os.path.basename(__file__)}] WARNING: Malformed OMIM entry in CLNDISDB for {record.chrom}:{record.pos} - '{pair}'. Skipping.", file=sys.stderr)
                         # else: # This entry is not an OMIM link, ignore

                    if found_match:
                        # If a match was found from any pair, we can stop checking pairs for this CLNDISDB value
                        break

                # If a match was found from any CLNDISDB value, we can stop checking CLNDISDB values for this record
                if found_match:
                    vcf_out.write(record)
                    filtered_count += 1

            # Optional: print progress indicator
            # if total_count % 10000 == 0:
            #    print(f"[{os.path.basename(__file__)}] Processed {total_count} records...")


        print(f"[{os.path.basename(__file__)}] Finished processing records.")
        print(f"[{os.path.basename(__file__)}] Processed {total_count} variants. Wrote {filtered_count} filtered variants to '{args.output}'.")

        # Close files (optional, done automatically on exit or garbage collection, but good practice)
        vcf_in.close()
        vcf_out.close()

        # --- Step 3: Index the output VCF ---
        # Indexing will be handled by the calling shell script after sorting.
        print(f"[{os.path.basename(__file__)}] Filtering complete. Output VCF written to '{args.output}'. Sorting and indexing required next.")

    except Exception as e:
        print(f"[{os.path.basename(__file__)}] An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)



if __name__ == "__main__":
    main()