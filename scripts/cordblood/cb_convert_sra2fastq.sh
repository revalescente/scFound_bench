#!/bin/bash

# Exit immediately if a command exits with a non-zero status (fails)
set -e

# Navigate to your base directory
cd ~/cite_seq_coord_blood_raw_data

# Loop through every directory that starts with "SRR"
for SRR_DIR in SRR*/; do
    
    # Remove the trailing slash to get the exact ID (e.g., "SRR5808750")
    SRR_ID=${SRR_DIR%/}

    echo "====================================================="
    echo "Processing $SRR_ID..."
    echo "====================================================="

    # Go into the specific SRR folder
    cd "$SRR_ID"

    # Check if the .sra file actually exists in this folder
    if [ -f "${SRR_ID}.sra" ]; then
        
        echo "Step 1: Extracting FASTQ files from ${SRR_ID}.sra..."
        fasterq-dump --split-files --progress "${SRR_ID}.sra"

        echo "Step 2: Compressing FASTQ files (this may take a few minutes)..."
        gzip "${SRR_ID}_1.fastq"
        gzip "${SRR_ID}_2.fastq"

        echo "Step 3: Validating output files..."
        FILE1="${SRR_ID}_1.fastq.gz"
        FILE2="${SRR_ID}_2.fastq.gz"

        # Check if both compressed files exist AND have a file size greater than zero
        if [ -s "$FILE1" ] && [ -s "$FILE2" ]; then
            echo "✅ Validation passed: $FILE1 and $FILE2 successfully created and are not empty."
            echo "🗑️ Deleting original ${SRR_ID}.sra to save space..."
            rm "${SRR_ID}.sra"
            echo "✅ Finished processing $SRR_ID"
        else
            # If we reach here, something went wrong but didn't trigger a hard crash
            echo "❌ ERROR: Validation failed for $SRR_ID. Output files missing or empty!"
            echo "⚠️ Keeping the original .sra file."
            # Exit the script to prevent further issues
            exit 1
        fi

    else
        echo "⚠️ Warning: ${SRR_ID}.sra not found in ${SRR_ID} directory. Skipping."
    fi

    # Go back up to the main directory to process the next folder
    cd ..
done

echo "====================================================="
echo "🎉 All SRA files have been successfully converted and verified!"
echo "====================================================="
