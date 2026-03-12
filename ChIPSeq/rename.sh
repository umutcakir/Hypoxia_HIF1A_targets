#!/bin/bash

# Read the TXT file line by line and process each line
while IFS=$'\t' read -r run_accession sample_title; do
    # Skip the header line
    if [[ "$run_accession" != "run_accession" ]]; then
        # Remove non-alphanumeric characters from sample title
        clean_sample_title=$(echo "$sample_title" | tr -cd '[:alnum:]_')
        
        # Rename the corresponding file
        old_file="$run_accession.fastq.gz"
        new_file="$run_accession"_"$clean_sample_title".fastq.gz
        
        # Check if the old file exists before renaming
        if [ -f "$old_file" ]; then
            mv "$old_file" "$new_file"
            echo "Renamed $old_file to $new_file"
        else
            echo "File $old_file not found"
        fi
    fi
done < filereport_read_run_PRJNA606242_tsv.txt
