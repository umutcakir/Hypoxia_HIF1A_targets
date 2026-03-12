#!/bin/bash
#SBATCH -p soeding
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 14-0:00:00
#SBATCH --mem=50G

awk 'BEGIN { FS=","; OFS="," } { gsub(/.*\(|\).*/, "", $3); gsub(/^ *| *$/, "", $1); gsub(/^ *| *$/, "", $2); gsub(/^ *| *$/, "", $3); $3 = "encTfChipPk" $3 ".txt.gz"; print }' encode_list_of_files.csv | sed 's/\//_/g' | tr ' ' '-' > encode_list.csv

# Output directory
output_dir="EncodeChIP"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Read each line of the CSV file
cat encode_list.csv | while IFS=',' read -r cell_type factor track_name; do
   # Construct the URL to download the file
   url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/${track_name}"
   echo "$url"
   
   # Construct the new file name
   new_file_name="${cell_type}_${factor}_${track_name}"

   # Download the file into the output directory
   wget -P "$output_dir" "$url" -O "$output_dir/$new_file_name"
done

mkdir EncodeData

cp EncodeChIP/*txt.gz EncodeData/

cd EncodeData/

for file in *; do
   # Check if the file is compressed
   if [ -f "$file" ]; then
       # Get the file extension
       extension="${file##*.}"
       
       # Uncompress based on the file extension
       case "$extension" in
           zip)
               unzip "$file" ;;
           gz)
               gunzip "$file" ;;
           tar)
               tar -xf "$file" ;;
           # Add more cases for other compression formats if needed
       esac
   fi
done

cd ..

# Directory containing the files
directory="EncodeData"

# Loop through each file in the directory
for filename in $directory/*.txt; do
   # Extract the second argument from the filename
   second_argument=$(basename "$filename" | cut -d '_' -f 2)
   
   # Create a directory to store merged files if it doesn't exist
   merged_dir="$directory/merged"
   mkdir -p "$merged_dir"
   
   # Append content of the file to the corresponding merged file
   cat "$filename" >> "$merged_dir/$second_argument.txt"
done

mkdir "ChIP_Peaks"

input_directory="EncodeData/merged/"

# Define the output directory where new files will be saved
output_directory="ChIP_Peaks/"

# Loop through each input file in the directory
for input_file in "$input_directory"/*.txt; do
   # Get the name of the input file without the extension
   input_file_name=$(basename "$input_file" .txt)
   
   # Initialize line counter
   line_number=1
   
   # Define the output file name
   output_file="$output_directory/${input_file_name}_summits.bed"
   
   # Process each line in the input file
   while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11; do
       # Calculate the values for the new columns
       second_column=$((col3 + col11)) #We are using col11 to get the peak
       third_column=$((second_column + 1))
       
       # Generate the name for the fourth column with line number
       fourth_column_name="${input_file_name}_${line_number}"
       
       # Write the new values to the output file
       echo -e "$col2\t$second_column\t$third_column\t$fourth_column_name\t$col8" >> "$output_file"
       
       # Increment line counter
       ((line_number++))
   done < "$input_file"
done


zcat Repeats_hg38.tsv.gz \
| awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 {next}  # skip header
    {
      chrom=$6; start0=$7; end0=$8; strand=$10;
      name=$11; cls=$12; fam=$13;

      # keep only main chromosomes
      if (chrom !~ /^chr([0-9]{1,2}|X|Y|M)$/) next;
 
      # exclude entries with "?" in any field
      # if (name ~ /\?/ || cls ~ /\?/ || fam ~ /\?/) next; #We are not removing them!

      # exclude unwanted classes
      if (cls ~ /^(Low_complexity|Simple_repeat|Satellite|rRNA|tRNA|scRNA|snRNA|RNA|RC|srpRNA)$/) next;

      # exclude too short
      # if ((end0-start0) < 100) next; #For enrichment, we do not remove based on length of TEs

      # output as BED
      print chrom,start0,end0,fam,cls,name,0,strand;
    }' > Repeats_hg38_reduced.bed




source activate r_env

mkdir "Enrichment_Scores"



# Get a list of all BED files in ChIP_Peaks directory
files=(ChIP_Peaks/*.bed)

# Determine the total number of files
total_files=${#files[@]}

# Set the number of files to process in each SLURM job
files_per_job=12

# Calculate the number of SLURM jobs needed
num_jobs=$(( (total_files + files_per_job - 1) / files_per_job ))

# Loop through each SLURM job
for ((job_id = 0; job_id < num_jobs; job_id++)); do
    # Calculate the start and end index for the current job
    start_index=$(( job_id * files_per_job ))
    end_index=$(( (job_id + 1) * files_per_job ))

    # Slice the files array to get the subset for the current job
    batch_files=("${files[@]:start_index:files_per_job}")

    # Create a temporary script file for SLURM job
    script_file="slurm_job_${job_id}.sh"
    printf '%s\n' '#!/bin/bash' > "$script_file"
    echo "#SBATCH -p soeding" >> "$script_file"
    echo "#SBATCH -N 1" >> "$script_file"
    echo "#SBATCH -n 24" >> "$script_file"
    echo "#SBATCH -t 14-0:00:00" >> "$script_file"
    echo "#SBATCH --mem=100G" >> "$script_file"
    echo "#SBATCH --constraint=inet" >> "$script_file"
    echo "" >> "$script_file"
    echo "source activate r_env" >> "$script_file"
    echo "" >> "$script_file"
    # Run R scripts in parallel
    echo "parallel -j 12 <<EOF" >> "$script_file"
    for file in "${batch_files[@]}"; do
        # Extract the filename without extension
        filename=$(basename "$file" .bed)
        # Run your R script with the current file (without extension) as an argument
        echo "Rscript enrichment.R 'Repeats_hg38_reduced.bed' '$file' 'Enrichment_Scores/${filename}.csv'" >> "$script_file"
    done
    echo "EOF" >> "$script_file"

    # Make the script file executable
    chmod +x "$script_file"
    
    # Submit SLURM job
    #sbatch "$script_file"
    #sleep 1
    # Remove temporary script file
    #rm "$script_file"
done
  

# Initialize an empty dataframe to store combined data
combined_dataframe=""

# Flag to check if the first column has been added
first_column_added=false

# Loop through each CSV file in the directory
for file in Enrichment_Scores/*.csv; do
    # Extract the second column from the CSV file
    second_column=$(cut -d',' -f2 "$file")
    
    # Extract the filename without extension
    filename=$(basename "$file" .csv)
    
    # If the first column hasn't been added yet, add it
    if ! $first_column_added; then
        first_column=$(cut -d',' -f1 "$file")
        combined_dataframe="$first_column"
        first_column_added=true
    fi

    # Extract the header of the second column
    second_column_header=$(echo "$second_column" | sed -n '1p')
    
    # Replace the header with the filename
    second_column=$(echo "$second_column" | sed "1s/$second_column_header/$filename/")
    
    # Append the second column to the combined dataframe
    combined_dataframe=$(paste -d, <(echo "$combined_dataframe") <(echo "$second_column"))
done

# Output the combined dataframe to a new CSV file
echo "$combined_dataframe" > combined_enrichment_score.csv




