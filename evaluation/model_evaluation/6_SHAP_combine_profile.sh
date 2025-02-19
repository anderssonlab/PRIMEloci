# Function to combine BED files based on prefix before "_chr"
# Takes two arguments: input directory and output directory
combine_bed_files() {
  local input_dir="$1"
  local output_dir="$2"

  # Ensure the output directory exists
  mkdir -p "$output_dir"

  # Step 1: Find all unique prefixes before "_chr"
  prefixes=$(ls "$input_dir"/*.csv | sed -E 's/(.*)(_chr[^_]+).*/\1/' | sort -u)

  # Step 2: For each unique prefix, concatenate all matching files
  for prefix in $prefixes; do
    # Find all files that match the prefix (with any chr)
    matching_files=$(ls "$input_dir"/$(basename "$prefix")_chr*.csv)

    # Combine those files into one file in the output directory
    combined_file="$output_dir/$(basename "$prefix")_combined.csv"

    # Handle headers: Include header from the first file, skip for others
    head -n 1 $(echo $matching_files | cut -d ' ' -f1) > "$combined_file"  # Extract header from the first file
    for file in $matching_files; do
      tail -n +2 "$file" >> "$combined_file"  # Skip header (start from line 2)
    done

    echo "Combined files into $combined_file"
  done

}

combine_bed_files /Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres_v2/K562_endres_profiles/profiles_subtnorm/K562_N /Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres_v2/K562_endres_profiles/profiles_subtnorm
combine_bed_files /Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres_v2/K562_endres_profiles/profiles_subtnorm/K562_C /Users/natsudanav/Documents/data_PRIMEloci_dev/K562_endres_v2/K562_endres_profiles/profiles_subtnorm