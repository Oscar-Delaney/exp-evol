#!/bin/bash

# Set the path to the breseq output directory
BRESEQ_DIR="output"

# Define the parent directory
PARENT_DIR="analyses"

# Define the sample names files
SAMPLES="data/samples.tsv"

# Set the paths to the main directories for AB3 and AB13
AB3_DIR="${PARENT_DIR}/AB3"
AB13_DIR="${PARENT_DIR}/AB13"

# Create the directory structure if it doesn't exist
mkdir -p "${AB3_DIR}/original"
mkdir -p "${AB3_DIR}/subtracted"
mkdir -p "${AB3_DIR}/results"
mkdir -p "${AB13_DIR}/original"
mkdir -p "${AB13_DIR}/subtracted"
mkdir -p "${AB13_DIR}/results"

# Function to get the target directory based on the sample name
get_target_directory() {
    SAMPLE_ID=$1
    # Extract the sample's Name column using awk and determine the reference genome based on its content
    SAMPLE_NAME=$(awk -v ID="$SAMPLE_ID" 'BEGIN{FS="\t"} $1==ID {print $2}' $SAMPLES)
    if [[ $SAMPLE_NAME == *"AB3"* ]]; then
        echo "$AB3_DIR"
    elif [[ $SAMPLE_NAME == *"AB13"* ]]; then
        echo "$AB13_DIR"
    else
        echo "UNKNOWN"
    fi
}

# Copy the annotated.gd file for each sample to the appropriate directory
for SAMPLE_DIR in "$BRESEQ_DIR"/*/; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    ANNOTATED_GD_FILE="$SAMPLE_DIR/data/annotated.gd"
    TARGET_DIR=$(get_target_directory "$SAMPLE_NAME")
    OUTPUT_FILE="${TARGET_DIR}/original/${SAMPLE_NAME}.gd"
    echo "Copying $ANNOTATED_GD_FILE to $OUTPUT_FILE..."
    cp "$ANNOTATED_GD_FILE" "$OUTPUT_FILE"
done

# # Run the mapping_success.sh script to generate a CSV file of read mapping success percentages
MAP_SUCCESS="${PARENT_DIR}/percent_aligned.csv"
code/mapping_success.sh "$BRESEQ_DIR" "$MAP_SUCCESS"
sed -i 's/\r//g' "$MAP_SUCCESS"
# tr -d '\r' < "$MAP_SUCCESS" > "$MAP_SUCCESS"

# Set the threshold for low read mapping success
THRESHOLD=0.5

# Initialize the list of samples with low read mapping success
LOW_SUCCESS_SAMPLES=""

# Read the percent_aligned.csv file and add the names of samples with low read mapping success to the list
while IFS=, read -r sample percent_aligned_bases; do
    if (( $(echo "$percent_aligned_bases < $THRESHOLD" | bc -l) )); then
        LOW_SUCCESS_SAMPLES="$LOW_SUCCESS_SAMPLES|$sample"
    fi
done < "$MAP_SUCCESS"

# Remove leading |
LOW_SUCCESS_SAMPLES=${LOW_SUCCESS_SAMPLES#|}

echo "$LOW_SUCCESS_SAMPLES"

# Run gdtools for each reference genome
for DIR in $AB3_DIR $AB13_DIR; do
    # Excluding low-success samples for intersect file creation
    ANNOTATED_GD_FILES=$(ls "$DIR/original"/*.gd | grep -v -E $LOW_SUCCESS_SAMPLES)
    OUTPUT_FILE="${DIR}/results/intersect_$(basename $DIR).gd"
    
    # Construct the gdtools intersect command
    GDTOOLS_CMD="gdtools intersect -o $OUTPUT_FILE $ANNOTATED_GD_FILES"

    echo "Running gdtools intersect command for $(basename $DIR)..."
    $GDTOOLS_CMD
    EXIT_STATUS=$?

    # Check the exit status and print an error message if it is non-zero
    if [ $EXIT_STATUS -ne 0 ]; then
        echo "Error running gdtools intersect command for $(basename $DIR): exit status $EXIT_STATUS"
    fi

    echo "gdtools intersect command for $(basename $DIR) complete"
done

# Define reference genomes
declare -A REF_GENOMES
REF_GENOMES["$AB3_DIR"]="genomes/AB3v2.2.gbk"
REF_GENOMES["$AB13_DIR"]="genomes/AB13v2.2.gbk"

for DIR in "$AB3_DIR" "$AB13_DIR"; do
    # Excluding low-success samples for intersect file creation
    ANNOTATED_GD_FILES=$(ls "$DIR/original"/*.gd | grep -v -E $LOW_SUCCESS_SAMPLES)
    INTERSECT_FILE="${DIR}/results/intersect_$(basename $DIR).gd"

    # Subtract the intersection from each gd file
    for SAMPLE in $ANNOTATED_GD_FILES; do
        SAMPLE_NAME=$(basename "$SAMPLE")
        # Exclude samples with low mapped percentages
        OUTPUT_FILE="${DIR}/subtracted/${SAMPLE_NAME}"
        GDTOOLS_CMD2="gdtools subtract -o $OUTPUT_FILE $SAMPLE $INTERSECT_FILE"
        $GDTOOLS_CMD2
    done
    
    OUTPUT_FILE2="${DIR}/results/comparison.html"
    GDTOOLS_CMD3="gdtools compare -o $OUTPUT_FILE2 -f html -r ${REF_GENOMES[$DIR]} ${DIR}/original/*"
    $GDTOOLS_CMD3

    OUTPUT_FILE3="${DIR}/results/simplified_comparison.html"
    GDTOOLS_CMD4="gdtools compare -o $OUTPUT_FILE3 -f html -r ${REF_GENOMES[$DIR]} ${DIR}/subtracted/*"
    $GDTOOLS_CMD4

    OUTPUT_FILE4="${DIR}/results/simplified_comparison.csv"
    GDTOOLS_CMD5="gdtools compare -o $OUTPUT_FILE4 -f csv -r ${REF_GENOMES[$DIR]} ${DIR}/subtracted/*"
    $GDTOOLS_CMD5
done
