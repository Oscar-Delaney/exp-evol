#!/bin/bash

# run this script using nohup ./scripts/breseq_script_2.sh &
# you can check what commands are running in the background using 

# Number of threads to use for breseq
NUM_THREADS=8

# Set the output directory name
OUTPUT_DIR="breseq_output_new"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Function to get the reference genome based on the sample name
get_ref_genome() {
    SAMPLE_ID=$1
    # Extract the sample's Name column using awk and determine the reference genome based on its content
    SAMPLE_NAME=$(awk -v ID="$SAMPLE_ID" 'BEGIN{FS="\t"} $1==ID {print $2}' /home/s4528540/liv/sequences_new/samples.tsv)
    if [[ $SAMPLE_NAME == *"_AB3_"* ]]; then
        echo "/home/s4528540/liv/AB3v2.2.gbk"
    elif [[ $SAMPLE_NAME == *"_AB13_"* ]]; then
        echo "/home/s4528540/liv/AB13v2.2.gbk"
    else
        echo ""
    fi
}

# Iterate over the R1 files in the sequences folder
for R1_FILE in /home/s4528540/liv/sequences_new/*_R1_001.fastq.gz; do
    # Extract the sample ID from the R1 filename
    SAMPLE_ID=$(basename "$R1_FILE" | cut -d "_" -f 1)

    # Get the appropriate reference genome for the sample ID
    REF_GENOME=$(get_ref_genome $SAMPLE_ID)
    
    # If no reference genome is found, continue to the next iteration
    if [ -z "$REF_GENOME" ]; then
        echo "No reference genome found for $SAMPLE_ID. Skipping."
        continue
    fi

    # Construct the path to the corresponding R2 file
    R2_FILE="${R1_FILE/_R1_/_R2_}"

    # Construct the breseq command
    BRESEQ_CMD="breseq -r $REF_GENOME $R1_FILE $R2_FILE -p -j $NUM_THREADS -n $SAMPLE_ID -o $OUTPUT_DIR/$SAMPLE_ID"

    # Run the breseq command and capture the exit status
    echo "Running breseq for $SAMPLE_ID..."
    echo "Running command $BRESEQ_CMD "
    $BRESEQ_CMD > $OUTPUT_DIR/$SAMPLE_ID.log
    EXIT_STATUS=$?

    # Check the exit status and print an error message if it is non-zero
    if [ $EXIT_STATUS -ne 0 ]; then
        echo "Error running breseq for $SAMPLE_ID: exit status $EXIT_STATUS"
    fi
done

echo "All breseq commands submitted"
