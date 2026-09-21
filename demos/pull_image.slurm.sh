#!/bin/bash -l

#SBATCH --job-name=pull-image
#SBATCH --partition=work
#SBATCH --reservation=ContainersTraining
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=pull-image-%j.out


# Stop if a command fails, if an unset variable is used, or if a pipeline fails.
set -euo pipefail


# Check the command-line arguments supplied by the launcher.
if [[ $# -ne 2 ]]; then
    echo "Usage: sbatch $0 IMAGE_REFERENCE OUTPUT_FILENAME" >&2
    exit 1
fi


# Read the image reference and output filename supplied by the launcher.
IMAGE_REFERENCE="$1"
OUTPUT_FILENAME="$2"


# Define the personal image library used throughout the training.
MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
OUTPUT_FILE="${MY_LOCAL_LIBRARY}/${OUTPUT_FILENAME}"


# Create the personal image library if it does not already exist.
mkdir -p "$MY_LOCAL_LIBRARY"


# Avoid downloading an image that is already present.
if [[ -f "$OUTPUT_FILE" ]]; then
    echo "Image already exists. Nothing to do:"
    echo "$OUTPUT_FILE"
    exit 0
fi


# Load the Singularity module used for these image downloads.
module load singularity/4.1.0-nompi


# Pull the image into the personal image library.
echo "Pulling image: $IMAGE_REFERENCE"
echo "Saving image:  $OUTPUT_FILE"

singularity pull "$OUTPUT_FILE" "$IMAGE_REFERENCE"


# Confirm that the image was created.
echo
echo "Pull completed successfully:"
ls -lh "$OUTPUT_FILE"
