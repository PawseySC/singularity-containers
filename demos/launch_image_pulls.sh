#!/bin/bash

# Stop if a command fails, if an unset variable is used, or if a pipeline fails.
set -euo pipefail


# List the images to pull.
#
# Each entry contains:
#   IMAGE_REFERENCE|OUTPUT_FILENAME|SHORT_NAME
#
# Add or remove entries only in this block.
IMAGES=(
    "docker://quay.io/pawsey/openfoam:v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04|openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif|openfoam"
    "docker://docker.io/trinityrnaseq/trinityrnaseq:2.8.6|trinityrnaseq--2.8.6.sif|trinity"
    "docker://quay.io/pawsey/pytorch:2.7.1-rocm6.3.3.sif|pytorch--2.7.1-rocm6.3.3.sif|pytorch"
)


# Set the personal image library used throughout the training.
export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"


# Locate the Slurm job script in the same demos directory as this launcher.
DEMOS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="${DEMOS_DIR}/pull_image.slurm.sh"


# Check that the Slurm job script is available.
if [[ ! -f "$JOB_SCRIPT" ]]; then
    echo "ERROR: Slurm job script not found: $JOB_SCRIPT" >&2
    exit 1
fi


# Create the personal image library if it does not already exist.
mkdir -p "$MY_LOCAL_LIBRARY"


# Initialise the dependency chain.
previous_job_id=""


# Submit one pull job for each image.
# Each job after the first waits for the preceding job to finish.
for image in "${IMAGES[@]}"; do
    IFS='|' read -r image_reference output_filename short_name <<< "$image"

    echo "Submitting pull job for: $image_reference"
    echo "Image will be saved as:  ${MY_LOCAL_LIBRARY}/${output_filename}"
    job_name="pulling:${short_name}"
    echo "Slurm job name:          $job_name"

    if [[ -z "$previous_job_id" ]]; then
        job_id="$(
            sbatch --parsable \
                --job-name="$job_name" \
                --partition="gpu" \
                --account="courses01-gpu" \
                -N 1 \
                --gres=gpu:1 \
                --reservation="ContainersTraining-gpu" \
                "$JOB_SCRIPT" \
                "$image_reference" \
                "$output_filename"
        )"
    else
        job_id="$(
            sbatch --parsable \
                --job-name="$job_name" \
                --partition="gpu" \
                --account="courses01-gpu" \
                -N 1 \
                --gres=gpu:1 \
                --reservation="ContainersTraining-gpu" \
                --dependency="afterany:${previous_job_id}" \
                "$JOB_SCRIPT" \
                "$image_reference" \
                "$output_filename"
        )"
    fi

    # Remove a cluster name if sbatch returns JOB_ID;CLUSTER_NAME.
    job_id="${job_id%%;*}"

    echo "Submitted Slurm job:     $job_id"

    previous_job_id="$job_id"

    echo
done


# Report that all pull jobs have been submitted.
echo "All image-pulling jobs have been submitted."
echo "Each job will start after the preceding image-pulling job finishes."
echo "Monitor them with: squeue --me"
