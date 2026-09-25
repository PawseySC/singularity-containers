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
    "docker://quay.io/pawsey/pytorch:2.7.1-rocm6.3.3|pytorch--2.7.1-rocm6.3.3.sif|pytorch"
)


# Display help for the launcher.
usage() {
    cat <<EOF_USAGE
Usage: $(basename "$0") [OPTIONS]

Options:
  -p, --partition PARTITION       Override the partition in the Slurm script
      --partition=PARTITION       Same as above
      --reservation RESERVATION   Use the specified Slurm reservation
      --reservation=RESERVATION   Same as above
  -h, --help                      Display this help and exit
EOF_USAGE
}


# Read the supported Slurm-style command-line options.
partition=""
reservation=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--partition)
            [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: $1 requires a value." >&2; exit 1; }
            partition="$2"
            shift 2
            ;;
        --partition=*)
            partition="${1#*=}"
            [[ -n "$partition" ]] || { echo "ERROR: --partition requires a value." >&2; exit 1; }
            shift
            ;;
        --reservation)
            [[ $# -ge 2 && -n "$2" ]] || { echo "ERROR: $1 requires a value." >&2; exit 1; }
            reservation="$2"
            shift 2
            ;;
        --reservation=*)
            reservation="${1#*=}"
            [[ -n "$reservation" ]] || { echo "ERROR: --reservation requires a value." >&2; exit 1; }
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unsupported option: $1" >&2
            echo >&2
            usage >&2
            exit 1
            ;;
    esac
done


# Build the optional arguments that will override or extend the Slurm script.
SBATCH_OPTIONS=()

if [[ -n "$partition" ]]; then
    SBATCH_OPTIONS+=(--partition="$partition")
fi

if [[ -n "$reservation" ]]; then
    SBATCH_OPTIONS+=(--reservation="$reservation")
fi


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

    submission_options=(
        --parsable
        --job-name="$job_name"
        "${SBATCH_OPTIONS[@]}"
    )

    if [[ -n "$previous_job_id" ]]; then
        submission_options+=(--dependency="afterany:${previous_job_id}")
    fi

    job_id="$(
        sbatch \
            "${submission_options[@]}" \
            "$JOB_SCRIPT" \
            "$image_reference" \
            "$output_filename"
    )"

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
