#!/bin/bash

set -euo pipefail

# Define the OpenFOAM image and tutorial case used in this episode.
imageName="openfoam--v2606-gcc13DPInt32Opt-mpich3.4.3-ubuntu24.04.sif"
caseName="periodicPlaneChannel"
tutorialPath="incompressible/pimpleFoam/LES/${caseName}"

# Check that the Pawsey software directory is available.
if [[ -z ${MYSOFTWARE:-} ]]; then
    echo "Error: MYSOFTWARE is not defined." >&2
    exit 1
fi

singularityImage="${MYSOFTWARE}/singularity/images/${imageName}"

# Check that all required commands and files are available.
if ! command -v singularity >/dev/null 2>&1; then
    echo "Error: singularity is not available. Load the required Singularity module first." >&2
    exit 1
fi

if [[ ! -f $singularityImage ]]; then
    echo "Error: OpenFOAM image not found: $singularityImage" >&2
    exit 1
fi

if [[ ! -f ./update-settings.sh ]]; then
    echo "Error: update-settings.sh was not found in the current directory." >&2
    exit 1
fi

if [[ ! -f ./mpi_openfoam_pawsey.slurm.sh ]]; then
    echo "Error: mpi_openfoam_pawsey.slurm.sh was not found in the current directory." >&2
    exit 1
fi

# Stop rather than overwrite a tutorial case prepared during an earlier run.
if [[ -e ./$caseName ]]; then
    echo "Error: ./$caseName already exists." >&2
    echo "Remove or rename the existing directory before running this script again." >&2
    exit 1
fi

# Copy the selected OpenFOAM tutorial from the container image to the host.
echo "Copying the OpenFOAM tutorial case to ./$caseName"
singularity exec "$singularityImage" \
    bash -c 'cp -r "$FOAM_TUTORIALS/$1" "$2"' \
    bash "$tutorialPath" "$PWD"

# Apply the settings required by this episode to the copied case and Slurm script.
echo "Applying the settings required for this episode"
bash ./update-settings.sh

# Report successful completion and the location of the prepared case.
echo "Tutorial preparation complete: $PWD/$caseName"
