#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#@@. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#AA. ../../RunFunctions        # Tutorial run functions
. ./RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

#@@
# Additional parameters for the openfoam commands
additionalPar="-fileHandler uncollated"

# Finally: decompose, solve the transient case in parallel & reconstruct for postprocessing
(
    cd wingMotion2D_pimpleFoam || exit

    # Decompose the domain for parallel execution:
    singularity exec $theImage decomposePar -force $additionalPar | tee log.decomposePar

    # Execute the solver in parallel:
    #@@mpirun -n 4 singularity exec $theImage pimpleFoam -parallel $additionalPar | tee log.pimpleFoam
    srun -n 4 singularity exec $theImage pimpleFoam -parallel $additionalPar | tee log.pimpleFoam

    # Reconstruct the decomposed results into full-domain results
    singularity exec $theImage reconstructPar $additionalPar | tee log.reconstructPar

    # Create dummy file for paraview loader
    touch wingMotion2D_pimpleFoam.foam
)

#------------------------------------------------------------------------------
