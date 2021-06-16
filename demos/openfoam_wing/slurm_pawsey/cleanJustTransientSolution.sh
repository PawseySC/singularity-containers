#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#@@. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions      # Tutorial clean functions
#AA. ../../CleanFunctions      # Tutorial clean functions
#AA. ../../RunFunctions        # Tutorial run functions
. ./CleanFunctions      # Tutorial clean functions
. ./RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

#@@
# Additional parameters for the openfoam commands
additionalPar="-fileHandler uncollated"

( cd wingMotion2D_pimpleFoam && cleanCase0 )

# Copy mesh from the steady state case, map the results to a mesh motion case.
(
    cd wingMotion2D_pimpleFoam || exit

    rm -rf constant/polyMesh
    \cp -r ../wingMotion2D_simpleFoam/constant/polyMesh constant
    restore0Dir
    singularity exec $theImage mapFields ../wingMotion2D_simpleFoam -sourceTime \
                                    latestTime -consistent $additionalPar | tee log.mapFields
    \mv 0/pointDisplacement.unmapped 0/pointDisplacement
)

#------------------------------------------------------------------------------
