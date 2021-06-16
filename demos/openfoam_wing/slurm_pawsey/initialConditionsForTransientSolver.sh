#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#@@. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#AA. ../../RunFunctions        # Tutorial run functions
. ./RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

#@@
# Additional parameters for the openfoam commands
additionalPar="-fileHandler uncollated"

# Make 3D mesh in slab of cells.
(
    cd wingMotion_snappyHexMesh || exit

    singularity exec $theImage blockMesh $additionalPar | tee log.blockMesh
    singularity exec $theImage snappyHexMesh -overwrite $additionalPar | tee log.snappyHexMesh
)

# Make a 2D mesh by extruding a patch and solve to steady state.
(
    cd wingMotion2D_simpleFoam || exit

    singularity exec $theImage extrudeMesh $additionalPar | tee log.extrudeMesh
    singularity exec $theImage createPatch -overwrite $additionalPar | tee log.createPatch
    restore0Dir
    singularity exec $theImage simpleFoam $additionalPar | tee log.simpleFoam
)

# Copy mesh from the steady state case, map the results to a mesh motion case,
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
