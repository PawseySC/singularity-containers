#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
#@@. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions      # Tutorial clean functions
#@@. ../../CleanFunctions      # Tutorial clean functions
. ./CleanFunctions      # Tutorial clean functions
#------------------------------------------------------------------------------

# Clean the mesh generation
( cd wingMotion_snappyHexMesh && cleanCase )
# Clean the potential solution case (initial conditions)
( cd wingMotion2D_simpleFoam && cleanCase0 )
# Clean the transient solution case
( cd wingMotion2D_pimpleFoam && cleanCase0 )

#------------------------------------------------------------------------------
