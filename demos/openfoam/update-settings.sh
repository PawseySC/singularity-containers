#!/bin/bash

# Gloabl settings
aeg_numberOfSubdomains="2" # must equal $SLURM_NTASKS in the main run
aeg_simpleCoeffs="2 1 1" # product must equal numberOfSubdomains above
# Applying global settings
sed -i -e "s,^ *numberOfSubdomains.*,numberOfSubdomains  $aeg_numberOfSubdomains;," \
       -e "/^ *simpleCoeffs/ { n; n; s,^ *n *(.*,    n           ($aeg_simpleCoeffs);,}" system/decomposeParDict

# Settings for main run
aeg_startFrom="startTime" # "latestTime" # for a restart
aeg_startTime="0"
aeg_endTime="20"
aeg_writeInterval="5"
aeg_purgeWrite="0"
# Applying settings for main run
sed -i -e "s,^ *startFrom.*,startFrom    $aeg_startFrom;," \
       -e "s,^ *startTime.*,startTime    $aeg_startTime;," \
       -e "s,^ *endTime.*,endTime    $aeg_endTime;," \
       -e "s,^ *writeInterval.*,writeInterval    $aeg_writeInterval;," \
       -e "s,^ *purgeWrite.*,purgeWrite    $aeg_purgeWrite;," system/controlDict

# Editing SLURM scripts as well
sed -i -e "s/ntasks=.*/ntasks=$aeg_numberOfSubdomains/" -e "s/ntasks\-per\-node=.*/ntasks\-per\-node=$aeg_numberOfSubdomains/" mpi_pawsey.sh
sed -i "s/NTASKS=.*/NTASKS=\"$aeg_numberOfSubdomains\"/" mpi_mpirun.sh
