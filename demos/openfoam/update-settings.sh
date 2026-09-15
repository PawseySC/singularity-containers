#!/bin/bash

set -euo pipefail

# New settings
newNumberOfSubdomains="8"  # Must equal Slurm --ntasks
newSimpleCoeffs="1 2 4"    # Product must equal newNumberOfSubdomains
newEndTime="50"
newWriteInterval="10"
newRunTimeModifiable="false"
newMGroup="4"              # I/O rank grouping used by the Slurm script

caseDir="./periodicPlaneChannel"
decomposeParDict="${caseDir}/system/decomposeParDict"
controlDict="${caseDir}/system/controlDict"
slurmScript="./mpi_openfoam_pawsey.slurm.sh"

backup_file() {
    local file=$1
    local index=0
    local backup

    while :; do
        printf -v backup '%s.original.%02d' "$file" "$index"
        if [[ ! -e $backup ]]; then
            cp -p -- "$file" "$backup"
            printf 'Backup created: %s\n' "$backup"
            return
        fi
        ((index += 1))
    done
}

# Validate all inputs before creating backups or changing files.
for file in "$decomposeParDict" "$controlDict" "$slurmScript"; do
    if [[ ! -f $file ]]; then
        printf 'Error: required file not found: %s\n' "$file" >&2
        exit 1
    fi
done

for file in "$decomposeParDict" "$controlDict" "$slurmScript"; do
    backup_file "$file"
done

# Update decomposition settings.
sed -E -i \
    -e "s|^[[:space:]]*numberOfSubdomains[[:space:]]+[^;]+;|numberOfSubdomains  ${newNumberOfSubdomains};|" \
    -e "s|^[[:space:]]*method[[:space:]]+[^;]+;|method          simple;|" \
    -e "s|^[[:space:]]*n[[:space:]]+\([^;]+\);|    n           (${newSimpleCoeffs});|" \
    "$decomposeParDict"

# Update execution settings.
sed -E -i \
    -e "s|^[[:space:]]*endTime[[:space:]]+[^;]+;|endTime         ${newEndTime};|" \
    -e "s|^[[:space:]]*writeInterval[[:space:]]+[^;]+;|writeInterval   ${newWriteInterval};|" \
    -e "s|^[[:space:]]*runTimeModifiable[[:space:]]+[^;]+;|runTimeModifiable ${newRunTimeModifiable};|" \
    "$controlDict"

# Update Slurm resources and OpenFOAM collated-I/O grouping.
sed -E -i \
    -e "s|^(#SBATCH[[:space:]]+--ntasks=).*|\1${newNumberOfSubdomains}|" \
    -e "s|^(#SBATCH[[:space:]]+--ntasks-per-node=).*|\1${newNumberOfSubdomains}|" \
    -e "s|^mGroup=.*|mGroup=${newMGroup}             #Size of the groups for collated fileHandling (32 is the initial recommendation for Setonix)|" \
    "$slurmScript"

printf 'Updated: %s\n' "$decomposeParDict" "$controlDict" "$slurmScript"
