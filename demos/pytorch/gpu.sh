#!/bin/bash --login

#SBATCH --job-name=pytorch-gpu
#SBATCH --partition=gpu
#SBATCH --account=<your-project>-gpu   # IMPORTANT: use your own project, with the -gpu suffix
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=pytorch_gpu.out

image="${MYSOFTWARE}/singularity/images/pytorch--2.7.1-rocm6.3.3.sif"
module load singularity/4.1.0-mpi
module load rocm   # provides rocm-smi on the host, for the hardware check below

# Cache the (small) FashionMNIST dataset inside this demo directory, so
# the example is self-contained -- for real work, prefer a persistent
# location such as $MYSOFTWARE.
export DATA_DIR=$(pwd)/mnist_data
mkdir -p $DATA_DIR

# Avoid clashing /tmp usage between concurrent jobs on the same node
export TMPDIR="/tmp/${USER}-${SLURM_JOB_ID}"
mkdir -p $TMPDIR

echo -e "\n\n#------------------------#"
echo "Printing from rocm-smi:"
srun -N 1 -n 1 -c 8 --gres=gpu:1 rocm-smi --showhw

echo -e "\n\n#------------------------#"
echo "Code execution:"
# Note: srun needs its own explicit resource flags for the job step; these
# are independent from the allocation flags above and are not inherited.
# With a single task and a single GPU there's no binding to fine-tune, so
# --gpus-per-task/--gpu-bind aren't needed here.  "-c 8" reserves a full
# CPU chiplet, matched to the GPU chiplet requested via --gres.
srun -l -u -N 1 -n 1 -c 8 --gres=gpu:1 \
    singularity exec --rocm "$image" python3 mnist.py

rm -rf ${TMPDIR}

echo -e "\n\n#------------------------#"
echo "Printing information of finished jobs steps using sacct:"
sacct -j ${SLURM_JOBID} -o jobid%20,Start%20,elapsed%20
