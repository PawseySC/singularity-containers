---
title: "Deep learning with GPU containers: a PyTorch example"
teaching: 5
exercises: 10
questions:
objectives:
- Run a PyTorch GPU training script inside a Singularity container
- Recognise container PATH/entrypoint quirks, and know how to work around them
keypoints:
- Deep learning frameworks such as PyTorch use the same `--rocm`/`--nv` container pattern as any other GPU application
- A container's `PATH` and command names may not match what you expect (*e.g.* `python` vs `python3`) — it's worth checking what's actually there
- Keep datasets and outputs outside the read-only container image, *e.g.* in a bind-mounted or current-working-directory path
---


> ## Note
>
> Setonix's GPU compute nodes are equipped with **AMD Instinct MI250X** GPUs, so we use the `--rocm` flag and a ROCm-enabled container.  If you're on Nvidia GPUs instead, swap `--rocm` for `--nv` and use a CUDA-enabled PyTorch container.
{: .callout}


### Why PyTorch?

Many scientific applications now ship GPU-enabled containers — from molecular dynamics packages like Gromacs to deep learning frameworks like PyTorch.  PyTorch in particular is the basis of a large and growing range of scientific computations, from image analysis to protein structure prediction tools such as AlphaFold.  Here, we'll use it to test the container workflow: a GPU-enabled container, a bind-mounted data directory, and a single-GPU Slurm allocation.  We'll run a a fully-connected network classifying *FashionMNIST* images to help us appreciate GPU-enabled container mechanics.


### Request a new interactive allocation

For this episode we need a compute node with a GPU, rather than the login node.  We do this by requesting an interactive allocation from the scheduler.

(If you already have an interactive allocation open from earlier in the workshop, exit it first by typing `exit`.)

Load the Singularity module first, so it's available once we're on the compute node:

```
$ module load singularity/4.1.0-mpi
```
{: .bash}

Now start a new interactive session on a compute node with:

```
$ salloc -p gpu -A courses01-gpu --gres=gpu:1 -N 1  --reservation=ContainersTraining-gpu -t 00:30:00
```
{: .bash}

```
salloc: Granted job allocation 3453895
salloc: Waiting for resource configuration
salloc: Nodes nid002928 are ready for job
```
{: .output}


### Prepare for the hands-on exercise

*(Skip this if you already have `$TUTO` set from earlier in the workshop.)*

If you haven't done so already, move to a suitable working directory and download the tutorial repository.  On Pawsey systems, use your scratch directory:

```
$ cd "$MYSCRATCH"
$ git clone https://github.com/PawseySC/singularity-containers
$ export TUTO="$PWD/singularity-containers"
$ cd "$TUTO"
```
{: .bash}

Now move to the working directory for this episode:

```
$ cd demos/pytorch
$ pwd
```
{: .bash}

The working directory should be something like:

```
/path/to/scratch/singularity-containers/demos/pytorch
```
{: .output}


### Choose a PyTorch image from Pawsey's registry

In your web browser, go to <https://quay.io/pawsey>.

Pawsey provides container images for several research applications, including PyTorch.  Find and select the `pawsey/pytorch` repository, then select the *Tags* icon on the left side of the page.

Find the PyTorch image with the tag `2.7.1-rocm6.3.3`.  Select the *Fetch Tag* icon on the right, then select *Docker Pull (by tag)*.  Copy only the image tag, not the complete Docker command, and close the window.

> ## Container images and GPU vendors
>
> A container built against CUDA (Nvidia) won't run on an AMD GPU, and vice versa — always match the container's platform to the GPU you're running on.  Pawsey's `pawsey/pytorch` images are built against ROCm, matching Setonix's AMD GPUs.
{: .callout}


### Pull the PyTorch image into your personal library

In the terminal running within the interactive allocation on Setonix, define your personal library directory and create it if it does not already exist:

```
$ export MY_LOCAL_LIBRARY="${MYSOFTWARE}/singularity/images"
$ mkdir -p "$MY_LOCAL_LIBRARY"
```
{: .bash}

Pull the image for this example into your personal library directory:

```
$ singularity pull \
  "${MY_LOCAL_LIBRARY}/pytorch--2.7.1-rocm6.3.3.sif" \
  docker://quay.io/pawsey/pytorch:2.7.1-rocm6.3.3
```
{: .bash}

Pulling the OCI image and converting it into a SIF image may take a few minutes.  For the rest of this episode, we'll refer to it as `$image`:

```
$ export image="${MY_LOCAL_LIBRARY}/pytorch--2.7.1-rocm6.3.3.sif"
```
{: .bash}


### A quick sanity check

Before running a full training job, it's worth checking interactively that the container can actually see the GPU:

```
$ singularity exec --rocm "$image" python3 -c '
import torch
print("PyTorch:", torch.__version__)
print("HIP:", torch.version.hip)
print("CUDA available:", torch.cuda.is_available())
print("GPU count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("GPU:", torch.cuda.get_device_name(0))
'
```
{: .bash}

```
PyTorch: 2.7.1a0+gite2d141d
HIP: 6.3.42134-a9a80e791
CUDA available: True
GPU count: 1
GPU: AMD Instinct MI250X
```
{: .output}

(PyTorch's ROCm build keeps the `cuda` naming throughout its API — `torch.cuda.*` and `device="cuda"` are correct here, they just address the ROCm/HIP backend under the hood.)

> ## Have you spotted the `python` vs `python3` gotcha yet?
>
> If you try `singularity exec --rocm "$image" python -c '...'` (**without** the `3`), you'll get:
> ```
> FATAL:   "python": executable file not found in $PATH
> ```
> {: .output}
> This image simply doesn't alias `python` to `python3` — only `python3` is on `PATH`.  Don't assume a command name inside a container, check what's actually there: `singularity exec "$image" ls /usr/bin | grep python` (or just trying it) will tell you.
{: .callout}

Once you've confirmed the GPU is visible, `exit` the interactive allocation — we'll submit the actual training as a batch job.


### Pre-fetching the dataset

The training script downloads *FashionMNIST* the first time it runs.  Setonix's compute nodes do have outbound internet access, but it's still good practice to fetch a dataset once rather than re-downloading it on every job run.  Let's fetch it now, while we still have our interactive allocation:

```
$ singularity exec --rocm "$image" python3 -c '
from torchvision import datasets
from torchvision.transforms import ToTensor
datasets.FashionMNIST(root="mnist_data", train=True, download=True, transform=ToTensor())
datasets.FashionMNIST(root="mnist_data", train=False, download=True, transform=ToTensor())
'
```
{: .bash}

This creates a `mnist_data` directory in the current folder (`demos/pytorch`) with the dataset already downloaded, so the batch job below won't need to reach out to the internet at all.


### Run the training as a batch job

The current directory has a training script, `mnist.py` — a small fully-connected network, trained for 10 epochs.  It reads where to find the dataset from the `DATA_DIR` environment variable, so the Slurm script controls that:

```
$ cat mnist.py
```
{: .bash}

And here's the Slurm batch script, `gpu.sh`:

```
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
module load singularity/4.1.0-nompi
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
srun -l -u -N 1 -n 1 -c 8 --gres=gpu:1 \
    singularity exec --rocm "$image" python3 mnist.py

rm -rf ${TMPDIR}

echo -e "\n\n#------------------------#"
echo "Printing information of finished jobs steps using sacct:"
sacct -j ${SLURM_JOBID} -o jobid%20,Start%20,elapsed%20
```
{: .bash}

This follows a familiar shape: `srun` combined with `singularity exec --rocm`, plus the account `-gpu` suffix required for GPU jobs on Setonix.  We can submit it with:

```
$ sbatch --account=courses01-gpu --reservation=ContainersTraining-gpu gpu.sh
```
{: .bash}

Once it completes, check `pytorch_gpu.out` for the training loss printed every 100 batches, decreasing epoch over epoch.


### Where to go from here

The pattern you've just used of a GPU-enabled container, `--rocm`, a bind-mounted/pre-fetched data directory, a single-GPU Slurm allocation  is exactly what larger PyTorch-based models (including AlphaFold3-style pipelines) need too.  What changes at that scale is mostly the container image itself (much larger, with model-specific dependencies), the input data (structure/sequence files instead of images), and possibly multi-GPU or multi-node settings — the container mechanics stay the same.
