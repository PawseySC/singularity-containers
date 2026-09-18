# MPI Mandelbrot container example

The image contains the MPI application and its runtime dependencies. Launch policy remains outside the image:

- `run-mandelbrot-docker.sh` performs a local, single-host Docker test using container-side `mpiexec`.
- `mpi_mandelbrot_pawsey.slurm.sh` uses host-side `srun`, which starts one `singularity exec` per Slurm task on Setonix.

Build with:

    docker build --platform linux/amd64 --file mandelbrot_mpi.dockerfile --tag mandelbrot-mpi:2026.09 .
