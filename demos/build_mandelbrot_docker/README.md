# MPI Mandelbrot container example

The workload and viewport are controlled through environment variables.

Local Docker defaults:

- `MPI_PROCESSES=4`
- `WIDTH=1200`
- `HEIGHT=800`
- `ITERATIONS=500`
- `CENTRE_REAL=-0.5`
- `CENTRE_IMAGINARY=0.0`
- `SCALE=3.0`

Run locally with the defaults:

    ./runMandelbrotDocker.sh

The Setonix script requests 16 tasks and defaults to a larger, zoomed workload:

- `WIDTH=6000`
- `HEIGHT=4000`
- `ITERATIONS=2000`
- `CENTRE_REAL=-0.743643887037151`
- `CENTRE_IMAGINARY=0.131825904205330`
- `SCALE=0.002`

Submit with its defaults:

    sbatch runMandelbrotSingularityPawsey.slurm.sh

Any value can be overridden through `sbatch --export`.
