---
title: "Computational Fluid Dynamics with MPI containers"
teaching: 10
exercises: 10
questions:
objectives:
- Learn the steps required to configure and run MPI applications from a container
keypoints:
- You need to build your application in the container with an MPI version which is ABI compatible with MPI libraries in the host
- Appropriate environment variables and bind mounts are required at runtime to make the most out of MPI applications (sys admins are your best friends)
- Singularity transparently interfaces with HPC schedulers such as Slurm
---


