---
title: "Build your own container image"
teaching: 5
exercises: 5
questions:
objectives:
- Learn what is a definition file (def file) and its basic syntax
- Learn how to build a container image and share it with others
- Learn the pros&cons of building with Singularity vs Docker
keypoints:
- Build images using `sudo singularity build`
- Use the remote builder with the flag `-r`, if you need to build images from a machine where you don't have sudo rights
- You can share you Singularity Image File with others, as you would do with any other (big) file
- Upload images to a web registry with `singularity push` (Sylabs account required)
---
