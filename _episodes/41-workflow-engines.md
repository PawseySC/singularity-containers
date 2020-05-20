---
title: "Reproducible scientific workflows"
teaching: 5
exercises: 10
questions:
objectives:
- Get an idea of the interplay between containers and workflow engines
keypoints:
- Some workflow engines offer transparent APIs for running containerised applications
- If you need to run data analysis pipelines, the combination of containers and workflow engines can really make your life easier!
---


### Scientific workflow engines

Let's see how Singularity containers can be used in conjunction with a popular workflow engine.

Scientific workflow engines are particularly useful for data-intensive domains including (and not restricted to) bioinformatics and radioastronomy, where data analysis and processing is made up of a number of tasks to be repeatedly executed across large datasets.  Some of the most popular ones, including [Nextflow](https://www.nextflow.io) and [Snakemake](https://snakemake.readthedocs.io), provide interfaces to container engines.  The combination of container and workflow engines can be very effective in enforcing reproducible, portable, scalable science.

Now, let's try and use Singularity and Nextflow to run a demo RNA sequencing pipeline based on [RNAseq-NF](https://github.com/nextflow-io/rnaseq-nf).


### Install Nextflow

First, if it's not already on your system, you'll need to install Nextflow.  You'll need to install a Java runtime and download the Nextflow executable.  It will take a few minutes to download all of the required dependencies, but the process is fairly automated.  This is a template install [script]({{ page.root }}/files/install-nextflow.sh) for a Linux box.

If you're running on the Pawsey Nimbus cloud, just run the above script via: `bash $TUTO/files/install-nextflow.sh`.

If you're running at Pawsey *e.g.* on Zeus, all you need is to `module load nextflow`.


### Run a workflow using Singularity and Nextflow

Let's `cd` into the appropriate directory:

```
$ cd $TUTO/demos/nextflow
```
{: .bash}

For convenience, the content of the pipeline [RNAseq-NF](https://github.com/nextflow-io/rnaseq-nf) is already made available in this directory.  There are two critical files in here, namely `main.nf`, that contains the translation of the scientific pipeline in the Nextflow language, and `nextflow.config`, that contains several profiles for running with different software/hardware setups. Here we are going to use the profile called `singularity`.

It's time to launch the pipeline with Nextflow:

```
$ nextflow run main.nf -profile singularity
```
{: .bash}

We'll get some information on the pipeline, along with the notice that the appropriate container is being downloaded:

```
N E X T F L O W  ~  version 19.10.0
Pulling marcodelapierre/rnaseq-nf ...
 downloaded from https://github.com/marcodelapierre/rnaseq-nf.git
Launching `marcodelapierre/rnaseq-nf` [hopeful_almeida] - revision: 91dd162c00 [master]
 R N A S E Q - N F   P I P E L I N E
 ===================================
 transcriptome: /data/work/.nextflow/assets/marcodelapierre/rnaseq-nf/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa
 reads        : /data/work/.nextflow/assets/marcodelapierre/rnaseq-nf/data/ggal/ggal_gut_{1,2}.fq
 outdir       : results

WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: /data/work/singularity-test/nxf/work/singularity
Pulling Singularity image docker://nextflow/rnaseq-nf:latest [cache /data/work/singularity-test/nxf/work/singularity/nextflow-rnaseq-nf-latest.img]
```
{: .output}

It will take a bunch of minutes to download the container image, then the pipeline will run:

```
[9e/a8a999] Submitted process > fastqc (FASTQC on ggal_gut)
[6a/4ec5ee] Submitted process > index (ggal_1_48850000_49020000)
[91/109c65] Submitted process > quant (ggal_gut)
[ab/081287] Submitted process > multiqc

Done! Open the following report in your browser --> results/multiqc_report.html

```
{: .output}

The final output of this pipeline is an HTML report of a quality control task, which you might eventually want to download and open up in your browser.  

However, the key question here is: how could the sole flag `-profile singularity` trigger the containerised execution?  This is the relevant snippet from the `nextflow.config` file:

```
  singularity {
    process.container = 'nextflow/rnaseq-nf:latest'
    singularity.enabled = true
    singularity.autoMounts = true
  }
```
{: .source}

The image name is specified using the `process.container` keyword.  Also, `singularity.autoMounts` is required to have the directory paths with the input files automatically bind mounted in the container.  Finally, `singularity.enabled` triggers the use of Singularity.

Based on this configuration file, Nextflow is able to handle all of the relevant Singularity commands by itself, *i.e.* `pull` and `exec` with the appropriate flags, such as `-B` for bind mounting host directories.  In this case, as a user you don't need to know in detail the Singularity syntax, but just the name of the container!

More information on configuring Nextflow to run Singularity containers can be found at [Singularity containers](https://www.nextflow.io/docs/latest/singularity.html).
