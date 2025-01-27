# Bulk_rna_Seq

**Introduction**

RNA sequencing (RNA-seq) is a powerful technique used to analyze the transcriptome of an organism, providing insights into gene expression patterns under various conditions. In this guide, we outline a comprehensive pipeline for processing RNA-seq data, specifically tailored raw data acquisition that we generated out lab to complete downstream analysis. This pipeline encompasses:

Data Acquisition: Retrieving raw sequencing data.
Quality Control: Assessing and ensuring the quality of the raw data.
Read Alignment: Mapping sequencing reads to a reference genome.
Quantification: Measuring gene expression levels.
Differential Expression Analysis: Identifying genes with significant expression changes between conditions.
Visualization: Creating informative plots to represent the data.


The RNA-Seq data preprocessing workflow was executed on the UB-HPC (CCR) cluster using the nf-core/rnaseq pipeline. The goal was to process raw FASTQ files into normalized expression matrices suitable for downstream differential expression gene analysis and pathway enrichment analyses.

###Job Submission on HPC
Below is the SLURM script used to submit the job on the HPC cluster:

### Job Submission on HPC

Below is the SLURM script used to submit the job on the HPC cluster:

```bash
#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=256G
#SBATCH --job-name=
#SBATCH --output=
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

nextflow run /../nf-core-rnaseq-3.12.0/workflow/ \
  --input /../samplesheet.csv \
  --outdir /.. \
  --genome GRCh37 \
  -profile singularity \
  --aligner star_rsem \
  --max_time 72.h \
  --max_memory 224.GB \
  --max_cpus 32
