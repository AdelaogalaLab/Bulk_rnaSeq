# Bulk_rna_Seq

**Introduction**

RNA sequencing (RNA-seq) is a powerful technique used to analyze the transcriptome of an organism, providing insights into gene expression patterns under various conditions. In this guide, we outline a comprehensive pipeline for processing RNA-seq data, specifically tailored from raw data acquisition to complete downstream analysis. This pipeline encompasses:

Data Acquisition: Retrieving raw sequencing data  
Quality Control: Assessing and ensuring the quality of the raw data  
Read Alignment: Mapping sequencing reads to a reference genome  
Quantification: Measuring gene expression levels  
Differential Expression Analysis: Identifying genes with significant expression changes between conditions  
Visualization: Creating informative plots to represent the data

The RNA-Seq data preprocessing workflow was executed on the UB-CCR cluster using the nf-core/rnaseq pipeline. The goal was to process raw FASTQ files into normalized expression matrices suitable for downstream differential expression gene analysis and pathway enrichment analyses.

# Job Submission on HPC
Below is the SLURM script used to submit the job on the UB-CCR cluster:

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

```
Explanation of SLURM Parameters\
--time: Total runtime of the job (72 hours).\
--nodes: Number of nodes requested (1 node with 32 cores).\
--mem: Memory allocation (256 GB).\
--job-name: Job name for easy identification in the job queue.\
--output: Output log file for pipeline progress and errors.\
--partition: HPC partition to use (general-compute).\
--max_time, --max_memory, --max_cpus: Specifies pipeline resource limits.

### Workflow Description

The preprocessing pipeline consists of the following steps:

1. Input Metadata

A samplesheet (samplesheet.csv) provides metadata for all RNA-Seq samples.
Required fields include sample IDs, file paths for paired-end or single-end FASTQ files, and experimental conditions.

2. Quality Control
Tool: FastQC
Assesses base quality scores, GC content, adapter contamination, and other metrics.
MultiQC: Aggregates QC reports for all samples into a single report for easy visualization.

3. Read Alignment
Tool: STAR
Aligns raw reads to the GRCh37 reference genome.
Outputs include:
BAM files (*.Aligned.sortedByCoord.out.bam).
Alignment statistics (*.Log.final.out).

4. Quantification
Tool: RSEM
Quantifies gene- and transcript-level expression.
Outputs:
Gene-level counts (*.rsem.genes.results).
Transcript-level counts (*.rsem.isoforms.results).

5. Normalization
Metric: Transcripts Per Million (TPM).
Ensures comparability across samples by normalizing for sequencing depth and gene length.

Pipeline Outputs
The pipeline generates the following outputs:

Quality Control Reports

FastQC: Individual QC reports for each sample.
MultiQC: Aggregated QC metrics.
Alignment Files

Sorted BAM files for all samples.
Alignment statistics (STAR .Log.final.out files).
Quantification Files
Gene-level and transcript-level expression matrices (*.rsem.genes.results and *.rsem.isoforms.results).

***

## Down stream Analysis
This repository contains the R code and analysis steps for the bulk rna sequencing project. The code performs differential gene expression analysis, pathway enrichment analysis, and visualization of key results. Below is a detailed explanation of the key scripts and their outputs.

## Overview
This project focuses on differential gene expression analysis, GSEA, bubble plot, and pathway enrichment analysis to identify the biological significance of MECOM regulation. Key methodologies include:
- *RNA-Seq normalization and filtering*
- *Visualization of differential expression results using volcano plots, bubble plots, heatmaps, and Venn diagrams*
- *Gene Set Enrichment Analysis (GSEA)*
- *Identification of overlapping upregulated and downregulated genes*

---

## Requirements
To reproduce this analysis, the following tools and R packages are required:

### Tools
- R (version 4.0 or higher)
- RStudio (recommended)

### R Packages
Ensure these packages are installed before running the scripts:
- `edgeR`
- `dplyr`
- `ggplot2`
- `ggrepel`
- `readxl`
- `pheatmap`
- `VennDiagram`
- `ggVennDiagram`
- `msigdbr`
- `clusterProfiler`

Install packages using:

install.packages(c("dplyr", "ggplot2", "ggrepel", "readxl", "pheatmap", "VennDiagram"))
BiocManager::install(c("edgeR", "msigdbr", "clusterProfiler", "ggVennDiagram"))
Data Input

Required Files
Gene expression counts: A TSV file containing normalized counts per gene.
Sample sheet: An Excel file (Sample_Sheet.xlsx) providing metadata for samples.
Pathway shortcuts file: An Excel file (pathways_shortcuts_final.xlsx) linking pathway names with abbreviations.
Ensure all input files are placed in the working directory.

Analysis Workflow

1. Preprocessing and Normalization
Filter out lowly expressed genes using edgeR's cpm and filterByExpr.
Normalize data using the calcNormFactors function.

2. Differential Gene Expression Analysis
Identify upregulated and downregulated genes across selected contrast IDs.
Save significant genes in CSV format for further visualization.

3. Visualization

* PCA plots
A PCA plot helps us see how similar or different the samples are. Samples with similar gene expression will be closer together on the plot. It uses an x-y coordinate system, where the two axes represent the main patterns (principal components) that explain most of the differences in the data.

* Differential Gene Heatmap
A differential gene heatmap shows the normalized values of the top upregulated genes that differ significantly between comparison groups (e.g., normal vs treatment). Genes are filtered by padj and ranked by the highest log2FoldChange. The values are scaled across samples to highlight differences.

* Volcano Plots
A volcano plot is a scatterplot showing log2FoldChange vs -log(padj) for each gene. It highlights genes that are highly upregulated, downregulated, or have significant p-values. We discussed different example for normal vs treatment comparison for different cell lines.Visualize log fold changes (logFC) and significance (-log10(p-value)) of genes for each contrast using ggplot2.

* Venn Diagrams
Identify overlapping genes among upregulated and downregulated groups using ggVennDiagram and VennDiagram.

* Gene Set Enrichment Analysis (GSEA)
GSEA identifies pathways enriched in one condition compared to another. For example, it can determine if genes in a hallmark pathway are collectively upregulated with treatment, helping validate the experiment and identify other affected pathways.
Gene sets from the msigdbr library and custom sets are combined into a single dataset (all_sets).
Genes are ranked by their log fold change (logFC) between selected groups (e.g., normal vs treatment) and filtered by specific contrasts of interest.
The ranked gene lists are analyzed using the GSEA() function from the clusterProfiler package with the combined gene sets as input.
Results are annotated with metadata and categorized into upregulated, downregulated, or unchanged pathways based on the normalized enrichment score (NES) and q-value.
This approach allows efficient pathway analysis across multiple contrasts, providing insights into pathways influenced under different conditions.

* Bubble Plots
Generate bubble plots for enriched pathways using GSEA results. The size of bubbles represents significance, and colors indicate normalized enrichment scores (NES).

* Outputs

Key Files

Normalized Counts: MECOM_Project_rna_mean_normalized-counts-per-million.csv

Volcano Plots: Generated for each contrast ID (e.g., Volcano_C42B_ER_shMECOM_v_control.pdf)

Venn Diagrams:

Up_regulated.pdf
Down_regulated.pdf
Heatmaps: Heatmap_top_genes.pdf

Bubble Plots: Pathway-level visualization for each contrast (e.g., Bubble_Plot_for_LREX_shMECOM_v_control.pdf).

Gene Lists: Overlapping and regulated genes in CSV format.

Example Output Snippets
Overlapping Genes: shMECOM_overlapping_genes_data.csv
All Regulated Genes: All_UpandDown_regulatedGenes_Selected_Contrasts.csv
Visualization Examples

PCA Plot, Volcano Plot, Heatmap, Bubble Plot, Enrichment Plot, and GSEA plot


```
***Whole Exome Sequencing***

Sarek is a robust bioinformatics workflow developed using Nextflow, a specialized framework designed for scalable and reproducible computational tasks. The intuitive structure of Nextflow facilitates easy interpretation, modification, and expansion of workflows. Sarek incorporates efficient error-handling mechanisms to detect and address issues such as software or hardware failures. Furthermore, the workflow supports seamless resumption of incomplete analyses, allowing users to restart from any stage without complications.
Nextflow’s versatility enables Sarek to operate across diverse computational environments, including high-performance computing (HPC) cluster systems (e.g., Slurm) and container orchestration platforms (e.g., Kubernetes). It also integrates with cloud services such as Google Cloud and AWS, making it highly adaptable for varied research needs. As part of the nf-core ecosystem (https://nf-co.re), Sarek benefits from extensive validation and a user-friendly interface, ensuring widespread applicability and reliability.
To enhance reproducibility and ease of deployment, Sarek packages all required tools using Conda and hosts pre-built container images on DockerHub (https://hub.docker.com). This setup allows users to run Sarek workflows seamlessly within Conda environments or using container technologies such as Docker and Singularity (Kurtzer et al., 2017).

**Resource Requirements**

Efficient execution of Sarek requires computational nodes with the following specifications:
•	At least 20 logical processors on a single physical CPU to handle parallelized tasks.
•	A minimum of 128 GB RAM.
•	At least 4 TB of free storage in the input/output working directory.
Storage usage includes approximately 1.4 TB for output files (e.g., BAM, CRAM, annotated VCF, and CNV files) and 2.3 TB for temporary data, which can be deleted after the workflow is complete unless re-runs from intermediate states are required.

**Computational Efficiency**
Sarek leverages available computational resources effectively by distributing tasks across 20 logical processors. The workflow divides the genome into smaller chunks, each processed as an independent core job, with results merged and sorted in a final step. Many tools utilized within Sarek are inherently parallelized, while others employ a scatter-gather approach to optimize CPU utilization and minimize wall-clock runtime. This design ensures efficient performance on HPC systems.
For this study, the complete Sarek workflow, including preprocessing, variant calling, and annotation, was executed in 7 hours and 58 minutes. The workflow required approximately three times the storage of the original input data, underscoring the need for adequate computational resources.

##Implementation:##

**Variant Calling and Filtering**
Somatic Single-Nucleotide Variants (SNVs) and Insertion/Deletion (Indel) Calling
Somatic SNVs and indels were identified using Mutect2 from the GATK toolkit (v4.3.0.0). Mutect2 was configured in tumor-only mode, utilizing its Bayesian model to distinguish somatic mutations from sequencing artifacts. Default settings were applied unless otherwise specified. Variants passing the default quality control filters (PASS) were retained for downstream analysis. To enhance the reliability of the results, stringent quality thresholds and post-calling filtering were employed to ensure high-confidence variant detection.
The validity of variant calls was evaluated by cross-comparing results with Strelka, SnpEff, and Control-FREECanalyses. These tools were executed within the nf-core/sarek (v3.4.0) pipeline and as standalone applications. This multi-tool validation approach ensured consistency and minimized potential biases in variant detection.
 
**Filtering Criteria**
To ensure high-confidence variants for downstream analysis, filtering was performed according to the following criteria:
1.	PASS Filter: Variants were required to pass the default quality filters applied by the respective variant callers.
2.	Chromosomal Location: Variants were restricted to chromosomes 1–22, X, or Y.
3.	Variant Allele Frequency (VAF): A minimum VAF threshold of 0.05 was applied.
4.	Alternate Allele Depth: Variants called by Mutect2 were required to have a minimum of three alternate allele (ALT) reads.
5.	Mutect2 Exclusivity: Only variants identified by Mutect2 were retained for further analysis.
These filtering steps ensured the inclusion of biologically and clinically relevant variants while reducing false positives.
 
**Variant Annotation**
Variant annotation was performed using SnpEff (v4.3t), which categorizes variants based on their functional impact using a four-tier classification system:
•	HIGH: Variants predicted to have disruptive effects on protein function (e.g., frameshift or nonsense mutations).
•	MODERATE: Variants likely to affect protein effectiveness without complete disruption (e.g., missense mutations).
•	LOW: Variants unlikely to significantly alter protein behavior (e.g., synonymous mutations).
•	MODIFIER: Variants located in non-coding regions or regulatory elements.
SnpEff annotations were used to prioritize variants for downstream functional and clinical interpretation.
 
**Copy Number Variation (CNV) Analysis**
Copy number variants were detected using CNVkit (v0.9.12). The analysis was conducted in tumor-only mode, employing default parameters for segmentation and gene-level copy number estimation. CNVkit leverages read depth and coverage metrics to identify regions of copy number gain or loss.
 
**Mutational Analysis**
Somatic mutations were identified using Mutect2 (GATK v4.3.0.0) and copy number variants (CNVs) were detected using CNVkit (v0.9.12). The identified mutations from Mutect2 were curated and validated by cross-referencing them with regions of copy number gain or loss from CNVkit. This integration ensured that mutations were not only present but also located within genomic regions exhibiting structural alterations, increasing their likelihood of biological relevance.
Validated mutations were further assessed for oncogenic potential using publicly available databases:
•	OncoKB: To determine the oncogenicity and clinical actionability of the mutations.
•	COSMIC: To validate the presence of mutations in cancer-specific datasets.
•	gnomAD: To exclude common germline variants based on population allele frequencies.
•	ClinVar: To identify mutations with known clinical significance.
Mutations classified as oncogenic or likely oncogenic based on these databases were selected for further analysis. Genes harboring these mutations were included in an oncoprint to visualize mutation patterns, co-occurrence, and frequency across samples.


