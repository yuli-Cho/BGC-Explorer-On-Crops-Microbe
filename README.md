### BGC-Explorer: Comprehensive BGC Gene Analysis on crops & Microbe

## Overview
**BGC-Explorer** is a bioinformatics pipeline designed to perform genome-wide analysis of **biosynthetic gene clusters (BGCs)** using the **deepbgc** tool. The project integrates genome-wide annotation with **eggNOG** and combines publicly available transcriptome data to explore the functional response of BGCs to various treatments.

## Features
**Genome-wide BGC Identification:** Utilizes deepbgc for large-scale identification of biosynthetic gene clusters (BGCs) from genomic sequences.\
**Genome Annotation:** Leverages eggNOG for comprehensive genome-wide functional annotation.\
**Transcriptomic Analysis:** Combines public transcriptome data to explore differential expression of BGCs under specific treatments.\
**Gene Interaction Analysis:** Investigates the potential interactions between BGCs, exploring whether they influence each other.\
**Gene Family Exploration:** Identifies key genes within BGC families and assesses their importance in regulating biosynthesis.\

## Getting Started
# Prerequisites
Before running the pipeline, ensure you have the following installed:\
**deepbgc**\
**eggNOG**\
**Python 3.6**\
**fastp** (for quality control)\
**hisat2** (for alignment)\
**samtools** (for BAM file processing)\
**featureCounts** (for gene counting)\
**eggNOG-mapper** (for functional annotation)\
**seqkit** (for sequence manipulation)\

**##Usage**
# In shell
1. Please download the eggNOG database before you start and save it in the data folder under the eggnog-mapper software directory.\
`wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz`\
`wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz`\
`wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz`\
`wget http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz`\
`wget http://eggnog5.embl.de/download/emapperdb-5.0.2/pfam.tar.gz`

2. Please download genome.fa file and species transcriptome data first\


**1. Functional Annotation and RNA-seq analysis:**

This part performs the following steps:

Functional Annotation using eggNOG.\
Genome Indexing with hisat2.\
Data Preprocessing with fastp for quality control and hisat2 for RNA-Seq read alignment.\
Gene Expression Quantification using featureCounts.  

**# Input**  
**$1:** Directory containing sample fastq files (paired-end fastq.gz files).\
**$2:** Directory where output files will be stored.\
**$3:** Reference genome in .fna format (FASTA format).\
**$4:** Directory containing eggNOG reference files.  


**# Output**  
**$2/anno/:** Directory containing annotation results, including functional annotation and CDS prediction.\
**$2/hisat2_index/:** Directory containing genome index files for hisat2.\
**$2/<sample>/:** Directories for each sample containing cleaned fastq files and corresponding aligned BAM files.\
**$2/all_counts.tsv:** Gene expression counts for all samples, generated by featureCounts.\

**Activate Conda Environment:** Ensure you have the necessary dependencies installed in your conda environment:  
`conda activate RNAseq`

**Run the Pipeline:** Run the script with the required arguments:  
`nohup bash ./rna_seq_pipeline.sh <input_directory> <output_directory> <genome_fna> <eggnog_reference_dir> > /log_files/rna_seq_pipeline.out 2>&1 &`



# In Rstudio
Before running the pipeline, you need to install and load the following packages：  

`
library(clusterProfiler)  
library(enrichplot)  
library(ggplot2)  
library(data.table)  
library(GO.db)  
library(DESeq2)  
library(pathview)  
library(GOSemSim)  
library(AnnotationDbi)  
library(org.Osativa.eg.db)  
library(rtracklayer)  
library(Rsubread)  
library(pheatmap)  
library(RColorBrewer)  
library(geneplotter)  
library(gplots)  
library(jsonlite)  
library(purrr)  
library(RCurl)  
library(stringr)  
library(dplyr)  
library(DOSE)  
library(enrichplot)  
library(globaltest)  
library(factoextra)  
library(FactoMineR)  
library(variancePartition)  
library(EnhancedVolcano)  
library(tidyr)  
library(pheatmap)  
library(gridExtra)  
library(lme4)  
library(variancePartition)  
library(Matrix)  
library(colorRamps)  
library(lmerTest)  
library(tibble)  
library(tidyverse)  
library(tidyquant)  
`



