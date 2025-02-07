# BGC-Explorer: Comprehensive BGC Gene Analysis on crops & Microbe

## Overview
**BGC-Explorer** is a bioinformatics pipeline designed to perform genome-wide analysis of **biosynthetic gene clusters (BGCs)** using the **deepbgc** tool. The project integrates genome-wide annotation with **eggNOG** and combines publicly available transcriptome data to explore the functional response of BGCs to various treatments.

## Features
**Genome-wide BGC Identification:** Utilizes deepbgc for large-scale identification of biosynthetic gene clusters (BGCs) from genomic sequences.  
**Genome Annotation:** Leverages eggNOG for comprehensive genome-wide functional annotation.  
**Transcriptomic Analysis:** Combines public transcriptome data to explore differential expression of BGCs under specific treatments.  
**Gene Interaction Analysis:** Investigates the potential interactions between BGCs, exploring whether they influence each other.  
**Gene Family Exploration:** Identifies key genes within BGC families and assesses their importance in regulating biosynthesis.  




**1. Functional Annotation and RNA-seq analysis:**

This part performs the following steps:

Functional Annotation using eggNOG.  
Genome Indexing with hisat2.  
Data Preprocessing with fastp for quality control and hisat2 for RNA-Seq read alignment.  
Gene Expression Quantification using featureCounts.  








## 1.Data Preparation:  
· Count data is extracted from the input dataset (data), with genes as rows and samples as columns.  
· A metadata table (database) is created to define the experimental conditions (e.g., HT0.5, HT1, LT0.5, etc.).  

## 2.DESeq2 Analysis:  
· A DESeqDataSet object is created from the count data and metadata.  
· Differential expression analysis is performed, with low-count genes filtered out, followed by running DESeq.  

## 3.Principal Component Analysis (PCA):  
· Variance Stabilizing Transformation (VST) is applied to the DESeq2 object, followed by PCA for visualizing the relationship between samples based on gene expression profiles.  

## 4.Pairwise Comparisons:
· Pairwise comparisons between conditions (e.g., HT0.5 vs HT1, HT24 vs LT0.5, etc.) are conducted, and results such as log2 fold change, p-values, and adjusted p-values are stored.

## 5.Normalized Counts:
· Normalized count data is extracted and merged with the differential expression results.

## 6.Output:
· The final results, including normalized counts and differential expression statistics, are saved in a comprehensive data table (resdata).




## 7.Gene Set Enrichment Analysis (GSEA)  and GO Analysis
This part performs Gene Set Enrichment Analysis (GSEA)  and GO Analysis to identify significant biological processes or pathways associated with different experimental conditions. It processes gene expression data (resdata) by filtering differentially expressed genes (DE genes) based on a log2FoldChange threshold and p-value. For each condition pair (e.g., HT0.5 vs. LT0.5), the script dynamically selects the relevant log2FoldChange and p-value columns and uses these to filter DE genes. The analysis is looped over multiple condition pairs (e.g., HT0.5_LT0.5, HT1_LT1), producing enrichment results for each.

Set the condition pairs in the condition_pairs variable and execute the script to perform GSEA and visualize the results. 



## 8.Use BGC annotation to enrich differential expression BGC
Use deepbgc to obtain the function-related annotations of each BGC. These annotations are used as a background library, similar to the annotation background in GO analysis. On the other hand, through the expression data, the BGCs with differential changes are found. The purpose is to find BGCs that have significant differences under specific conditions and will be enriched.
![image](https://github.com/user-attachments/assets/0385c23a-c8fe-48b0-8fd1-7a72a528c78d)



## 9.Principle of Globaltest
Globaltest is a gene set analysis method that evaluates whether a predefined set of genes is associated with a particular phenotype or condition. Unlike traditional differential expression analysis, which examines individual genes, Globaltest assesses the collective behavior of a gene set, such as genes within a biosynthetic gene cluster (BGC).

The method is based on a random effects model, where gene expression values within a set are regressed against the phenotype of interest. Instead of testing each gene separately, Globaltest considers the correlation structure among genes and determines whether the variation in gene expression is significantly linked to the condition. The null hypothesis assumes that the genes in the set have no systematic association with the phenotype, and the statistical significance is determined using likelihood ratio tests.

From a biological perspective, Globaltest is particularly useful for identifying coordinated gene regulation in pathways or clusters, such as BGCs in microbial or plant genomes. Since BGCs encode enzymes responsible for producing secondary metabolites, their expression is often co-regulated under specific conditions. By applying Globaltest, researchers can detect whether the entire cluster is differentially regulated after treatment, rather than relying on the expression changes of individual genes. This holistic approach provides insights into the functional relevance of BGCs under different environmental or experimental conditions.


```
> result_table
                          BGC_id        Gene      p_value significance
                          <char>      <char>        <num>       <char>
    1: Chr10_17545673-17877467.1     Chr10_0 0.7517224700           ns
    2: Chr10_17545673-17877467.1 Chr10_10018 0.6547764784           ns
    3: Chr10_17545673-17877467.1 Chr10_10333 0.0767277593           ns
    4: Chr10_17545673-17877467.1 Chr10_10405 0.7172450233           ns
    5: Chr10_17545673-17877467.1 Chr10_10612 0.0002512284  significant
```


## 10.GENIE3-Based Regulatory Network Analysis for BGCs
Secondary metabolite biosynthetic gene clusters (BGCs) play a crucial role in producing bioactive compounds, including antibiotics, anticancer agents, and plant defense molecules. However, understanding how genes within a BGC are regulated remains a challenge. This project utilizes GENIE3, a machine-learning-based approach, to infer gene regulatory networks within BGCs based on expression data.

By applying GENIE3, we can systematically predict key regulatory genes within a BGC, revealing potential transcriptional control mechanisms. This approach allows researchers to identify critical regulatory hubs that may govern the production of secondary metabolites, providing a foundation for rational metabolic engineering and functional validation.

Why GENIE3?
GENIE3 reconstructs regulatory networks by predicting directed interactions between genes using a Random Forest-based feature selection method. Unlike traditional correlation-based methods, which only capture general co-expression patterns, GENIE3 identifies potential causal relationships between transcription factors and their target genes.

![image](https://github.com/user-attachments/assets/2c328bb3-5814-46bc-92e1-7b664cabcce2)


## 11.Visualizing Differential BGC Distribution on Chromosomes
we leverage the factoextra package in R to visualize the distribution of Biosynthetic Gene Clusters (BGCs) that exhibit differential patterns across various chromosomes. The factoextra package, known for enhancing the visualization of multivariate data analyses such as PCA and clustering methods,we aim to gain insights into the genomic landscape of BGCs under specific conditions. 
![image](https://github.com/user-attachments/assets/a26005e3-b647-4d01-977d-1768865bf2ff)


## 12.Gene Product Class Distribution  
Biosynthetic gene clusters (BGCs) are responsible for the production of a variety of secondary metabolites that can include polyketides, terpenes, non-ribosomal peptides (NRPs), and others. This part of analyze the distribution of product classes from gene clusters in a genomic dataset. The primary goal of this analysis is to visualize the abundance of different biosynthetic products across the gene clusters, particularly focusing on understanding the diversity of biosynthetic gene clusters (BGCs) involved in various metabolic pathways.  


![image](https://github.com/user-attachments/assets/80e4c48f-a0f6-4a56-8833-0035b784cc8f)




