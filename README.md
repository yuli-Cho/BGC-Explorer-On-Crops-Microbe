# x-BGC: Comprehensive BGC Gene Analysis on crops & Microbe

## Overview  
**x-BGC** is an innovative bioinformatics pipeline that enables researchers to analyze **biosynthetic gene clusters (BGCs)** across a wide range of species in a flexible and reference-free approach. Whether you are studying plants, microbes, or fungi, x-BGC allows you to predict and analyze BGC dynamics with reference genome and transcriptome files. Thus, xPAN-BGC can provide a solution for analyzing non-model organisms with incomplete or no genome annotations with dynamic time series, allowing you to track the dynamics of BGCs over time.  

**Origin of the name x-BGC**  
X: an abstract concept symbolizing flexibility – indicating that x-BGC can be run without solid reference annotation and can be dynamically analyzed at different time points, providing continuous insights into their behavior  


## Quick Start  
The whole process is divided into three parts:  

1. RNASEQ
---shell
bash  bgc_rnaseq.sh $1 $2 ./dir/to/output/
----

2. FASDFAS
---shell
deepbgc / / 
---

3. R for de-BGC (differ) and dynamic
avalibale in `git_example.R` 

## implementary
**The idea of ​​the shell part:** the core input is fasta file, and the output is gff file predicted and edited by eggnog.  
**The operation of the shell part:** Please run the `bgc_rnaseq.sh` script. Before using it, please make sure that the necessary software has been installed: `deepbgc, eggNOG, Python, fastp, hisat2, samtools, featureCounts, eggNOG-mapper, seqkit`, and `build_em_genbank.py` needs to be saved before running.  

**The operation of the R script part:** the core input is `all_counts.tsv`, a series of files generated by eggnog, and `bgc.tsv` generated by deepbgc containing the position information of each BGC  

 

## Example Function and Output  
## 1.Gene Set Enrichment Analysis (GSEA)  and GO Analysis  
This part performs Gene Set Enrichment Analysis (GSEA)  and GO Analysis to identify significant biological processes or pathways associated with different experimental conditions. It processes gene expression data (resdata) by filtering differentially expressed genes (DE genes) based on a log2FoldChange threshold and p-value. For each condition pair (e.g., HT0.5 vs. LT0.5), the script dynamically selects the relevant log2FoldChange and p-value columns and uses these to filter DE genes. The analysis is looped over multiple condition pairs (e.g., HT0.5_LT0.5, HT1_LT1), producing enrichment results for each.  

Set the condition pairs in the condition_pairs variable and execute the script to perform GSEA and visualize the results.   
![GSEAplot](https://github.com/user-attachments/assets/f766f5da-8d79-4c69-8a3c-7918a6d2740e)



## 2.Use BGC annotation to enrich differential expression BGC  
Use deepbgc to obtain the function-related annotations of each BGC. These annotations are used as a background library, similar to the annotation background in GO analysis. On the other hand, through the expression data, the BGCs with differential changes are found. The purpose is to find BGCs that have significant differences under specific conditions and will be enriched.  
![image](https://github.com/user-attachments/assets/0385c23a-c8fe-48b0-8fd1-7a72a528c78d)



## 3.Principle of Globaltest  
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


## 4.GENIE3-Based Regulatory Network Analysis for BGCs  
Secondary metabolite biosynthetic gene clusters (BGCs) play a crucial role in producing bioactive compounds, including antibiotics, anticancer agents, and plant defense molecules. However, understanding how genes within a BGC are regulated remains a challenge. This project utilizes GENIE3, a machine-learning-based approach, to infer gene regulatory networks within BGCs based on expression data.  

By applying GENIE3, we can systematically predict key regulatory genes within a BGC, revealing potential transcriptional control mechanisms. This approach allows researchers to identify critical regulatory hubs that may govern the production of secondary metabolites, providing a foundation for rational metabolic engineering and functional validation.  

Why GENIE3?  
GENIE3 reconstructs regulatory networks by predicting directed interactions between genes using a Random Forest-based feature selection method. Unlike traditional correlation-based methods, which only capture general co-expression patterns, GENIE3 identifies potential causal relationships between transcription factors and their target genes.  

![image](https://github.com/user-attachments/assets/2c328bb3-5814-46bc-92e1-7b664cabcce2)


## 5.Visualizing Differential BGC Distribution on Chromosomes  
we leverage the factoextra package in R to visualize the distribution of Biosynthetic Gene Clusters (BGCs) that exhibit differential patterns across various chromosomes. The factoextra package, known for enhancing the visualization of multivariate data analyses such as PCA and clustering methods,we aim to gain insights into the genomic landscape of BGCs under specific conditions.  
![image](https://github.com/user-attachments/assets/a26005e3-b647-4d01-977d-1768865bf2ff)


## 6.Gene Product Class Distribution  
Biosynthetic gene clusters (BGCs) are responsible for the production of a variety of secondary metabolites that can include polyketides, terpenes, non-ribosomal peptides (NRPs), and others. This part of analyze the distribution of product classes from gene clusters in a genomic dataset. The primary goal of this analysis is to visualize the abundance of different biosynthetic products across the gene clusters, particularly focusing on understanding the diversity of biosynthetic gene clusters (BGCs) involved in various metabolic pathways.  


![image](https://github.com/user-attachments/assets/80e4c48f-a0f6-4a56-8833-0035b784cc8f)

## 7.BGC Time-series analysis  
We also pay attention to the changes in BGC over time.  
![image](https://github.com/user-attachments/assets/3a267266-2f56-406b-a60a-f2066f2e9319)





