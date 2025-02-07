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


```
#
countData <- as.matrix(data[, 7: ncol(data)])
rownames(countData) <- data$Geneid


# add group meta 
database <- data.frame(name = colnames(data)[7: ncol(data)], 
                       condition=c(rep(c("HT0.5","HT1","HT24","LT0.5","LT1","LT24"), each = 2)))

row.names(database) <- colnames(data)[7: ncol(data)]

database$condition <- as.factor(database$condition)

## select 
# https://www.biostars.org/p/336298/
dds <- DESeqDataSetFromMatrix(countData, colData=database, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- as.data.frame(results(dds))


##
# PCA
vsd <- vst(dds, blind=FALSE)
p_pca <- plotPCA(vsd, intgroup=c("condition")) + theme_bw()
# ggsave(paste0(out_dir, "/PCA_all.png"), p_pca, width = 5, height = 5)



###
##
compa_pair <- combn(c("HT0.5","HT1","HT24","LT0.5","LT1","LT24"), 2, simplify = F)
de_sum <- function(dds, cond_pair){
  #
  cond_1 <- cond_pair[1]
  cond_2 <- cond_pair[2]
  
  #
  result_dt <- as.data.frame(results(dds, contrast = c("condition", cond_1, cond_2)))
  
  #
  setDT(result_dt, keep.rownames = "Geneid")
  result_dt <- result_dt[, .(Geneid, log2FoldChange, pvalue, padj)]
  colnames(result_dt) <- paste0(cond_1, "_", cond_2, "_", colnames(result_dt))
  colnames(result_dt)[1] <- "Geneid"
  
  #
  return(result_dt)
}
res <- Reduce(function(x, y) merge(x, y, all=TRUE), lapply(compa_pair, function(x) de_sum(dds, x)))

##
norm_data <- as.data.table(counts(dds, normalized=TRUE), keep.rownames = "Geneid")
colnames(norm_data)[2:ncol(norm_data)] <- paste0("norm_", colnames(norm_data)[2:ncol(norm_data)])


## merge result and normalized counts
resdata <- merge(res, norm_data, all.x = T, by = "Geneid")

```

## 7.Gene Set Enrichment Analysis (GSEA)  and GO Analysis
This part performs Gene Set Enrichment Analysis (GSEA)  and GO Analysis to identify significant biological processes or pathways associated with different experimental conditions. It processes gene expression data (resdata) by filtering differentially expressed genes (DE genes) based on a log2FoldChange threshold and p-value. For each condition pair (e.g., HT0.5 vs. LT0.5), the script dynamically selects the relevant log2FoldChange and p-value columns and uses these to filter DE genes. The analysis is looped over multiple condition pairs (e.g., HT0.5_LT0.5, HT1_LT1), producing enrichment results for each.

Set the condition pairs in the condition_pairs variable and execute the script to perform GSEA and visualize the results. 

```
decorated_gff <- gff_data <- import("yourpath/eggnog.emapper.decorated.gff")

genepred_gff <- gff_data <- import("yourpath/eggnog.emapper.genepred.gff")

## need edit before load !!
emapper_annotations <- fread("yourpath/eggnog.emapper.annotations",  fill=TRUE)


gene2GO <- emapper_annotations[, .(GOs, query)]
gene2Description <- emapper_annotations[, .(GOs, Description)]


### GSEA analysis
# Loop over each condition pair to perform the analysis
# List of condition pairs to analyze
condition_pairs <- list(
  c("HT0.5", "LT0.5"),
  c("HT1", "LT1"),
  c("HT24", "LT24")
)

for (pair in condition_pairs) {
  # Generate the log2FoldChange column name dynamically
  condition_col <- paste0(pair[1], "_", pair[2], "_log2FoldChange")
  
  # Debugging: Print column names in resdata
  print("Available columns in resdata:")
  print(condition_col)
  
  # Prepare GSEA data by selecting Geneid and the dynamic condition column using .SD
  gsea_dt <- resdata[, .SD, .SDcols = c("Geneid", condition_col)]
  
  # Debugging: Check the first few rows of gsea_dt
  print(paste("Checking data for condition:", condition_col))
  print(head(gsea_dt))
  
  # Rank the values in the specified condition column
  gsea_dt[, rank := rank(get(condition_col), ties.method = "random")]
  
  # Order by rank
  gsea_dt <- gsea_dt[order(-rank)]
  
  # Create the ranked list for GSEA
  glist <- gsea_dt$rank
  names(glist) <- gsea_dt$Geneid
  
  # Sort the list in decreasing order
  glist_sorted <- sort(glist, decreasing = TRUE)
  
  # Perform GSEA
  ego <- GSEA(
    gene = glist_sorted, 
    TERM2GENE = gene2GO, 
    TERM2NAME = gene2Description, 
    scoreType = "pos"
  )
  
  # Visualize the results
  enrichplot::dotplot(ego, showCategory = 10, font.size = 14)
}



######
genome_anno = "yourpath/eggnog.anno"
gene_id = "yourpath/geneid_match.tab"

egg <- fread(genome_anno, sep  ="\t")

idmatch <- fread(gene_id, header = F,  col.names = c("query"))
egg <- merge(egg, idmatch, all.x = T)


anno_raw <- merge(egg[, .(query, Preferred_name, Description, COG_category)], 
                  data[, c(1, 7:ncol(data)), with = FALSE],
                  by.x = "query", by.y = "Geneid",
                  all.y = T)
anno_all <- merge(anno_raw, norm_data,
                  by.x = "query", by.y = "Geneid",all.x = T)


## 
# prepare background GO
eggGO <- egg[GOs != '-', .(query, GOs)]
eggGO <- eggGO[, .(GOs = unlist(tstrsplit(GOs, ",", type.convert = TRUE))), by = "query"]
colnames(eggGO) <- c("gene", "go_id")
term_name <- go2term(eggGO$go_id)
ont <- go2ont(eggGO$go_id)
eggGO <- merge(eggGO, ont, by = "go_id")
term_name <- merge(term_name, ont, by = "go_id")

# Loop over each condition pair to perform the analysis
for (pair in condition_pairs) {
  # Generate the log2FoldChange and pvalue column names dynamically
  condition_col <- paste0(pair[1], "_", pair[2], "_log2FoldChange")
  condition_col_pvalue <- paste0(pair[1], "_", pair[2], "_pvalue")
  
  # Ensure the relevant columns are numeric
  resdata[[condition_col]] <- as.numeric(resdata[[condition_col]])
  resdata[[condition_col_pvalue]] <- as.numeric(resdata[[condition_col_pvalue]])
  
  # Filter DE genes based on log2FoldChange and p-value
  resdata_DE <- resdata[abs(resdata[[condition_col]]) >= 1 & resdata[[condition_col_pvalue]] <= 0.05]
  
  # Debugging: Check if there are any DE genes
  print(paste("Number of DE genes for", condition_col, ":", nrow(resdata_DE)))
  
  if (nrow(resdata_DE) == 0) {
    next  # Skip if no DE genes for the current condition pair
  }
  
  # Perform GSEA
  ego <- enricher(
    gene = resdata_DE$Geneid, 
    TERM2GENE = eggGO[, .(go_id, gene)], 
    TERM2NAME = term_name
  )
  
  # Convert result to data.table and plot
  ego_dt <- as.data.table(ego)
  barplot(ego)
  
  # Dotplot visualization
  enrichplot::dotplot(ego, showCategory = 10, font.size = 14)
  
  # Pairwise term similarity for enrichment result
  ego <- pairwise_termsim(ego)
  
  # EMAP plot
  emapplot(ego, showCategory = 10)
}

```



## Principle of Globaltest
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


## GENIE3-Based Regulatory Network Analysis for BGCs
Secondary metabolite biosynthetic gene clusters (BGCs) play a crucial role in producing bioactive compounds, including antibiotics, anticancer agents, and plant defense molecules. However, understanding how genes within a BGC are regulated remains a challenge. This project utilizes GENIE3, a machine-learning-based approach, to infer gene regulatory networks within BGCs based on expression data.

By applying GENIE3, we can systematically predict key regulatory genes within a BGC, revealing potential transcriptional control mechanisms. This approach allows researchers to identify critical regulatory hubs that may govern the production of secondary metabolites, providing a foundation for rational metabolic engineering and functional validation.

Why GENIE3?
GENIE3 reconstructs regulatory networks by predicting directed interactions between genes using a Random Forest-based feature selection method. Unlike traditional correlation-based methods, which only capture general co-expression patterns, GENIE3 identifies potential causal relationships between transcription factors and their target genes.






