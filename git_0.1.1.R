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
library(ggthemes)




##
out_dir <- "./BGC/RNAseq/Nip_try2"
dir.create(out_dir, recursive = T)

## load data
all_counts <- fread("C:/your/output/path/all_counts.tsv")

data <- all_counts


## modify file names
colnames(data) <- gsub("*.bam$", "", colnames(data))
colnames(data) <- gsub('/project/egg_rnaseq/nip/outputtry1//', "", colnames(data))

## set name
# HT
setnames(data,"HT-0.5-1_sorted","HT0.5-1")
setnames(data,"HT-0.5-2_sorted","HT0.5-2")


# LT
setnames(data,"LT-0.5-1_sorted","LT0.5-1")
setnames(data,"LT-0.5-2_sorted","LT0.5-2")


###
#
#
countData <- as.matrix(data[, 7: ncol(data)])
rownames(countData) <- data$Geneid


# add group meta 
# Replace with your own experimental conditions
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
# Replace with your own experimental conditions
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



###### 
decorated_gff <- gff_data <- import("C:/path/to/ref/eggnog.emapper.decorated.gff")

genepred_gff <- gff_data <- import("C:/path/to/ref/eggnog.emapper.genepred.gff")



## need edit before load !!
emapper_annotations <- fread("C:/path/to/ref/eggnog.emapper.annotations",  fill=TRUE)


gene2GO <- emapper_annotations[, .(GOs, query)]
gene2Description <- emapper_annotations[, .(GOs, Description)]


### GSEA analysis
# Loop over each condition pair to perform the analysis
# List of condition pairs to analyze
# Replace with your own experimental conditions
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
genome_anno = "C:/path/to/ref/eggnog.anno"
gene_id = "C:/path/to/ref/geneid_match.tab"

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








#####
deepbgc <- fread("C:/path/to/ref/all_gbk_2.bgc.tsv")

genepred_gff_dt <- as.data.table(genepred_gff)

############  gene in BGC 
## According to the start and end of each gene, 
##as well as the start and end of each BGC, map the gene to the corresponding BGC
smlr_rangers <- rbind(genepred_gff_dt[, .(name = ID, chr = seqnames, start = start, end = end, strand)] 
)


# big_rangers : deepbgc
big_rangers <- deepbgc[, .(sequence_id, bgc_candidate_id  , nucl_start , nucl_end, protein_ids)]


## Convert the data.tables to GRanges objects
smlr_granges <- GRanges(seqnames = smlr_rangers$chr,
                        ranges = IRanges(start = smlr_rangers$start,
                                         end = smlr_rangers$end),
                        names = smlr_rangers$name)


big_granges <- GRanges(seqnames = big_rangers$sequence_id,
                       ranges = IRanges(start = big_rangers$nucl_start,
                                        end = big_rangers$nucl_end),
                       names = big_rangers$bgc_candidate_id)


overlaps <- findOverlaps(smlr_granges, big_granges, type = "within")

gene_in_BGC <- data.table(Geneid = smlr_rangers$name[queryHits(overlaps)],
                          BGC_id = big_rangers$bgc_candidate_id[subjectHits(overlaps)],
                          chr = smlr_rangers$chr[queryHits(overlaps)],
                          gene_start = smlr_rangers$start[queryHits(overlaps)], 
                          gene_end = smlr_rangers$end[queryHits(overlaps)], 
                          BGC_start = big_rangers$nucl_start[subjectHits(overlaps)],
                          BGC_end = big_rangers$nucl_end[subjectHits(overlaps)],
                          protein_id = big_rangers$protein_ids[subjectHits(overlaps)])

gene_in_BGC <- unique(gene_in_BGC)

decorated_gff_dt <- as.data.table(decorated_gff)

gene_in_BGC <- merge(gene_in_BGC, decorated_gff_dt[, .(ID, em_target)], by.x = "Geneid", by.y = "ID")






###### 
### geneid mapping to MSU anno
## According to the existing gtf file (note: you need to select the gtf annotation file from the same source as the fasta file), 
## map the newly annotated genes to the existing genes
gtf_data <- import("ref/all.MSU.read.gtf")
gtf_data <- as.data.table(gtf_data)

#
gtf_data <- gtf_data[type == "transcript"]

## 
smlr_rangers <- rbind(genepred_gff_dt[, .(name = ID, chr = seqnames, start = start, end = end, strand)])

# big_rangers : MSU id
big_rangers <- gtf_data[, .(gene_id, seqnames, start, end)]


## Convert the data.tables to GRanges objects
smlr_granges <- GRanges(seqnames = smlr_rangers$chr,
                        ranges = IRanges(start = smlr_rangers$start,
                                         end = smlr_rangers$end),
                        names = smlr_rangers$name)


big_granges <- GRanges(seqnames = big_rangers$seqnames,
                       ranges = IRanges(start = big_rangers$start,
                                        end = big_rangers$end),
                       # strand = big_rangers$strand,
                       names = big_rangers$gene_id)


overlaps <- findOverlaps(smlr_granges, big_granges, type = "any")

gene_in_MSU <- data.table(geneid_new = smlr_rangers$name[queryHits(overlaps)],
                          geneid_MSU = big_rangers$gene_id[subjectHits(overlaps)],
                          chr = smlr_rangers$chr[queryHits(overlaps)],
                          gene_new_start = smlr_rangers$start[queryHits(overlaps)], 
                          gene_new_end = smlr_rangers$end[queryHits(overlaps)], 
                          MSU_start = big_rangers$start[subjectHits(overlaps)],
                          MSU_end = big_rangers$end[subjectHits(overlaps)])

gene_in_MSU <- unique(gene_in_MSU)




###
gene_in_BGC_expr <- merge(gene_in_BGC, resdata)

gene_in_BGC_expr <- gene_in_BGC_expr[, -c("protein_id")]

gene_count <- gene_in_BGC_expr[, .N, by = BGC_id]


# Filter out BGC_id with a count greater than 1
bgc_ids_to_keep <- gene_count[N > 1, BGC_id]

# Delete rows where BGC_id appears only once
gene_in_BGC_expr_filtered <- gene_in_BGC_expr[BGC_id %in% bgc_ids_to_keep]

gene_in_BGC_expr_filtered[, .N, by = BGC_id]




###### define DEBGC
# Get all BGC IDs
BGC_ids <- unique(gene_in_BGC_expr_filtered$BGC_id)

# Create a new data table to save the results
result_table <- data.table(BGC_id = character(),
                           Gene = character(),
                           p_value = numeric(),
                           significance = character())


# Loop through each BGC_id
for (BGC_id_in in BGC_ids) {
  
  print(paste("Processing BGC_id:", BGC_id_in))
  
  # Extract data corresponding to the current BGC_id
  X <- t(gene_in_BGC_expr_filtered[BGC_id == BGC_id_in, 
                                   .(`norm_HT0.5-2`, `norm_HT0.5-3`, `norm_HT1-3`, `norm_HT1-4`, `norm_HT24-2`, `norm_HT24-3`,
                                     `norm_LT0.5-2`, `norm_LT0.5-3`, `norm_LT1-3`, `norm_LT1-4`, `norm_LT24-2`, `norm_LT24-4`)])
  
  # Set the column names of X to be the gene IDs corresponding to the current BGC_id
  colnames(X) <- gene_in_BGC_expr_filtered[BGC_id == BGC_id_in, Geneid]
  
  # Loop through each column of X and perform the GlobalTest
  for (gene_index in 1:ncol(X)) {
    
    # Expression data for the current gene
    Y <- X[, gene_index, drop = FALSE]  # Use drop=FALSE to retain the data frame structure
    
    # Set Y as the response variable, and the other columns of X as covariates
    X_others <- X[, -gene_index, drop = FALSE]  # Select columns other than the current gene as covariates
    
    # Get the current gene's name
    gene_name <- colnames(X)[gene_index]
    
    # Print debug information to check X and Y
    print(paste("Running globaltest for BGC_id:", BGC_id_in, "Gene:", gene_name))
    
    # Run globaltest (assumes X_others are covariates and Y is the response variable)
    result <- gt(Y, X_others)
    
    # Extract p-value
    p_value <- p.value(result)
    
    # Print the current p-value to check if it's valid
    print(paste("p_value for", gene_name, ":", p_value))
    
    # Save the result, adding BGC_id, gene name, and p-value
    if (!is.na(p_value) && p_value < 0.05) {
      significance <- "significant"
    } else {
      significance <- "ns"  # "ns" for non-significant
    }
    
    result_table <- rbind(result_table, 
                          data.table(BGC_id = BGC_id_in, 
                                     Gene = gene_name, 
                                     p_value = p_value, 
                                     significance = significance))
  }
  
  print(paste0("End processing BGC_id:", BGC_id_in))
}

# Print the final result table
print(result_table)

# Get the unique BGC_ids with significant results
unique(result_table[significance == "significant", BGC_id])



######  GENIE3-Based Regulatory Network Analysis for BGCs
library(GENIE3)

BGC_ids <- unique(gene_in_BGC_expr_filtered$BGC_id)

# Initialize a list to store the linkList_filtered for each BGC_id
all_linkLists <- list()

# Function to plot the correlation heatmap for each BGC_id
plot_regulatory_heatmap <- function(BGC_id_in) {
  
  # Get the regulatory links for the selected BGC_id
  linkList_filtered <- all_linkLists[[BGC_id_in]]
  
  # Ensure the filtered links are not empty and have at least two rows
  if (nrow(linkList_filtered) < 2) {
    print(paste("Not enough regulatory links for BGC_id:", BGC_id_in))
    return(NULL)
  }
  
  # Create a correlation matrix based on the regulatory links
  link_matrix <- as.data.table(linkList_filtered)[, .(regulatoryGene, targetGene, weight)]
  
  # Reshape the data into a wide format matrix
  link_matrix_wide <- dcast(link_matrix, regulatoryGene ~ targetGene, value.var = "weight", fill = 0)
  
  
  rownames <- link_matrix_wide$regulatoryGene  # Assign rownames before dropping the column
  link_matrix_wide <- link_matrix_wide[, -1]
  link_matrix_wide <- as.data.frame(link_matrix_wide)
  rownames(link_matrix_wide) <- rownames
  
  # Print to verify the output
  print("Final link_matrix_wide:")
  print(link_matrix_wide)
  rownames(link_matrix_wide)
  
  # Plot heatmap
  pheatmap(as.matrix(link_matrix_wide), 
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           display_numbers = TRUE, 
           show_rownames = TRUE,  # Ensure row names are displayed
           show_colnames = TRUE, 
           angle_col = 45,  # Rotate column names by 45 degrees
           fontsize = 14,  # Increase overall font size
           fontsize_row = 14,  # Increase row name font size
           fontsize_col = 14,  # Increase column name font size
           main = paste("Regulatory Network Heatmap for", BGC_id_in))
}


# Loop through each BGC_id for analysis
for (BGC_id_in in BGC_ids) {
  
  print(paste("Processing BGC_id:", BGC_id_in))
  
  # Extract the gene expression matrix for the current BGC_id
  exprMatrix <- as.matrix(gene_in_BGC_expr_filtered[BGC_id == BGC_id_in, 
                                                    .(`norm_HT0.5-2`, `norm_HT0.5-3`, `norm_HT1-3`, `norm_HT1-4`, `norm_HT24-2`, `norm_HT24-3`,
                                                      `norm_LT0.5-2`, `norm_LT0.5-3`, `norm_LT1-3`, `norm_LT1-4`, `norm_LT24-2`, `norm_LT24-4`)])
  
  # Set row names to be gene IDs
  rownames(exprMatrix) <- gene_in_BGC_expr_filtered[BGC_id == BGC_id_in, Geneid]
  
  # Ensure that the expression matrix is correct
  if (nrow(exprMatrix) == 0 || ncol(exprMatrix) == 0) {
    print(paste("Empty expression matrix for BGC_id:", BGC_id_in))
    next  # Skip if the matrix is empty
  }
  
  # Set the seed to ensure reproducibility
  set.seed(123)  
  
  # Get the number of genes
  num_genes <- nrow(exprMatrix)
  
  # Ensure there are at least two genes to be used as regulators
  if (num_genes < 2) {
    print(paste("Not enough genes for BGC_id:", BGC_id_in))
    next  # Skip this BGC_id if there are not enough genes
  }
  
  # Use all genes as regulators
  regulators <- rownames(exprMatrix)  
  
  # Run GENIE3 to infer the gene regulatory network
  weightMatrix <- GENIE3(exprMatrix, regulators = regulators)
  
  # Extract regulatory links using getLinkList
  linkList_filtered <- getLinkList(weightMatrix, threshold = 0.5)
  
  # Add the filtered linkList for each BGC_id to all_linkLists
  all_linkLists[[BGC_id_in]] <- linkList_filtered
  
  # Check if there are enough regulatory links
  if (nrow(linkList_filtered) < 2) {
    print(paste("Skipping BGC_id:", BGC_id_in, "due to insufficient regulatory links"))
    next  # Skip if there are too few regulatory links
  }
  
  # Print end message for current BGC_id
  print(paste("BGC_id:", BGC_id_in, "down"))
  
  # Plot the regulatory heatmap for the current BGC_id
  plot_regulatory_heatmap(BGC_id_in)
}

# Print all regulatory links for all BGC_ids
print(all_linkLists)





### Visualizing Differential BGC Distribution on Chromosomes
#######  

fxtra <- gene_in_BGC_expr_filtered[ , .(`BGC_id`,`norm_HT0.5-2`, `norm_HT0.5-3`, `norm_HT1-3`, `norm_HT1-4`, `norm_HT24-2`, `norm_HT24-3`,
                                        `norm_LT0.5-2`, `norm_LT0.5-3`, `norm_LT1-3`, `norm_LT1-4`, `norm_LT24-2`, `norm_LT24-4`)]

##
# Extract the numeric portion between "Chr" and "_" in the BGC_id column
fxtra[, chr := as.numeric(gsub("Chr([0-9]+)_.*", "\\1", BGC_id))]

# Print the updated data.table with the new "chr" column
print(fxtra)

fxtra_scale <- fxtra[chr == "1"]

fxtra_scale <- as.data.table(fxtra_scale)

fxtra_scale[, avg_LT := rowMeans(.SD, na.rm = TRUE), .SDcols = c("norm_LT0.5-2", "norm_LT0.5-3", "norm_LT1-3", "norm_LT1-4", "norm_LT24-2", "norm_LT24-4")]
fxtra_scale[, avg_HT := rowMeans(.SD, na.rm = TRUE), .SDcols = c("norm_HT0.5-2", "norm_HT0.5-3", "norm_HT1-3", "norm_HT1-4", "norm_HT24-2", "norm_HT24-3")]



fxtra_scale$Geneid <- gene_in_BGC_expr_filtered$Geneid
fxtra_scale$BGC_id <- gene_in_BGC_expr_filtered$BGC_id


fxtra_scale_avg <- fxtra_scale[, .(
  avg_BGC_LT = sum(avg_LT) / .N,   # Sum avg_LT and divide by the number of Geneids
  avg_BGC_HT = sum(avg_HT) / .N    # Sum avg_HT and divide by the number of Geneids
), by = BGC_id]   # Group by BGC_id


fxtra_scale_avg <- fxtra_scale_avg[, .(avg_BGC_LT = log10(avg_BGC_LT),
                                       avg_BGC_HT = log10(avg_BGC_HT))]


fxtra_scale_avg <- as.data.frame(fxtra_scale_avg)

rownames(fxtra_scale_avg) <- fxtra_scale_avg$BGC_id

fxtra_scale_avg$BGC_id <- NULL

fxtra_scale_avg



set.seed(123)
km.res <- kmeans(scale(fxtra_scale_avg), 3, nstart = 25)

# 3. Visualize
library("factoextra")
fviz_cluster(km.res, data = fxtra_scale_avg,
             palette = c("#1B9E77", "#D95F02", "#7570B3"),
             axes = c(1, 2),
             ggtheme = theme_minimal(),
             main = "Cluster plot for Chr1 BGC",
             pointsize = 1.5,
             labelsize = 16
)



# Compute hierarchical clustering and cut into 4 clusters
res <- hcut(fxtra_scale_avg, k = 5, stand = TRUE)

# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                       "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#33A02C",
                       "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6",
                       "#6A3D9A", "#FFFF99", "#B15928", "#2C7BB6", "#D32F2F",
                       "#6A1B9A", "#F57C00", "#388E3C", "#0288D1", "#FBC02D"))




# Optimal number of clusters for k-means

my_data <- scale(fxtra_scale_avg)
fviz_nbclust(my_data, kmeans, method = "gap_stat")








########  Variance Partitioning
#### try 1

# vp_dt <- gene_in_BGC_expr_filtered[ , .(Geneid, BGC_id,
#                                         `norm_HT0.5-2`, `norm_HT0.5-3`, `norm_HT1-3`, `norm_HT1-4`, `norm_HT24-2`, `norm_HT24-3`,
#                                         `norm_LT0.5-2`, `norm_LT0.5-3`, `norm_LT1-3`, `norm_LT1-4`, `norm_LT24-2`, `norm_LT24-4`)]


vp_dt <- gene_in_BGC_expr_filtered[ , .(Geneid, BGC_id,
                                        `norm_HT0.5-2`, `norm_HT0.5-3`,
                                        `norm_LT0.5-2`, `norm_LT0.5-3`)]

vp_dt$Geneid <- as.factor(vp_dt$Geneid)
vp_dt$BGC_id <- as.factor(vp_dt$BGC_id)


vp_dt

#
vp_df <- as.data.frame(vp_dt)

rownames(gene_in_BGC_expr_filtered) <- gene_in_BGC_expr_filtered$Geneid

vp_df


# turn to long format 
long_dt <- melt(vp_dt, 
                id.vars = c("Geneid", "BGC_id"), 
                measure.vars = patterns("norm_"),
                variable.name = "Condition", 
                value.name = "Expression")

# Make sure Condition is a factor variable
long_dt$Condition <- as.factor(ifelse(grepl("HT", long_dt$Condition), "HT", "LT"))




library(lmerTest)
# Fit the model and output results
model_interaction <- lmer(Expression ~ Condition * BGC_id + (1 | BGC_id:Geneid), data = long_dt)
summary(model_interaction)


# Estimate the model using maximum likelihood method
# Switch to ML (maximum likelihood) estimation method
model_ml <- lmer(Expression ~ Condition * BGC_id + (1 | BGC_id:Geneid), data = long_dt, REML = FALSE)
summary(model_ml)


# Assume p_values is a vector storing the ConditionLT p-values for each BGC_id
p_values <- summary(model_interaction)$coefficients[,"Pr(>|t|)"]

# Apply FDR correction
p_adjusted_FDR <- p.adjust(p_values, method = "fdr")

# Or apply Bonferroni correction
p_adjusted_Bonferroni <- p.adjust(p_values, method = "bonferroni")

# View the adjusted p-values
head(p_adjusted_FDR)
head(p_adjusted_Bonferroni)


# Select significant BGCs (e.g., FDR-adjusted p-values < 0.05)
significant_BGC_FDR <- names(p_adjusted_FDR)[p_adjusted_FDR < 0.05]
significant_BGC_Bonferroni <- names(p_adjusted_Bonferroni)[p_adjusted_Bonferroni < 0.05]

# View the differentially expressed BGCs
significant_BGC_FDR
significant_BGC_Bonferroni





#####
#### GO anno to BGC anno
BGC_anno <- emapper_annotations[, .(query, Description, GOs)]

## ----TERNAME (GOid, Term, Ontology)

TER2NAME_BGC <- gene_in_BGC[, .(BGC_id, BGC_id)]


## ----TERM2GENE (GOid, gene)

TERM2GENE_BGC <- gene_in_BGC[, .(BGC_id, Geneid)]


##
ego <- enricher(gene = DEgene$Geneid, 
                TERM2GENE = TERM2GENE_BGC, 
                TERM2NAME = TER2NAME_BGC,
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,)

as.data.table(ego)
barplot(ego)




####### 
##classify product
gene_in_BGC_product <- merge(gene_in_BGC_expr_filtered, deepbgc[, .(bgc_candidate_id, product_class)],
                             by.x = "BGC_id", by.y = "bgc_candidate_id", all.x =T)


head(gene_in_BGC_product)

gene_in_BGC_product_sum <- gene_in_BGC_product[, .N, by = product_class]


## product barplot
gene_in_BGC_product_sum$text_position <- ifelse(gene_in_BGC_product_sum$N > 100, "inside", "outside")

gene_in_BGC_product_sum <- gene_in_BGC_product_sum[!is.na(product_class) & product_class != ""]



ggplot(data = gene_in_BGC_product_sum, aes(x = product_class, y = N, fill = product_class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  geom_text(aes(
    label = N,
    hjust = ifelse(text_position == "inside", 1.2, -0.2)
  ), color = ifelse(gene_in_BGC_product_sum$text_position == "inside", "white", "black"), size = 6) +  # Increased text size
  scale_fill_viridis_d(option = "G", end = 0.9) +
  labs(
    x = "Product Class",
    y = "Count of Product Class",
    fill = "Product Class"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Enlarge x-axis text
    axis.text.y = element_text(size = 14),  # Enlarge y-axis text
    axis.title.x = element_text(size = 16),  # Enlarge x-axis label
    axis.title.y = element_text(size = 16),  # Enlarge y-axis label
    legend.title = element_text(size = 16),  # Enlarge legend title
    legend.text = element_text(size = 14),  # Enlarge legend text
    plot.title = element_text(size = 18, hjust = 0.5),  # Enlarge plot title
    legend.position = "top"
  )





########
### Gene expression changes within BGCs
### Use gggenes to display gene arrangement within BGC

gene_in_BGC_expr_filtered

gene_in_BGC_expr_arrangement <- merge(gene_in_BGC_expr_filtered, genepred_gff_dt[, .( width, strand, em_target )],
                                      by = "em_target", all.x = T)



#### loop
out_dir <- "BGC/heatmap/"
# dir.create(out_dir)


# 
gene_in_BGC_expr_arrangement <- as.data.table(gene_in_BGC_expr_arrangement)

# get all  BGC_id
BGC_ids <- unique(gene_in_BGC_expr_arrangement$BGC_id)


# generate heatmap for each BGC_id
for (BGC_id_in in BGC_ids) {
  # get a BGC_id 
  BGC_data <- gene_in_BGC_expr_arrangement[BGC_id == BGC_id_in]
  
  print(BGC_id_in)
  
  # get expression col
  selected_data <- BGC_data[, .(Geneid, norm_Control1, norm_Control2, norm_treat1, norm_treat2)]
  
  # convert to long
  long_data <- pivot_longer(selected_data, 
                            cols = starts_with("norm_"),
                            names_to = "Condition",
                            values_to = "Expression")
  
  # convert to matrix
  long_data <- unique(long_data)
  
  heatmap_data <- reshape2::dcast(long_data, Geneid ~ Condition, value.var = "Expression")
  
  # check NA
  heatmap_data[is.na(heatmap_data)] <- 0  #  NA to 0
  heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]  # delete col include NA
  
  print( heatmap_data$Geneid)
  
  
  # row normalization
  normalized_data <- t(apply(heatmap_data[-1], 1, function(x) {
    if (max(x) - min(x) == 0) {
      return(rep(0, length(x)))  
    } else {
      return((x - min(x)) / (max(x) - min(x)))  
    }
  }))
  
  # set colnames and rownames
  # print(heatmap_data$Geneid)
  
  rownames(normalized_data) <- heatmap_data$Geneid
  colnames(normalized_data) <- colnames(heatmap_data)[2:ncol(heatmap_data)]
  
  print(colnames(normalized_data) )
  print(rownames(normalized_data) )
  
  # Check if there are duplicate values 
  if (any(duplicated(normalized_data))) {
    normalized_data <- normalized_data + matrix(runif(length(normalized_data), 0, 1e-10),
                                                nrow=nrow(normalized_data), ncol=ncol(normalized_data))
  }
  
  # Transpose the data
  normalized_data_transposed <- t(normalized_data)
  
  # 
  p_heatmap <- pheatmap(normalized_data_transposed, 
                        cluster_rows = FALSE, 
                        cluster_cols = FALSE, 
                        display_numbers = FALSE, 
                        fontsize_number = 10,
                        show_rownames = TRUE, 
                        show_colnames = TRUE,
                        main = paste("Heatmap of Gene Expressions for BGC:", BGC_id_in),
                        cellwidth = 15,  
                        cellheight = 15,  
                        fontsize_col = 10,  
                        angle_col = 45, 
                        silent = TRUE)  
  
  # gggenes
  plot_layout <- ggplot(BGC_data) + 
    geom_gene_arrow(
      aes(
        xmin = gene_start, 
        xmax = gene_end, 
        y = seq_along(em_target), 
        fill = strand
      ),
      arrowhead_width = unit(3, "mm"),
      arrowhead_height = unit(3, "mm"),
      arrow_body_height = unit(2, "mm")
    ) +
    scale_y_continuous(breaks = seq(1, nrow(BGC_data), by = 1)) +  
    labs(y = "Genes") + 
    theme_void() +
    theme(legend.position = "none")
  
  
  # Create your combined plot
  combined_plot <- grid.arrange(
    grobs = list(ggplotGrob(plot_layout), p_heatmap$gtable),
    nrow = 2,
    heights = c(1/3, 2/3)
  )
  
  # Set dimensions for saving
  custom_width <- 12  # Width in inches
  custom_height <- 8  # Height in inches
  
  # Save the combined plot
  ggsave(paste0(out_dir, "/", BGC_id_in, "combined_heatmap.png"), 
         plot = combined_plot, width = custom_width, height = custom_height, dpi = 300)
  
  
  
  # heatmap_data <- heatmap_data[0]
  
  
  
}





####### 
# Time-series analysis
# Merge the data
gene_in_BGC_expr <- merge(gene_in_BGC, resdata)

# Remove 'protein_id' column
gene_in_BGC_expr <- gene_in_BGC_expr[, !"protein_id"]

# Filter BGC_ids with more than one entry and calculate mean log2FoldChange
gene_in_BGC_expr_filtered <- gene_in_BGC_expr[, .(
  mean_log2FoldChange_HT1_LT1 = mean(HT1_LT1_log2FoldChange, na.rm = TRUE),
  mean_log2FoldChange_HT0.5_LT0.5 = mean(HT0.5_LT0.5_log2FoldChange, na.rm = TRUE),
  mean_log2FoldChange_HT24_LT24 = mean(HT24_LT24_log2FoldChange, na.rm = TRUE)
), by = BGC_id][abs(mean_log2FoldChange_HT1_LT1) > 1 | 
                  abs(mean_log2FoldChange_HT0.5_LT0.5) > 1 | 
                  abs(mean_log2FoldChange_HT24_LT24) > 1]

# Merge the filtered data with expression values
gene_in_BGC_expr_filtered <- merge(gene_in_BGC_expr_filtered, gene_in_BGC_expr, by = "BGC_id")

# Melt the data to long format
tidy_dt <- melt(
  gene_in_BGC_expr_filtered,
  id.vars = c("Geneid", "BGC_id"),
  measure.vars = patterns("^norm_"),
  variable.name = "condition",
  value.name = "expression"
)

# Extract the time from the condition names and convert to date
tidy_dt[, time := paste0("2024-1-", gsub(".*_(HT|LT)(\\d+(\\.\\d+)?)-.*", "\\2", condition))]
tidy_dt[, date := as.Date(time, format = "%Y-%m-%d")]



# Apply tq_transmute for each BGC_id
result_list <- split_by_bgc %>%
  map(~ .x %>%
        tq_transmute(
          select = expression,        # Column to apply the function to
          mutate_fun = apply.daily,   # Apply daily aggregation
          FUN = colMeans,             # Function to calculate means
          na.rm = TRUE,               # Remove NA values
          col_rename = "mean_expression" # Rename the output column
        ) %>%
        mutate(BGC_id = unique(.x$BGC_id))  # Add the BGC_id back to the result
  )



# Convert to data.table
mean_tidyverse_downloads_w <- as.data.table(mean_tidyverse_downloads_w)

# Reorder 'date' column using data.table syntax (keeping factor order intact)
mean_tidyverse_downloads_w[, date := factor(date, levels = c("2024-01-05", "2024-01-01", "2024-01-24"))]

# Select 10 BGC_ids
selected_bgc_ids <- unique(mean_tidyverse_downloads_w$BGC_id)[c(70, 71, 72, 80, 81, 82, 86, 40, 68)]


# Filter the data to include only the first 10 BGC_ids
filtered_data <- mean_tidyverse_downloads_w[BGC_id %in% selected_bgc_ids]

filtered_data[, date := factor(format(as.Date(date), "%d"), levels = c("05", "01", "24"))]

# View the updated data
filtered_data


# Plot using ggplot2 with facet_wrap and enlarged font size
ggplot(filtered_data, aes(x = date, y = mean_expression, color = BGC_id)) +
  geom_point() +  # Plot points
  geom_smooth(aes(x = as.numeric(date)), method = "loess", se = FALSE) +  # Smooth curve
  expand_limits(y = 0) +
  scale_color_tq() +
  scale_x_discrete(limits = c("05", "01", "24")) +  # Correct order of dates
  theme_tq() +
  theme(
    legend.position = "none", 
    axis.title.x = element_blank(),  # Axis titles
    axis.title.y = element_text(size = 16), 
    axis.text = element_text(size = 16),   # Axis tick labels
    strip.text = element_text(size = 16),  # Facet labels
    plot.title = element_text(size = 16),  # Plot title (if you add one)
    legend.text = element_text(size = 16), # Legend text (if used)
    legend.title = element_text(size = 16) # Legend title (if used)
  ) +
  facet_wrap(~ BGC_id, ncol = 3, scales = "free_y")  # Facet by BGC_id















