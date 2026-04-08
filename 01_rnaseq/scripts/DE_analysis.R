library(rtracklayer)
library(tidyverse)
library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggridges)

#Connecting to snakeake inputs and outputs:
files <- unlist(snakemake@input)
out_csv <- snakemake@output[[1]]
out_pdf <- snakemake@output[[2]]
out_session <- snakemake@output[[3]]


#Opening pdf with graphs
pdf(out_pdf, height = 10, width = 15)

### 1. File and reference download ###

## Downloading the reference: hg38 version 49
gtf <- import("references/gencode.v49.annotation.gtf.gz")
#mapping ensembl id and hugo symbols
gtf_genes <- gtf[gtf$type == "gene"]
gene_map <- data.frame(ensembl_id = gtf_genes$gene_id,
                       hugo_symbol = gtf_genes$gene_name)


## Creating the gene-transcript correspondence
#Selecting only "transcripts
gtf_tx <- gtf[gtf$type == "transcript"]
#Creating the table (Using id, not names, because they are unique)
tx2gene <- data.frame(transcript = gtf_tx$transcript_id,   
                      gene = gtf_tx$gene_id)
#Removing duplicates (if any)
tx2gene <- distinct(tx2gene)


## Downloading the count files
#files <- unlist(snakemake@input). Ottenuti all'inizio
names(files) <- paste0("S", 1:12)



#Importing counts
tx_matrix <- tximport(files,
                      type = "salmon",
                      tx2gene = tx2gene, # dataframe: txID, geneID
                      txOut = FALSE, 
                      ignoreAfterBar = TRUE)

## Creating the metadata table
metadata <- data.frame(infection = c("TB", "CTRL", "TB", "CTRL", "TB", "CTRL", "TB", "CTRL", "TB", "CTRL", "TB", "CTRL"),
                       replicate = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6"))
rownames(metadata) <- paste0("S", 1:12)
metadata$infection <- as.factor(metadata$infection)   #CTRL is first
metadata$replicate <- factor(metadata$replicate)




### 2. DE analisys: DESeq2 ###
dds <- DESeqDataSetFromTximport(tx_matrix,
                                colData = metadata,
                                design = ~ replicate + infection)
#Filtering: keeping only genes with at least 10 counts in half of the samples
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

##DESeq2:
dds <- DESeq(dds)




### 3. Explorative analysis ###

## For this experiment there are 6 samples in a paired set-up. For initial visualization only PCA will show clustering and outliers (if any). 
## t-SNE and UMAP need more sample to start being informative so they will not be implemented. 
## Heatmap will confirm clustering and add information on expression. 
## Then boxplots of overall expression for each sample will show whether any samples have a general over/under expression, and will verify that normalization was successful

vst_dds <- vst(dds, blind = FALSE)      #Normalization: variance stabilizing transformation
counts_norm <- assay(vst_dds)           #Extracting counts
counts_norm <- counts_norm[rowVars(counts_norm) > 0, ]  #Deleting rows with variance = 0
counts_norm_t <- t(counts_norm)           #The matrix is needed transposed for pca

##PCA:
pca <- prcomp(counts_norm_t)
summary(pca)

##PCA graph:
#Creating the plot dataframe
pca_df <- as.data.frame(pca$x[, 1:2])
#Adding metadata information
pca_df <- merge(pca_df, metadata, by = "row.names")
pca_df$Row.names <- NULL

#ggplot
ggplot(pca_df, aes(x = PC1, y = PC2, shape = replicate, color = infection)) +
  geom_point(size = 3) +                                                                                    #Adds the data as points
  labs(title = "Explorative analysis - PCA",
       x = paste0("PCA1 (", round(100 * summary(pca)$importance[2,1], 1), "% of the variance)"),    #Title and labels
       y = paste0("PCA2 (", round(100 * summary(pca)$importance[2,2], 1), "% of the variance)")) +
  theme_bw() +                                                                                      #Setting white background
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) +                   #Formatting title
  scale_shape_manual(values= c(15:19, 8)) +                                                            #Setting points shapes. From 15 on shapes are full inside
  scale_color_manual(values = c("CTRL" = "red", "TB" ="darkblue")) +                                #Setting colors
  guides(fill = guide_legend(override.aes = list(shape = 21)))                                      #Adding legend


## Heatmap
#Again using the counts after transposition: counts_norm_t
counts_norm_matrix <- scale(counts_norm_t)
heatmap_matrix <- t(counts_norm_matrix)

Conditions <- data.frame(
  Infection  = metadata$infection,
  Replicate = metadata$replicate
)
row.names(Conditions) <- colnames(heatmap_matrix)

pheat <- pheatmap(heatmap_matrix,
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   annotation_col = Conditions,
                   show_rownames = FALSE,
                   show_colnames = TRUE,
                   main = "Explorative analysis - Heatmap")


##Boxplot for Sample
#For the boxplot using the counts are needed in the "long" format: transforming
counts_norm_long <- counts_norm %>% 
                    as.data.frame() %>% 
                    rownames_to_column("Gene") %>%
                    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "expr")
counts_norm_long$Sample <- factor(counts_norm_long$Sample, levels = paste0("S", 1:12))

#Boxplot
ggplot(counts_norm_long, aes(x = Sample, y = expr)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Explorative analysis - Expression boxplots",
       x = "Sample",
       y = "VST gene expression")

## Results
## PCA -> Strong Clustering for condition (infection). There is also a slight clustering for sample (on PCA2): paired design is correct.
## Heatmap -> Confirms clustering seen for PCA
## boxplot -> From the boxplot profiles it appears that: 1) Normalization was correctly done. 2) Overall expression is coherent among different samples: no sample is an anomaly. 
## 3) All samples have an expected distribution of gene expression, (with a peak for housekeeping genes)  



### 4. Significant genes extraction and visualization ###

##Extracting Differentially expressed genes (by adj_pvalue)
#Considering the contrast on the condition: infection
DE_res <- as.data.frame(results(dds, contrast = c("infection", "TB", "CTRL")))
#Mapping the gene names as HUGO symbols
DE_res$ENSEMBL <- row.names(DE_res)
DE_res <- merge(DE_res, gene_map,
                by.x = "ENSEMBL",
                by.y = "ensembl_id",
                all.x = TRUE)

#Removing duplicates (there are some, as shown by "sum(duplicated(sign_DE_res$hugo_symbol))")
DE_res <- DE_res %>%
  group_by(hugo_symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()

sign_DE_res <- DE_res %>% filter(padj < 0.0001 & abs(log2FoldChange) > 1)
#5265 DE genes for very stringent conditions: there is very strong variability
write.csv(sign_DE_res, out_csv, row.names = FALSE)


##Visualization: volcano plot
#Applying shrinkage:
res_shrunk <- lfcShrink(dds, coef = "infection_TB_vs_CTRL", type = "apeglm")
res_shrunk_df <- as.data.frame(res_shrunk)

ggplot(res_shrunk_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.0001), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(data = res_shrunk_df[!is.na(res_shrunk_df$padj) & 
                                    res_shrunk_df$padj < 0.0001 & 
                                    abs(res_shrunk_df$log2FoldChange) > 1, ], color = "red") +
  theme_minimal() +
  labs(title = "Volcano plot: infection", x = "log2 Fold Change", y = "-log10(adj pvalue)")




### 5. Functional enrichment analysis ###

##DE analysis resulted in a large number of extremely significant genes, and the number of sample of this experiment isn't very high. 
##Therefore, Over Representation Analysis (ORA) is expected to be the most coherent choice for FE analysis, and will be the primary tool in this case.
##A GSEA will also be performed only as a confirmatory tool, and to add information on weaker phenomena. It must therefore be considered to be less significant.


##ORA (ClusterProfiler)
#Defining the background
background <- DE_res$hugo_symbol
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
mapped_symbols <- intersect(background, all_symbols)
length(mapped_symbols)
#Only 17809/23637. So 5794 genes missing genes.

#Now trying to retrieve the missing genes. A possibility is that these are codes for novel/ncRNAs! 
#In that case we don't need to retrieve them, cause they simply don't have associated GO terms. Let's check:
not_in_symbol <- DE_res %>% filter(!hugo_symbol %in% mapped_symbols)   
novel <-grepl("^ENSG", not_in_symbol$hugo_symbol)
sum(novel)
#5738/5794 have an "ENSG" hugo symbol. This generally means these are novel/unmapped genes. Therefore they can be discarded. 
#The other 56 genes must be further analyzed
missing <- not_in_symbol[!novel ,]
print(missing$hugo_symbol)
#These are all ncRNAs. Therefore they can be discarded too.
#Now checking how many of the significantly DE genes remain after discarding novel genes and ncRNAs
sign_GO_genes <- intersect(sign_DE_res$hugo_symbol, all_symbols) #4499. This will be the gene set for ORA
background <- mapped_symbols                                     # This will be the background for ORA
##So the FE analysis will work on 4499/5265 of the significant genes, and 17809/23637 of the background genes.

#Dividing for upregulated and downregulated.
up_sign_DE_res <- sign_DE_res %>% filter(log2FoldChange > 0 & hugo_symbol %in% sign_GO_genes)     #2285 genes
down_sign_DE_res <- sign_DE_res %>% filter(log2FoldChange < 0 & hugo_symbol %in% sign_GO_genes)   #2214 genes

#ORA all
ORA_all_BP <-  enrichGO(gene = sign_GO_genes,       
                        universe = background,        
                        OrgDb = org.Hs.eg.db,       
                        keyType = "SYMBOL",             
                        ont = "BP",                 
                        #pAdjustMethod = "none",         
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)


#ORA up
ORA_up_BP <-  enrichGO(gene = up_sign_DE_res$hugo_symbol,       
                       universe = background,        
                       OrgDb = org.Hs.eg.db,       
                       keyType = "SYMBOL",             
                       ont = "BP",                 
                       #pAdjustMethod = "none",         
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)


#ORA down
ORA_down_BP <-  enrichGO(gene = down_sign_DE_res$hugo_symbol,       
                         universe = background,        
                         OrgDb = org.Hs.eg.db,       
                         keyType = "SYMBOL",             
                         ont = "BP",                 
                         #pAdjustMethod = "none",         
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)


barplot(ORA_all_BP, showCategory = 10, title = "Biological Process - ORA: all genes")
dotplot(ORA_all_BP, showCategory = 10, title = "Biological Process - ORA: all genes")

barplot(ORA_up_BP, showCategory = 10, title = "Biological Process - ORA: upregulated genes")
dotplot(ORA_up_BP, showCategory = 10, title = "Biological Process - ORA: upregulated genes")

barplot(ORA_down_BP, showCategory = 10, title = "Biological Process - ORA: downregulated genes")
dotplot(ORA_down_BP, showCategory = 10, title = "Biological Process - ORA: downregulated genes")

#Results:
length(ORA_all_BP@result$ID)           #5629 terms

# Finding common and unique terms between up and downregulated
up_terms   <- ORA_up_BP@result %>% filter(p.adjust < 0.05) %>% pull(ID)
down_terms <- ORA_down_BP@result %>% filter(p.adjust < 0.05) %>% pull(ID)
#common terms
common_terms <- intersect(up_terms, down_terms)  #5100 terms
#unique terms
unique_up <- setdiff(up_terms, down_terms)    #315 terms
unique_down <- setdiff(down_terms, up_terms)  #214 terms

#Visualizing only unique terms:
ORA_up_unique <- ORA_up_BP
ORA_up_unique@result <- ORA_up_BP@result %>% filter(ID %in% unique_up)
ORA_down_unique <- ORA_down_BP
ORA_down_unique@result <- ORA_down_BP@result %>% filter(ID %in% unique_down)
barplot(ORA_up_unique, showCategory = 10, title = "Biological Process - ORA: unique upregulated genes")
barplot(ORA_down_unique, showCategory = 10, title = "Biological Process - ORA: unique downregulated genes")


##GSEA
#Using the "stat" list from DESeq2 as metric for the ranking
DE_res_GSEA <- DE_res %>% filter(DE_res$hugo_symbol %in% background)
stat_list <- DE_res_GSEA$stat
names(stat_list) <- DE_res_GSEA$hugo_symbol
stat_list <- sort(stat_list, decreasing = TRUE)

set.seed(10)
gsea_res <- gseGO(
  geneList = stat_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  verbose = FALSE
)

#Plot
dotplot(gsea_res, showCategory = 15, title = "Biological Process - GSEA")
p <- ridgeplot(gsea_res, showCategory = 15)         #library(ggridges)
p + ggtitle("Biological Process - GSEA")


dev.off()

### Obtaining session info:
sink(out_session)
sessionInfo()
sink()