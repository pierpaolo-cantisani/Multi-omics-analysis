library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggridges)
#library(rGREAT)

setwd("C:/Users/pierp/Desktop/THESIS PROJECT/Integration analysis")
#pdf("C:/Users/pierp/Desktop/THESIS PROJECT/Integration analysis/Integration analysis (GO) Graphs.pdf")



### 0. Importing and prefiltering data

#Importing RNA-Seq DE genes
sign_DE_res <- read.csv("DE_results.csv")
sign_DE_res <- sign_DE_res %>% dplyr::rename(SYMBOL = hugo_symbol)  #renaming for coherence
#Deleting duplicate genes (keeping the one with the higher log2FC)
sign_DE_res <- sign_DE_res %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup()
#Importing BS-Seq DM sites
DM_sites <- read.csv("DM_sites.csv")

## Obtaining the universe N: 
#Importing RNA-seq universe (all genes considered for the DESeq2 analysis)
DE_res <- read.csv("RNAseq_universe.csv", col.names = c("SYMBOL", "log2FoldChange", "stat"))
#Importing WGBS universe (all genes that had methylated sites with coverage > 5)
WGBS_universe <- read.csv("WGBS_universe.csv", col.names = "SYMBOL")

# Universe N: all genes considered in both differential analyses
universe <- intersect(na.omit(DE_res$SYMBOL), 
                      na.omit(WGBS_universe$SYMBOL))


# Final genes ready for comparison are:
DE_genes <- unique(DE_res$SYMBOL)  # 4144 DE genes
DM_genes <- unique(DM_sites$SYMBOL)    # 277 genes with DM sites

# All genes annotated in the database:
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")





### 1. RNA-Seq. Functional enrichment analysis ###

##DE analysis resulted in a large number of extremely significant genes, and the number of sample of this experiment isn't very high. 
##Therefore, Over Representation Analysis (ORA) is expected to be the most coherent choice for FE analysis, and will be the primary tool in this case.
##A GSEA will also be performed only as a confirmatory tool, and to add information on weaker phenomena. It must therefore be considered to be less significant.


##ORA (ClusterProfiler)
#Defining the background
background <- DE_res$SYMBOL
mapped_symbols <- intersect(background, all_symbols)

#Now trying to retrieve the missing genes. A possibility is that these are codes for novel/ncRNAs! 
#In that case we don't need to retrieve them, cause they simply don't have associated GO terms. Let's check:
not_in_symbol <- DE_res %>% filter(!SYMBOL %in% mapped_symbols)   
novel <-grepl("^ENSG", not_in_symbol$SYMBOL)
sum(novel)
#The ones that have an "ENSG" hugo symbol are novel/unmapped genes. Therefore they can be discarded. 
#The other genes must be further analyzed
missing <- not_in_symbol[!novel ,]
print(missing$SYMBOL)
#If these are all ncRNAs, they can be discarded too.

#Now checking how many of the significantly DE genes remain after discarding novel genes and ncRNAs
sign_GO_genes <- intersect(sign_DE_res$SYMBOL, all_symbols)      # This will be the gene set for ORA
background <- mapped_symbols                                     # This will be the background for ORA

#ORA all
ORA_all_BP <-  enrichGO(gene = sign_GO_genes,       
                        universe = background,        
                        OrgDb = org.Hs.eg.db,       
                        keyType = "SYMBOL",             
                        ont = "BP",                 
                        #pAdjustMethod = "none",         
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

ORA_all_BP_simp <- simplify(ORA_all_BP, 
                            cutoff = 0.7, 
                            by = "p.adjust")


#Visualizations:
bar_all_RNA <- barplot(ORA_all_BP_simp, showCategory = 10, title = "RNA-Seq BP - ORA: all genes")
dot_all_RNA <- dotplot(ORA_all_BP_simp, showCategory = 10, title = "RNA-Seq BP - ORA: all genes")
plot_grid(bar_all_RNA, dot_all_RNA, nrow = 1, ncol = 2)

#Results:
DE_terms  <- ORA_all_BP_simp@result %>% filter(p.adjust < 0.05) %>% dplyr::pull(ID)





### 2. WGBS: Functional enrichment analysis ###

## All DM
background_BS <- unique(WGBS_universe$SYMBOL)
mapped_symbols_BS <- intersect(background_BS, all_symbols)
#Final lists:
DM_GO_genes <- intersect(DM_sites$SYMBOL, all_symbols)           # This will be the gene set for ORA
background_BS <- mapped_symbols_BS                                     # This will be the background for ORA

ORA_BS_all_BP <-  enrichGO(gene = DM_GO_genes,       
                        universe = background_BS,        
                        OrgDb = org.Hs.eg.db,       
                        keyType = "SYMBOL",             
                        ont = "BP",                 
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)

ORA_BS_all_BP_simp <- simplify(ORA_BS_all_BP, 
                                  cutoff = 0.7, 
                                  by = "p.adjust")

bar_all_BS <- barplot(ORA_BS_all_BP_simp, showCategory = 10, title = "WGBS BP - ORA: all genes")
dot_all_BS <- dotplot(ORA_BS_all_BP_simp, showCategory = 10, title = "WGBS BP - ORA: all genes")
plot_grid(bar_all_BS, dot_all_BS, nrow = 1, ncol = 2)

#Results:
DM_terms  <- ORA_BS_all_BP_simp@result %>% filter(p.adjust < 0.05) %>% pull(ID)





### 3. Integration: DE terms vs DM terms ###

common_DM_DE <- intersect(DM_terms, DE_terms)

## Graph 1: Euler diagram - Significant GO terms intersection ##

#All
fit <- euler(c(
  "DE terms"           = length(DE_terms) - length(common_DM_DE),
  "DM terms"           = length(DM_terms) - length(common_DM_DE),
  "DE terms&DM terms"        = length(common_DM_DE)
))
plot(fit,
     quantities = list(cex = 1.1),
     labels = list(cex = 1.1),
     fills = list(fill = c("#4393C3", "#D6604D", "#92C592"), alpha = 0.6),
     edges = TRUE,
     main = paste0("Overlap: DE vs DM terms")
)



## Graph 2: Dotplot of top ranked shared terms
shared_results <- bind_rows(
  ORA_all_BP_simp@result   %>% filter(ID %in% common_DM_DE) %>% mutate(omic = "DE"),
  ORA_BS_all_BP_simp@result %>% filter(ID %in% common_DM_DE) %>% mutate(omic = "DM")
)

# Cleaning GeneRatio (From string "k/n", to [0,1])
shared_results <- shared_results %>%
  mutate(GeneRatio_num = sapply(GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  }))
# Ranking genes according to -log10(p.adjust). So genes highly significantly in BOTH are ranked higher
top_terms <- shared_results %>%
  group_by(Description) %>%
  summarise(combined_score = sum(-log10(p.adjust))) %>%
  slice_max(combined_score, n = 20) %>%
  pull(Description)

shared_results_top <- shared_results %>%
  filter(Description %in% top_terms)

# Plot
ggplot(shared_results_top, 
       aes(x = omic, y = Description, size = GeneRatio_num, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#D6604D", high = "#4393C3",
                       name = "p.adjust") +
  scale_size_continuous(name = "Gene Ratio", range = c(2, 8)) +
  labs(
    title    = "Shared GO terms: DE vs DM",
    subtitle = paste0(length(common_DM_DE), " terms in common"),
    x        = NULL,
    y        = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )



#dev.off()








### +: Intergenic regions analysis

gr_intergenic <- as(peakAnno_df, "GRanges")[peakAnno_df$annotation == "Distal Intergenic"]
job <- submitGreatJob(gr_intergenic, species = "hg38")
# Results table
res <- getEnrichmentTables(job)

# Associated genes:
inter_genes_df <- as.data.frame(getRegionGeneAssociations(job))
# Associated GO terms (filtered by significance)
go_bp <- res[["GO Biological Process"]]
go_bp_sig <- go_bp[go_bp$Binom_Adjp_BH < 0.05, ]
go_bp_sig <- go_bp_sig[order(go_bp_sig$Binom_Adjp_BH), ]  # ordina per p-value