library(dplyr)
library(ggplot2)
library(eulerr)
library(tidyr)
library(cowplot)


setwd("C:/Users/pierp/Desktop/THESIS PROJECT/Integration analysis")
pdf("C:/Users/pierp/Desktop/THESIS PROJECT/Integration analysis/Integration analysis Graphs.pdf")


### !!! Summary:
#This analysis answer the following questions (p-value < 0.05 means YES):

## Section 1
# (Significance of intersection): Are the DM genes also differentially expressed (significantly)?
## Section 2
# All DM (p-value): Do DM genes have a trend towards up/down regulation?
# DM ∩ DE DM (p-value): Do DM ∩ DE genes have a trend towards up/down regulation?
# DM ∩ DE vs All DM (p-value): Do DM ∩ DE have a significantly higher trend towards up/down regulation compared to All DM?
# Hypomethylated DM: Do hypomethylatedDM genes have a trend towards up/down regulation?
# Hypermethylated DM: Do hypermethylated DM genes have a trend towards up/down regulation?
## Section 3
# (Trend towards region type): Are the genes with DM sites in this region more differentially expressed compared to other regions?
## Section 4
# (meth vs expr): Is the magnitude of the methylation (in the category) correlated to the genes up/down regulation?
## Section 5
# (multi DM genes) : Are genes with 2 or more (non intergenic) methylation sites more likely to be DE than those with 1?
## Section 6
# (sites vs log2FC): Do the number of DM sites in a gene correlate with its expression (making it more likely DE)?
# (sites vs |log2FC|): Do the number of DM sites in a DE gene correlate with the magnitude of its differential expression?

## -- final_table -- contains results for Section 1, 3, 4, 5, 6
## -- final_expr_table -- contains results for Section 2




### 0. Importing and prefiltering data

#Importing RNA-Seq DE genes
DE_results <- read.csv("DE_results.csv")
DE_results <- DE_results %>% rename(SYMBOL = hugo_symbol)  #renaming for coherence
#Importing BS-Seq DM sites
DM_sites <- read.csv("DM_sites.csv")

## Obtaining the universe N: 
#Importing RNA-seq universe (all genes considered for the DESeq2 analysis)
RNAseq_universe <- read.csv("RNAseq_universe.csv", col.names = c("SYMBOL", "log2FoldChange", "padj"))
#Importing WGBS universe (all genes that had methylated sites with coverage > 5)
WGBS_universe <- read.csv("WGBS_universe.csv")
# Universe N: all genes considered in both differential analyses
universe <- intersect(na.omit(RNAseq_universe$SYMBOL), 
                      na.omit(WGBS_universe$SYMBOL))


## Filtering the DE_results and DE_sites for the Universe. 
# This is needed because only genes that were present in both analyses can be compared.
DE_results <- DE_results %>% filter(SYMBOL %in% universe)
DM_sites <- DM_sites  %>% filter(SYMBOL %in% universe)

# Final genes ready for comparison are:
DE_genes <- unique(DE_results$SYMBOL)
DM_genes <- unique(DM_sites$SYMBOL)

## Dividing DM_sites in 3 groups:
#Group A: DMR on promoters only: the most relevant, but few. Additional specific analysis
#Group B: DMR on promoters/exons/introns/UTR: All the relevant DMR. Will be used for the statistical analyses, since numerous.
#Group C: DMR on intergenic regions: excluded from the intersection analyses. Will be investigated separately with GREAT

#Group A:
prom_DM <- DM_sites[(DM_sites$annotation == 'Promoter (<=1kb)') | (DM_sites$annotation == 'Promoter (1-2kb)'), ]
#Group B:
rel_DM <- DM_sites[DM_sites$annotation != 'Distal Intergenic', ]
#Group C:
interg_DM <- DM_sites[DM_sites$annotation == 'Distal Intergenic', ]





### 1. Analysis of the intersection between DE genes and DM sites ###

rel_intersect_genes <- unique(rel_DM$SYMBOL[rel_DM$SYMBOL %in% DE_results$SYMBOL])
all_intersect_df <- merge(DM_sites[, c("seqnames", "pvalue", "qvalue", "meth.diff", "annotation", "SYMBOL")], DE_results[, c("log2FoldChange", "padj", "SYMBOL")], by = "SYMBOL")
rel_intersect_df <- all_intersect_df %>% filter(SYMBOL %in% rel_intersect_genes)
rel_genes <- unique(rel_DM$SYMBOL)

# Getting the percentage of intersection for Group B
rel_percent <- length(rel_intersect_genes)*100/length(universe)


## Hypergeometric test
#How statistically significant is this intersection result? Is the intersection obtained by chance?

N <- as.numeric(length(universe))                             # universe: all genes considered for both differential analyses
m <- as.numeric(length(DE_genes))                             # all DE_genes
k <- as.numeric(length(rel_genes))                            # all DM genes non-intergenic
q <- as.numeric(length(rel_intersect_genes))                  # intersecting genes
## Ho un urna con N palline. Ce ne sono m rosse e n=N-m bianche. 
## Qual è la probabilità che pescandone k ne ottengo q rosse?

# Formula is:
# P(q = observed success, m = success in the universe, n = insuccess in the universe, k = # of extractions,)
# lower.tail converts P(X <= q) into P(X > q). "q-1" is needed because it's > and not >=.
hyp_intersect <- phyper(q-1, m, N-m, k, lower.tail = FALSE)


#Output on table: # of intersection, percentages
n_rel_DM <- as.numeric(length(unique(rel_DM$SYMBOL)))
n_rel_int_DM <- as.numeric(length(rel_intersect_genes))
rel_DM_percent <- as.numeric(n_rel_int_DM*100/n_rel_DM)
rel_percent <- length(rel_intersect_genes)*100/length(universe)
final_table <- data.frame("All relevant" = c(n_rel_DM, n_rel_int_DM, paste0(round(rel_DM_percent, 1), "%"), paste0(rel_percent, "%"), hyp_intersect, "/"))
row.names(final_table) <- c("DM genes", "DM ∩ DE genes", "fraction of intersected", "fraction of universe", "p-value (Significance of intersection)", "p-value (Trend towards region type)")




## Graph 1: Euler diagram: DE / DM / Intersection ##

fit <- euler(c(
  "DE"           = length(DE_genes) - length(rel_intersect_genes),
  "DM"           = n_rel_DM - length(rel_intersect_genes),
  "DE&DM"        = length(rel_intersect_genes)
))

plot(fit,
     quantities = list(cex = 1.1),
     labels = list(cex = 1.1),
     fills = list(fill = c("#4393C3", "#D6604D", "#92C592"), alpha = 0.6),
     edges = TRUE,
     main = paste0("Overlap: DE and DM genes (universe = ", length(universe), " genes)")
)





### Exploring trends ###


### 2: Do intersected genes have a trend towards up/down regulation? ###

#2.1: UPREGULATION
DE_up <- DE_results %>% filter(log2FoldChange > 0)

#Obtaining up and down intersecting genes
rel_intersect_up <- rel_intersect_df %>% filter(log2FoldChange > 0)
rel_intersect_up_genes <- unique(rel_intersect_up$SYMBOL)     # upregolated intersecting genes
rel_intersect_down <- rel_intersect_df %>% filter(log2FoldChange < 0)
rel_intersect_down_genes <- unique(rel_intersect_down$SYMBOL) # downregolated intersecting genes

# universe restricted to DE_genes only, because we are testing 
# directionality of expression within the DE subset, not across all genes
N <- as.numeric(length(DE_genes))                             # universe: all DE_genes
m <- as.numeric(length(DE_up$SYMBOL))                         # all DE up genes
k <- as.numeric(length(rel_intersect_genes))                  # all intersected DM genes
q <- as.numeric(length(rel_intersect_up_genes))               # upregulated intersected DM genes
## Ho un urna con N palline. Ce ne sono m rosse e n=N-m bianche. 
## Qual è la probabilità che pescandone k ne ottengo q rosse?

hyp_up <- phyper(q-1, m, N-m, k, lower.tail = FALSE)


## But what if all DM genes are significantly enriched for upregulation? If that is the case this result is less interesting. 

##Let's see:
#Obtaining all up and down non_DM genes (also non DE)
non_DM_genes <- setdiff(universe, DM_genes)
non_DM_up <- RNAseq_universe %>% filter(log2FoldChange > 0 & SYMBOL %in% non_DM_genes)
non_DM_down <- RNAseq_universe %>% filter(log2FoldChange < 0 & SYMBOL %in% non_DM_genes)
#Obtaining all up genes (also non DE)
all_up <- RNAseq_universe %>% filter(log2FoldChange > 0 & SYMBOL %in% universe)
DM_up_genes <- all_up$SYMBOL[all_up$SYMBOL %in% DM_genes]

## Hypergeometric test on all rel DM sites
#How significant is the trend of all DM genes towards being upregulated
N <- as.numeric(length(universe))                             # universe: all genes used for analyses
m <- as.numeric(length(all_up$SYMBOL))                        # all upregulated genes (also non DE)
k <- as.numeric(length(DM_genes))                             # all DM genes
q <- as.numeric(length(DM_up_genes))                          # upregulated DM genes
## Ho un urna con N palline. Ce ne sono m rosse e n=N-m bianche. 
## Qual è la probabilità che pescandone k ne ottengo q rosse?

# Formula is:
# P(q = observed success, m = success in the universe, n = insuccess in the universe, k = # of extractions,)
# lower.tail converts P(X <= q) into P(X > q). "q-1" is needed because it's > and not >=.
hyp_DM_up <- phyper(q-1, m, N-m, k, lower.tail = FALSE)


## The next question is: is the DM ∩ DE case significantly more upregulated than all DM?

# Up/down for relevant DM genes (Group B, non-intergenic)
rel_DM_up_genes   <- all_up$SYMBOL[all_up$SYMBOL %in% rel_genes]

## Hypergeometric test on all rel DM sites
#How significant is the trend of all DM genes towards being upregulated
N <- as.numeric(length(rel_genes))                            # all rel DM genes
m <- as.numeric(length(unique(rel_DM_up_genes)))              # all DM relevant upregulated genes (also non DE)
k <- as.numeric(length(rel_intersect_genes))                  # all intersected DM genes
q <- as.numeric(length(rel_intersect_up_genes))               # upregulated intersected DM genes
## Ho un urna con N palline. Ce ne sono m rosse e n=N-m bianche. 
## Qual è la probabilità che pescandone k ne ottengo q rosse?

# Formula is:
# P(q = observed success, m = success in the universe, n = insuccess in the universe, k = # of extractions,)
# lower.tail converts P(X <= q) into P(X > q). "q-1" is needed because it's > and not >=.
hyp_fin_up <- phyper(q-1, m, N-m, k, lower.tail = FALSE)



##2.2: DOWNREGULATION:

DE_down <- DE_results %>% filter(log2FoldChange < 0)

# All down genes in universe (also non DE)
all_down <- RNAseq_universe %>% filter(log2FoldChange < 0 & SYMBOL %in% universe)
DM_down_genes <- all_down$SYMBOL[all_down$SYMBOL %in% DM_genes]
rel_DM_down_genes <- all_down$SYMBOL[all_down$SYMBOL %in% rel_genes]

# Hypergeometric test: are all DM genes enriched for downregulation?
N <- as.numeric(length(universe))                 # universe: all genes common to the two analyses
m <- as.numeric(length(all_down$SYMBOL))          # all downregulated genes in universe
k <- as.numeric(length(DM_genes))                 # all DM genes
q <- as.numeric(length(DM_down_genes))            # downregulated DM genes
hyp_DM_down <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

# Hypergeometric test: are DM ∩ DE genes enriched for downregulation?
N <- as.numeric(length(DE_genes))                 # universe: all DE genes
m <- as.numeric(length(DE_down$SYMBOL))           # all DE down genes
k <- as.numeric(length(rel_intersect_genes))      # all intersected DM genes
q <- as.numeric(length(rel_intersect_down_genes)) # downregulated intersected DM genes
hyp_down <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

# Hypergeometric test: is DM ∩ DE more downregulated than all relevant DM?
N <- as.numeric(length(rel_genes))                # all relevant (on gene body/promoter) DM genes
m <- as.numeric(length(rel_DM_down_genes))        # all DM relevant downregulated genes (also non DE)
k <- as.numeric(length(rel_intersect_genes))      # all relevant intersecting genes
q <- as.numeric(length(rel_intersect_down_genes)) # all downregulated intersected genes
hyp_fin_down <- phyper(q-1, m, N-m, k, lower.tail = FALSE)



#Output on table
final_expr_table <- data.frame("Upregulated" = c(hyp_DM_up, hyp_up, hyp_fin_up), "Downregulated" = c(hyp_DM_down, hyp_down, hyp_fin_down))
row.names(final_expr_table) <- c("p-value (All DM)", "p-value (DM ∩ DE)", "p-value (DM ∩ DE vs All DM)")


## Graph 2: Up/Down regulation across DE, DM, and DE∩DM genes ##

upreg_summary <- data.frame(
  group = c("All DE genes", "All DM genes", "DM (only gene body/promoter)", "DM ∩ DE genes"),
  up   = c(length(DE_up$SYMBOL),
           length(DM_up_genes),
           length(rel_DM_up_genes),
           length(rel_intersect_up_genes)),
  down = c(length(DE_genes) - length(DE_up$SYMBOL),
           length(DM_down_genes),
           length(rel_DM_down_genes),
           length(rel_intersect_down_genes))
)

upreg_long <- upreg_summary %>%
  pivot_longer(cols = c(up, down),
               names_to  = "direction",
               values_to = "n") %>%
  group_by(group) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100,
    group = factor(group, levels = c("All DE genes", "All DM genes", "DM (only gene body/promoter)", "DM ∩ DE genes")),
    direction = factor(direction, levels = c("up", "down"),
                       labels = c("Upregulated", "Downregulated"))
  )

totals_upreg <- upreg_long %>%
  distinct(group, total)

p_upreg <- ggplot(upreg_long, aes(x = group, y = pct, fill = direction)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(
    data = totals_upreg,
    aes(x = group, y = 103, label = paste0("n = ", total)),
    inherit.aes = FALSE,
    size = 4, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Upregulated" = "#D6604D", "Downregulated" = "#4393C3"),
    name   = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 108),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    title    = "Expression direction across gene groups",
    subtitle = "Percentage of up- vs downregulated genes",
    x        = NULL,
    y        = "Percentage of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    plot.title       = element_text(face = "bold"),
    axis.text.x      = element_text(size = 11)
  )
print(p_upreg)



##2.3: Association: UP/DOWN DE vs Hyper/Hypo DM

#Obtaining hypo/hyper lists
hypo_DM <- rel_DM %>% filter(meth.diff < 0)
hyper_DM <- rel_DM %>% filter(meth.diff > 0)
hypo_int_DM <- hypo_DM %>% filter(SYMBOL %in% DE_results$SYMBOL)
hyper_int_DM <- hyper_DM %>% filter(SYMBOL %in% DE_results$SYMBOL)
#I'm not filtering for unique genes. So there will be duplicates: this analysis will be done for site, rather than gene
hypo_int_up_DM <- hypo_DM %>% filter(SYMBOL %in% DE_up$SYMBOL)
hyper_int_up_DM <- hyper_DM %>% filter(SYMBOL %in% DE_up$SYMBOL)
hypo_int_down_DM <- hypo_DM %>% filter(SYMBOL %in% DE_down$SYMBOL)
hyper_int_down_DM <- hyper_DM %>% filter(SYMBOL %in% DE_down$SYMBOL)

# Hypergeometric test: Are hypo/hyper DM sites associated with strong up/downregulation?
# hypo-up
N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))        # All intersecting genes
m <- as.numeric(length(hypo_int_up_DM$SYMBOL) + length(hyper_int_up_DM$SYMBOL))  # All up intersecting genes
k <- as.numeric(length(hypo_int_DM$SYMBOL))                                      # All hypo intersecting genes
q <- as.numeric(length(hypo_int_up_DM$SYMBOL))                                   # hypo-up intersecting genes
hyp_up_hypo <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
# hypo-down
N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))            # All intersecting genes
m <- as.numeric(length(hypo_int_down_DM$SYMBOL) + length(hyper_int_down_DM$SYMBOL))  # All down intersecting genes
k <- as.numeric(length(hypo_int_DM$SYMBOL))                                          # All hypo intersecting genes
q <- as.numeric(length(hypo_int_down_DM$SYMBOL))                                     # hypo-down intersecting genes
hyp_down_hypo <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
# hyper-up
N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))        # All intersecting genes
m <- as.numeric(length(hypo_int_up_DM$SYMBOL) + length(hyper_int_up_DM$SYMBOL))  # All up intersecting genes
k <- as.numeric(length(hyper_int_DM$SYMBOL))                                     # All hyper intersecting genes
q <- as.numeric(length(hyper_int_up_DM$SYMBOL))                                  # hyper-up intersecting genes
hyp_up_hyper <- phyper(q-1, m, N-m, k, lower.tail = FALSE)
# hyper-down
N <- as.numeric(length(hypo_int_DM$SYMBOL) + length(hyper_int_DM$SYMBOL))            # All intersecting genes
m <- as.numeric(length(hypo_int_down_DM$SYMBOL) + length(hyper_int_down_DM$SYMBOL))  # All down intersecting genes
k <- as.numeric(length(hyper_int_DM$SYMBOL))                                         # All hypo intersecting genes
q <- as.numeric(length(hyper_int_down_DM$SYMBOL))                                    # hypo-down intersecting genes
hyp_down_hyper <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#Output on table:
final_expr_table <- rbind(final_expr_table, "p-value (Hypomethylated DM)" = c(hyp_up_hypo, hyp_down_hypo))
final_expr_table <- rbind(final_expr_table, "p-value (Hypermethylated DM)" = c(hyp_up_hyper, hyp_down_hyper))

##! Duplicates were not cleaned: this is wanted, as different sites on the same gene may have different trends.
##  In this way all information is kept. Therefore the analysis is for "site", not for "gene".



## Graph 3: Up/Down regulation in relation to hypo/hypermethylation ##

upreg_summary <- data.frame(
  group = c("Hypomethylated", "Hypermethylated"),
  up   = c(length(hypo_int_up_DM$SYMBOL),
           length(hyper_int_up_DM$SYMBOL)),
  down = c(length(hypo_int_down_DM$SYMBOL),
           length(hyper_int_down_DM$SYMBOL))
)

upreg_long <- upreg_summary %>%
  pivot_longer(cols = c(up, down),
               names_to  = "direction",
               values_to = "n") %>%
  group_by(group) %>%
  mutate(
    total = sum(n),
    pct   = n / total * 100,
    group = factor(group, levels = c("Hypomethylated", "Hypermethylated")),
    direction = factor(direction, levels = c("up", "down"),
                       labels = c("Upregulated", "Downregulated"))
  )

totals_upreg <- upreg_long %>%
  distinct(group, total)

p_upreg <- ggplot(upreg_long, aes(x = group, y = pct, fill = direction)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(
    data = totals_upreg,
    aes(x = group, y = 103, label = paste0("n = ", total)),
    inherit.aes = FALSE,
    size = 4, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Upregulated" = "#D6604D", "Downregulated" = "#4393C3"),
    name   = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 108),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(
    title    = "Expression direction relation to hypo/hypermethylation sites",
    subtitle = "DM (only gene body/promoter) ∩ DE gene set",
    x        = NULL,
    y        = "Percentage of sites"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    plot.title       = element_text(face = "bold"),
    axis.text.x      = element_text(size = 11)
  )
print(p_upreg)


# Graph 4:Volcano plot meth.diff vs log2FC
all_DM_DE_df <- merge(RNAseq_universe, DM_sites, by = "SYMBOL")
ggplot(all_DM_DE_df, aes(x = meth.diff, y = log2FoldChange)) +
         geom_point(alpha = 0.5) +
         geom_point(data = all_DM_DE_df[all_DM_DE_df$padj < 0.05 & abs(all_DM_DE_df$log2FoldChange) > 1, ], color = "red") +
         theme_minimal() +
         geom_hline(yintercept = 0, color = "blue") +
         geom_vline(xintercept = 0, color = "blue") +
         xlim(-max(abs(all_DM_DE_df$meth.diff), na.rm = TRUE), max(abs(all_DM_DE_df$meth.diff), na.rm = TRUE)) +
         ylim(-max(abs(all_DM_DE_df$log2FoldChange), na.rm = TRUE), max(abs(all_DM_DE_df$log2FoldChange), na.rm = TRUE)) +
         labs(title = "Volcano plot: DM meth.diff vs DE log2FC", x = "meth.diff", y = "log2 Fold Change")





### 3: Are intersected genes specifically associated to exons/introns/promoters/3UTR? ###

# DMR on promoters only

#Intersecting with DE genes
prom_intersect_genes <- unique(prom_DM$SYMBOL[prom_DM$SYMBOL %in% DE_results$SYMBOL])
prom_intersect_df <- all_intersect_df %>% filter(SYMBOL %in% prom_intersect_genes)
# 
n_prom_DM <- as.numeric(length(unique(prom_DM$SYMBOL)))
n_prom_int_DM <- as.numeric(length(prom_intersect_genes))

## These genes must be valued as the most relevant for this study. More considerations on them later

## DMR on exon/intron/UTR genes

exon_genes <- unique(rel_DM$SYMBOL[grepl("^Exon", rel_DM$annotation)])
intron_genes <- unique(rel_DM$SYMBOL[grepl("^Intron", rel_DM$annotation)])
UTR3_genes <- unique(rel_DM$SYMBOL[grepl("^3' UTR", rel_DM$annotation)])
#Obtaining only exon/intron/UTR intersecting genes
exon_inter_genes <- unique(rel_intersect_df$SYMBOL[grepl("^Exon", rel_intersect_df$annotation)])
intron_inter_genes <- unique(rel_intersect_df$SYMBOL[grepl("^Intron", rel_intersect_df$annotation)])
UTR3_inter_genes <- unique(rel_intersect_df$SYMBOL[grepl("^3' UTR", rel_intersect_df$annotation)])

## The exon and 3'UTR cases have numbers very low: must be careful with the test interpretation. 

## Hypergeometric tests
#Are intersecting genes associated to methylation specific to introns/exons/3UTR/promoters?

#exon
N <- as.numeric(length(unique(rel_DM$SYMBOL)))                # all DM genes
m <- as.numeric(length(exon_genes))                           # exon DM genes
k <- as.numeric(length(rel_intersect_genes))                  # all intersecting DM genes
q <- as.numeric(length(exon_inter_genes))                     # intersecting exon DM genes
hyp_ex <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#intron
N <- as.numeric(length(unique(rel_DM$SYMBOL)))                # all DM genes
m <- as.numeric(length(intron_genes))                         # intron DM genes
k <- as.numeric(length(rel_intersect_genes))                  # all intersecting DM genes
q <- as.numeric(length(intron_inter_genes))                   # intersecting intron DM genes
hyp_int <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#UTR3
N <- as.numeric(length(unique(rel_DM$SYMBOL)))                # all DM genes
m <- as.numeric(length(UTR3_genes))                           # 3UTR DM genes
k <- as.numeric(length(rel_intersect_genes))                  # all intersecting DM genes
q <- as.numeric(length(UTR3_inter_genes))                     # intersecting 3UTR on DM genes
hyp_UTR <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#Promoter
N <- as.numeric(length(unique(rel_DM$SYMBOL)))                # all DM genes
m <- as.numeric(n_prom_DM)                                    # prom DM genes
k <- as.numeric(length(rel_intersect_genes))                  # all intersecting DM genes
q <- as.numeric(n_prom_int_DM)                                # intersecting prom on DM genes
hyp_prom_trend <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

## Ho un urna con N palline. Ce ne sono m rosse e n=N-m bianche. 
## Qual è la probabilità che pescandone k ne ottengo q rosse?

perc_prom <- as.numeric(n_prom_int_DM*100/n_prom_DM)
perc_prom_univ <- length(prom_intersect_genes)*100/length(universe)
perc_ex = as.numeric(length(exon_inter_genes))*100/as.numeric(length(exon_genes))
perc_ex_univ = as.numeric(length(exon_inter_genes))*100/as.numeric(length(universe))
perc_int = as.numeric(length(intron_inter_genes))*100/as.numeric(length(intron_genes))
perc_int_univ = as.numeric(length(intron_inter_genes))*100/as.numeric(length(universe))
perc_UTR = as.numeric(length(UTR3_inter_genes))*100/as.numeric(length(UTR3_genes))
perc_UTR_univ = as.numeric(length(UTR3_inter_genes))*100/as.numeric(length(universe))
#Output on table
final_table <- cbind(final_table, Promoter = c(n_prom_DM, n_prom_int_DM, paste0(round(perc_prom, 1), "%"), paste0(perc_prom_univ, "%"), "/", hyp_prom_trend))
final_table <- cbind(final_table, exon = c(as.numeric(length(exon_genes)), as.numeric(length(exon_inter_genes)) , paste0(round(perc_ex, 1), "%"), paste0(perc_ex_univ, "%"), "/", hyp_ex))
final_table <- cbind(final_table, intron = c(as.numeric(length(intron_genes)), as.numeric(length(intron_inter_genes)) , paste0(round(perc_int, 1), "%"), paste0(perc_int_univ, "%"), "/", hyp_int))
final_table <- cbind(final_table, "3' UTR" = c(as.numeric(length(UTR3_genes)), as.numeric(length(UTR3_inter_genes)) , paste0(round(perc_UTR, 1), "%"), paste0(perc_UTR_univ, "%"), "/", hyp_UTR))



## Graph 5: Region type ##
# Building counts
region_summary <- data.frame(
  region = c("Promoter", "Exon", "Intron", "3' UTR"),
  all_DM = c(n_prom_DM, length(exon_genes), length(intron_genes), length(UTR3_genes)),
  inter_DM = c(n_prom_int_DM, length(exon_inter_genes), length(intron_inter_genes), length(UTR3_inter_genes))
)

# Obtaining long format for ggplot
region_long <- region_summary %>%
  pivot_longer(cols = c(all_DM, inter_DM),
               names_to = "group",
               values_to = "n") %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100,
         group = factor(group,
                        levels = c("all_DM", "inter_DM"),
                        labels = c("All DM genes", "DM ∩ DE genes")),
         region = factor(region, levels = c("Promoter", "Exon", "Intron", "3' UTR")))

region_colors <- c(
  "Promoter" = "#E41A1C",
  "Exon"     = "#FF7F00",
  "Intron"   = "#4393C3",
  "3' UTR"   = "#984EA3"
)

# Tot of annotation on bar
totals <- region_long %>%
  group_by(group) %>%
  summarise(total = sum(n), .groups = "drop")

p_region <- ggplot(region_long, aes(x = group, y = pct, fill = region)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(
    data = totals,
    aes(x = group, y = 102, label = paste0("n = ", total)),
    inherit.aes = FALSE,
    size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = region_colors, name = "Genomic region") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "Genomic region distribution",
    subtitle = "All DM genes vs DM ∩ DE genes",
    x = NULL,
    y = "Percentage of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 12)
  )

print(p_region)





### 4: The more the sites are differentially methylated (in magnitude), the more the gene is up/down regulated? (Correlation between magnitude of methyl and up/down regulation) ###

#Considering unique genes: meth will be the mean of the sites
gene_level_df <- all_intersect_df %>%
  group_by(SYMBOL) %>%
  summarise(
    mean_meth = mean(meth.diff),
    log2FoldChange = unique(log2FoldChange))

prom_intersect_df <- prom_intersect_df %>%
  group_by(SYMBOL) %>%
  summarise(
    mean_meth = mean(meth.diff),
    log2FoldChange = unique(log2FoldChange))

#Doing it for intersected sites, and then for all sites
exon_intersect_df <- gene_level_df %>% filter(SYMBOL %in% exon_inter_genes)
intron_intersect_df <- gene_level_df %>% filter(SYMBOL %in% intron_inter_genes)

#Intersected
#all
all_cor <- cor.test(gene_level_df$mean_meth, gene_level_df$log2FoldChange, method = "spearman")
#prom
prom_cor <- cor.test(prom_intersect_df$mean_meth, prom_intersect_df$log2FoldChange, method = "spearman")
#exon
exon_cor <- cor.test(exon_intersect_df$mean_meth, exon_intersect_df$log2FoldChange, method = "spearman")
#intron
intron_cor <- cor.test(intron_intersect_df$mean_meth, intron_intersect_df$log2FoldChange, method = "spearman")

#Adding output to table:
final_table <- rbind(final_table, "Rho corr (meth vs expr)" = c(all_cor$estimate, prom_cor$estimate, exon_cor$estimate, intron_cor$estimate, "/"))
final_table <- rbind(final_table, "p-value corr (meth vs expr)" = c(all_cor$p.value, prom_cor$p.value, exon_cor$p.value, intron_cor$p.value, "/"))

## Graphs 6: Scatter plot ##
# Scatter plot
all_scatt <- ggplot(gene_level_df, aes(x = mean_meth, y = log2FoldChange)) +
                geom_point() +
                geom_smooth(method = "lm") +
                labs(title = "All Relevant methylation vs Expression")
prom_scatt <- ggplot(prom_intersect_df, aes(x = mean_meth, y = log2FoldChange)) +
                geom_point() +
                geom_smooth(method = "lm") +
                labs(title = "Promoter methylation vs Expression")
exon_scatt <- ggplot(exon_intersect_df, aes(x = mean_meth, y = log2FoldChange)) +
                geom_point() +
                geom_smooth(method = "lm") +
                labs(title = "Exon methylation vs Expression")
intron_scatt <- ggplot(intron_intersect_df, aes(x = mean_meth, y = log2FoldChange)) +
                geom_point() +
                geom_smooth(method = "lm") +
                labs(title = "Intron methylation vs Expression")
plot_grid(all_scatt, prom_scatt, exon_scatt, intron_scatt, nrow = 2, ncol = 2)





### 5: Effects of number of DM sites on Expression ###

##5.1: The more DM sites a gene has, the more probability of being DE it has?

#Extracting the genes with more than 1 DM site
#Choosing to exclude intergenic DM from this. So the universe of this analysis is only "Relevant genes". 
DM_num <- as.data.frame(table(rel_DM$SYMBOL))
colnames(DM_num) <- c("gene", "frequency")
multi_DM <- rel_DM %>% filter(SYMBOL %in% DM_num$gene[DM_num$frequency >= 2])
multi_DM_genes <- unique(multi_DM$SYMBOL)
#Intersecting:
multi_rel_intersect_genes <- intersect(multi_DM_genes, rel_intersect_genes)

## Hypergeometric test: Are genes with 2 or more (non intergenic) methyl sites more likely to be DE than those with 1?
N <- as.numeric(length(rel_genes))                              # all rel DM genes
m <- as.numeric(length(rel_intersect_genes))                    # all intersecting DM/DE genes
k <- as.numeric(length(multi_DM_genes))                         # genes with 2 or more DM sites
q <- as.numeric(length(multi_rel_intersect_genes))              # intersecting genes with 2 or more DM sites
hyp_multi <- phyper(q-1, m, N-m, k, lower.tail = FALSE)

#Adding to table:
final_table <- rbind(final_table, "p-value (multi DM genes)" = c(hyp_multi, "/", "/", "/", "/"))



## Graph 7: Probability of being DE depending on number of DM sites ##
# Summary dataframe
single_DM <- rel_DM %>% filter(SYMBOL %in% DM_num$gene[DM_num$frequency < 2])
single_DM_genes <- unique(single_DM$SYMBOL)
single_rel_intersect_genes <- setdiff(rel_intersect_genes, multi_rel_intersect_genes)
dm_summary <- data.frame(
  group    = c("All DM genes", "DM ∩ DE genes"),
  single   = c(length(single_DM_genes),  length(single_rel_intersect_genes)),
  multi    = c(length(multi_DM_genes),   length(multi_rel_intersect_genes))
)

# Long format
dm_long <- dm_summary %>%
  pivot_longer(cols = c(single, multi),
               names_to  = "dm_class",
               values_to = "n") %>%
  group_by(group) %>%
  mutate(
    total    = sum(n),
    pct      = n / total * 100,
    group    = factor(group, levels = c("All DM genes", "DM ∩ DE genes")),
    dm_class = factor(dm_class,
                      levels = c("single", "multi"),
                      labels = c("1 DM site", "≥2 DM sites"))
  )

# Totals for bar annotation
totals_dm <- dm_long %>%
  distinct(group, total)

# Colors
dm_colors <- c(
  "1 DM site"   = "#4393C3",
  "≥2 DM sites" = "#D6604D"
)

# Plot
p_prob <- ggplot(dm_long, aes(x = group, y = pct, fill = dm_class)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(
    data = totals_dm,
    aes(x = group, y = 102, label = paste0("n = ", total)),
    inherit.aes = FALSE,
    size = 4, fontface = "bold"
  ) +
  scale_fill_manual(values = dm_colors, name = "DM sites per gene") +
  scale_y_continuous(
    limits = c(0, 108),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "DM site multiplicity",
    subtitle = "All DM genes vs DM ∩ DE genes",
    x        = NULL,
    y        = "Percentage of genes"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    axis.text.x      = element_text(size = 12)
  )

print(p_prob)



##5.2: Correlation between number of DM sites and log2FC

# Counting DM sites per gene (non-intergenic only, consistent with 2.4.a)
DM_sites_count <- as.data.frame(table(rel_DM$SYMBOL))
colnames(DM_sites_count) <- c("SYMBOL", "n_DM_sites")

# Joining expression values from the RNA-seq universe
sites_fc_df <- DM_sites_count %>%
  left_join(RNAseq_universe, by = "SYMBOL") %>%
  filter(!is.na(log2FoldChange))

# Spearman correlation: number of DM sites vs signed log2FC
sites_fc_cor <- cor.test(
  sites_fc_df$n_DM_sites,
  sites_fc_df$log2FoldChange,
  method = "spearman"
)

# Spearman correlation: number of DM sites vs magnitude of regulation
sites_abs_fc_cor <- cor.test(
  sites_fc_df$n_DM_sites,
  abs(sites_fc_df$log2FoldChange),
  method = "spearman"
)

# Adding to table
final_table <- rbind(final_table,
                     "Rho corr (sites vs |log2FC|)" = c(sites_abs_fc_cor$estimate, "/", "/", "/", "/", "/"),
                     "p-value corr (sites vs |log2FC|)" = c(sites_abs_fc_cor$p.value, "/", "/", "/", "/", "/"),
                     "Rho corr (sites vs log2FC)" = c(sites_fc_cor$estimate, "/", "/", "/", "/", "/"),
                     "p-value corr (sites vs log2FC)" = c(sites_fc_cor$p.value, "/", "/", "/", "/", "/"))

## Graph 8: Scatter plot of number of DM sites vs log2FC ##
#Magnitude
p_sites_fc_abs <- ggplot(sites_fc_df, aes(x = n_DM_sites, y = abs(log2FoldChange))) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
  geom_smooth(method = "lm") +
  labs(title = "Correlation between number of DM sites and expression: magnitude",
       x = "Number of DM sites",
       y = "|log2FoldChange|") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))
#Direction
p_sites_fc <- ggplot(sites_fc_df, aes(x = n_DM_sites, y = log2FoldChange)) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.7) +
  geom_smooth(method = "lm") +
  labs(title = "Correlation between number of DM sites and expression: direction",
       x = "Number of DM sites",
       y = "log2FoldChange") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

plot_grid(p_sites_fc_abs, p_sites_fc, nrow = 2, ncol = 1)


dev.off()