library(data.table)
library(methylKit)
library(cowplot)
library(ggplot2)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomeInfoDb)
library(karyoploteR)


#Connecting to snakeake inputs and outputs:
files <- unlist(snakemake@input)
out_csv <- snakemake@output[[1]]
out_pdf <- snakemake@output[[2]]
out_session <- snakemake@output[[3]]
out_universe <- snakemake@output[[4]]
out_pdf_2 <- snakemake@output[[5]]

#Opening pdf with graphs
pdf(out_pdf, width = 12, height = 8)

### 1. File downloading, and MethylKit object building ###

#The data files are not in a typical bismark output format, but they were processed and simplified.
#Files only have 4 columns: "chromosome", "position", "counts M", "counts M+U (coverage)".

#Naming samples:
sample_names = c("rep_1_TB", "rep_1_CTRL", "rep_2_TB", "rep_2_CTRL", 
                 "rep_3_TB", "rep_3_CTRL", "rep_4_TB", "rep_4_CTRL", 
                 "rep_5_TB", "rep_5_CTRL", "rep_6_TB", "rep_6_CTRL")


##Using MethylSeq.This tool expects a different format, so creating the methylKit object:
methyl_obj <- new("methylRawList",                          #methylraw class is defined in the methylKit library
                  lapply(seq_along(files), function(i) {
                    df <- fread(files[i], col.names = c("chr", "pos", "M", "M+U"))
                    # removing spike-ins
                    df <- df[grepl("^chr", chr)]
                    new("methylRaw",
                        data.frame(chr = df$chr,
                                   start = df$pos,
                                   end = df$pos,          #all(start(BS_obj) == end(BS_obj)) is TRUE --> start and end are the same, cause it's only a position.
                                   strand = "*",               #No strand specifics in the original files
                                   coverage = df$'M+U',
                                   numCs = df$M,
                                   numTs = df$'M+U' - df$M
                        ),
                        sample.id = sample_names[i],
                        assembly = "hg19",                    #Data are from 2015. They used hg19 ref
                        context = "CpG",
                        resolution = "base")
                  }
                  ),
                  treatment = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0)    # 1 = "TB", 0 = "CTRL"
)


#Stats visualization
#Methylation
par(mfrow = c(3, 4))
for (i in 1:12) {
  getMethylationStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}

#Coverage
par(mfrow = c(3, 4))
for (i in 1:12) {
  getCoverageStats(methyl_obj[[i]], plot = TRUE, both.strands = FALSE)
}

#Resetting layout to normal
par(mfrow = c(1, 1))


##Filtering
#Analyzing the data:
# 1. Checking means of coverage
sapply(methyl_obj, function(x) mean(getData(x)$coverage)) 
# 2. Percentiles
cov_max <- do.call(pmax, lapply(methyl_obj, function(x) getData(x)$coverage))
quantile(cov_max, c(0.99, 0.999, 0.9999, 0.99999))

#Based on these quantiles, the coverage filter must be decided. STandard is 10, but if the medians are lower it must be changed
#PCR artifacts abundance will instead set ceiling threshold.

filtered_methyl_obj=filterByCoverage(methyl_obj,lo.count=5,lo.perc=NULL,
                                     hi.count=NULL,hi.perc=99.99)

rm(methyl_obj)
##Merging
meth=unite(filtered_methyl_obj, destrand=FALSE)

rm(filtered_methyl_obj)
##Explorative analysis on the merged:
clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)
PCASamples(meth)

#If clustering doesn't follow the expected paired design, a batch effect may be present.




### 3. Differential Analysis ###
covariates <- data.frame(replicate = factor(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)))

#Doing the analysis with the correction for overdispersion: "MN".
myDiff <- calculateDiffMeth(meth,
                            covariates = covariates,
                            overdispersion = "MN",
                            test = "Chisq")            #With correction the default is the F test: must force this for the comparison


##Finally: selecting differentially methylated bases:
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)


##Volcano plot
myDiff_df <- as.data.frame(as(myDiff, "GRanges"))
#Too many sites for the plot. Filtering high qvalue and low meth.diff values:
myDiff_df_plot     <- myDiff_df[!is.na(myDiff_df$qvalue) & 
                                  myDiff_df$qvalue < 0.6 & 
                                  abs(myDiff_df$meth.diff) > 2, ]


print (ggplot(myDiff_df_plot, aes(x = meth.diff, y = -log10(qvalue))) +
         geom_point(alpha = 0.5) +
         geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
         geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "blue") +
         geom_point(data = myDiff_df_plot[!is.na(myDiff_df_plot$qvalue) & 
                                     myDiff_df_plot$qvalue < 0.01 & 
                                     abs(myDiff_df_plot$meth.diff) > 25, ], color = "red") +
         theme_minimal() +
         labs(title = "Volcano plot: infection", x = "meth.diff", y = "-log10(adj pvalue)")
)



### 4. Annotation ###

## Before doing the annotation, it is important to note that the methyl calling was obtained with the h19 genome
## But the RNA-Seq data quantification was done using the h38. FOr consistency it is good to switch to the h38 genome.
#Doing this with liftOver()

#Download chain file
chain <- import.chain("02_wgbs/references/hg19ToHg38.over.chain")

#Converting to GRanges and liftover
myDiff25p_GR <- as(myDiff25p, "GRanges")
seqlevelsStyle(myDiff25p_GR) <- "UCSC"
myDiff25p_GR_hg38 <- unlist(liftOver(myDiff25p_GR, chain))

#Now using ChipSeeker (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#Annotation:
peakAnno <- annotatePeak(myDiff25p_GR_hg38,
                         tssRegion = c(-2000, 200),   #This is the standard definition for Promoter. Can be arbitrarly changed
                         TxDb      = txdb,
                         annoDb    = "org.Hs.eg.db")

## Visualizations
#barplot
print(plotAnnoBar(peakAnno) +
        ggtitle("DM sites annotation - barplot"))
#pie chart
plotAnnoPie(peakAnno)
#Distance plot
print(plotDistToTSS(peakAnno) +
        ggtitle("Distance of DM sites to nearest TSS"))


## Graph: Position of methylation sites on all chromosomes
## Checking if their position is clusterized around centromeres. Then considering filtering

# Recreating the GRanges object, this time with gene names
myDiff_25_GR_full <- makeGRangesFromDataFrame(peakAnno_df,
                                              keep.extra.columns = TRUE,
                                              seqnames.field = "seqnames",
                                              start.field    = "start",
                                              end.field      = "end")
seqlevelsStyle(myDiff_25_GR_full) <- "UCSC"

dev.off()
pdf(out_pdf_2, width = 12, height = 8)

kp <- plotKaryotype(genome="hg38")
kp <- kpPlotDensity(kp, myDiff_25_GR_full)

dev.off()


#Creating the final dataframe for export:
GR_df <- as.data.frame(myDiff25p_GR_hg38)
peakAnno_df <- as.data.frame(peakAnno)

final_df <- cbind(                                                             #order is garanteed, as peakAnno comes from mydiff25p_GR
  GR_df[, c("seqnames", "start", "end", "pvalue", "qvalue", "meth.diff")],
  peakAnno_df[, c("annotation", "distanceToTSS", "SYMBOL", "ENSEMBL", "geneLength", "geneStrand")]
)

#Exporting results
write.csv(final_df, out_csv, row.names = FALSE)



## Universe annotation
#Annotating all the genes considered for the analyses (from "myDiff" object)

#Again, liftover from hg19 to hg38
myDiff_GR <- as(myDiff, "GRanges")
seqlevelsStyle(myDiff_GR) <- "UCSC"
myDiff_GR_hg38 <- unlist(liftOver(myDiff_GR, chain))

#Now the annotation
myDiff_Anno  <- annotatePeak(myDiff_GR_hg38,
                             tssRegion = c(-2000, 200),
                             TxDb      = txdb,
                             annoDb    = "org.Hs.eg.db")

myDiff_Anno_df <- as.data.frame(myDiff_Anno)
WGBS_universe <- na.omit(unique(data.frame("SYMBOL" = myDiff_Anno_df$SYMBOL)))

#Exporting the universe:
write.csv(WGBS_universe, out_universe, row.names = FALSE)


### Obtaining session info:
sink(out_session)
sessionInfo()
sink()