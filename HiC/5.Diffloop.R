library(diffloop)
library(diffloopdata)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)

#use loop_counts.bedpe data
setwd("diffloop_NTC_HOXD13")
bed_dir <- list.files(path = "/path/to/folder/containing/hichipper/bedpe")
bed_dir
samples <- c("A375-shH2-K27ac", "A375-shH3-K27ac", "A375-shNTC-K27ac1", "A375-shNTC-K27ac2")
full <- loopsMake(bed_dir, samples)
celltypes <- c("A375-shHOXD13", "A375-shHOXD13", "A375-shNTC", "A375-shNTC")
full <- updateLDGroups(full, celltypes)
head(full, 8)

# remove and loops that merged together from import
full <- subsetLoops(full, full@rowData$loopWidth >= 5000)

cm <- full@counts
k_dis <- ((cm[,3]>=5 &cm[,4]==0)|(cm[,4]>=5&cm[,3]==0))
m_dis <- ((cm[,4]>=5 &cm[,4]==0)|(cm[,4]>=5&cm[,4]==0))
qc_filt <- subsetLoops(full, !(k_dis | m_dis))
qc_filt <- removeSelfLoops(qc_filt)
qc_filt <- filterLoops(qc_filt, width = 0, nreplicates = 2, nsamples = 1)
dim(qc_filt)

pdf(file.path(bed_dir, "p1.pdf"))
p1 <- loopDistancePlot(qc_filt)
p1
dev.off()

loopMetrics(qc_filt)

pcp1dat <- qc_filt
pcp1dat@colData$sizeFactor <- 1

pdf(file.path(bed_dir, "PCA1.pdf"))
pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) +
  scale_x_continuous(limits = c(-500, 500)) + ggtitle("PC Plot with no Size Factor Correction") +
  theme(legend.position="none")
pcp1
dev.off()

pdf(file.path(bed_dir, "PCA2.pdf"))
pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)) +
  scale_x_continuous(limits = c(-500, 500)) + ggtitle("PC Plot with Size Factor Correction") +
  theme(legend.position="none")
pcp2
dev.off()

#Differential Loop Calling
km_res <- quickAssoc(qc_filt)
head(km_res)
loopMetrics(km_res)
plotTopLoops(km_res, n = 100, PValue = 0.05, FDR = 1, organism = "h", colorLoops = FALSE)
#computeBoundaryScores(km_res, samples = 0, windowSize = 5e+05)
summary <- summary(km_res)
write.csv(summary, "Diffloop_shHOXD13_vs_shNTC.csv")
saveRDS(km_res, "km_res_HOXD13.rds")
