library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(arrow)
library(AUCell)
library(GSEABase)
library(harmony)
library(ggpubr)
library(dittoSeq)
setwd("ST_folder")
DirRes <-"ST_folder"


# Pozniak
Melanoma2_1 <- Load10X_Spatial('Visium_outs/Melanoma2_1/outs', 
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "tissue_lowres_image.png")
Melanoma2_1 <- SCTransform(Melanoma2_1, assay = "Spatial", verbose = FALSE)
Melanoma2_1@meta.data$sample_ID <- 'Melanoma2_1'

Pt09_K3000 <- Load10X_Spatial('Visium_outs/Pt09_K3000/outs', 
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "tissue_lowres_image.png")
Pt09_K3000 <- SCTransform(Pt09_K3000, assay = "Spatial", verbose = FALSE)
Pt09_K3000@meta.data$sample_ID <- 'Pt09_K3000'

Pt10_K28284 <- Load10X_Spatial('Visium_outs/Pt10_K28284/outs', 
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "tissue_lowres_image.png")
Pt10_K28284 <- SCTransform(Pt10_K28284, assay = "Spatial", verbose = FALSE)
Pt10_K28284@meta.data$sample_ID <- 'Pt10_K28284'

Pt11_K28439 <- Load10X_Spatial('Visium_outs/Pt11_K28439/outs', 
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "tissue_lowres_image.png")
Pt11_K28439 <- SCTransform(Pt11_K28439, assay = "Spatial", verbose = FALSE)
Pt11_K28439@meta.data$sample_ID <- 'Pt11_K28439'

Pt12_K29368 <- Load10X_Spatial('Visium_outs/Pt12_K29368/outs', 
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "tissue_lowres_image.png")
Pt12_K29368 <- SCTransform(Pt12_K29368, assay = "Spatial", verbose = FALSE)
Pt12_K29368@meta.data$sample_ID <- 'Pt12_K29368'

met <- merge(Melanoma2_1, y = c(Pt09_K3000, Pt10_K28284, Pt11_K28439, Pt12_K29368))
met <- SCTransform(met, assay = "Spatial", verbose = T)
met <- RunPCA(met, assay = "SCT", verbose = FALSE)
met <- RunHarmony(met, group.by.vars = "sample_ID", assay.use="SCT")
met <- FindNeighbors(met, reduction = "harmony", dims = 1:30)
met <- FindClusters(met, reduction = "harmony", resolution = 0.5)
met <- RunUMAP(met, reduction = "harmony", dims = 1:30)
DimPlot(met)
# Define HOXD13 positivity
expr <- FetchData(met, vars = "HOXD13")
met$HOXD13_pos <- ifelse(expr$HOXD13 > 0, "positive", "negative")

dittoBarPlot(
  met,
  var = "HOXD13_pos",
  group.by = "sample_ID",
  scale = "percent"
)

genes<-c("MIA", "TYR", "SLC45A2", "CDH19", "PMEL", "SLC24A5", "MAGEA6", "GJB1", "PLP1", "PRAME", "CAPN3", "ERBB3", 
         "GPM6B", "S100B", "FXYD3", "PAX3", "S100A1", "MLANA", "SLC26A2", "GPR143", "CSPG4", "SOX10", "MLPH", "LOXL4", 
         "PLEKHB1", "RAB38", "QPCT", "BIRC7", "MFI2", "LINC00473", "SEMA3B", "SERPINA3", "PIR", "MITF", "ST6GALNAC2", 
         "ROPN1B", "CDH1", "ABCB5", "QDPR", "SERPINE2", "ATP1A1", "ST3GAL4", "CDK2", "ACSL3", "NT5DC3", "IGSF8", "MBP", 
         "LINC00518", "LINC00520", "SAMMSON")
geneSets <- GeneSet(genes, setName="Tirosh_malignant")
geneSets
cells_rankings <- AUCell_buildRankings(met@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
Tirosh_malignant<-getAUC(cells_AUC)
Tirosh_malignant<-t(Tirosh_malignant)
met@meta.data<-cbind(met@meta.data, Tirosh_malignant)
VlnPlot(met, features = 'Tirosh_malignant', group.by = 'seurat_clusters')
met$cell_type <- ifelse(met$Tirosh_malignant > 0.18, 
                               "malignant", 
                               "non-malignant")
VlnPlot(met, features = 'Tirosh_malignant', group.by = 'cell_type')
genes<-c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG")
geneSets <- GeneSet(genes, setName="CTL")
geneSets
cells_rankings <- AUCell_buildRankings(met@assays[["SCT"]]@counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=1, assign=TRUE)
CTL<-getAUC(cells_AUC)
CTL<-t(CTL)
met@meta.data<-cbind(met@meta.data, CTL)

met_malignant <- subset(met, subset = Tirosh_malignant > 0.18)
VlnPlot(met_malignant, features = c('HOXD13', 'CTL'), group.by = 'sample_ID')
met_malignant$HOX <- ifelse(met_malignant$sample_ID %in% c('Melanoma2_1', 'Pt11_K28439', 'Pt12_K29368'), 'HOXD13high',
                  ifelse(met_malignant$sample_ID %in% c('Pt09_K3000', 'Pt10_K28284'), 'HOXD13low', NA))
VlnPlot(met_malignant, features = c('HOXD13', 'CTL'), group.by = 'HOX')  + NoLegend()
p = VlnPlot(met_malignant, features = 'HOXD13', group.by = 'HOX', y.max = 1.2)  + NoLegend()
p + stat_compare_means(
  comparisons = list(c("HOXD13high", "HOXD13low")),
  method = "t.test",
  label = "p.format"
)

p = VlnPlot(met_malignant, features = 'CTL', group.by = 'HOX', y.max = 0.4)  + NoLegend()
p + stat_compare_means(
  comparisons = list(c("HOXD13high", "HOXD13low")),
  method = "t.test",
  label = "p.format"
)

met <- merge(Melanoma2_1, y = c(Pt09_K3000))
met <- SCTransform(met, assay = "Spatial", verbose = T)

pdf('HOXD13_samples.pdf', width = 20, height = 17)
SpatialFeaturePlot(met, features = c('HOXD13', 'CTL', 'Tirosh_malignant'), 
                   keep.scale = "feature",
                   pt.size.factor = 3)
dev.off()



# Define HOXD13-positive spots
met$HOXD13_pos <- ifelse(met@assays$SCT@data["HOXD13", ] > 0, "positive", "negative")
Idents(met) <- 'HOXD13_pos'
VlnPlot(met, features = c("CTL", "HOXD13_regulon", "Tirosh_malignant"), group.by = "HOXD13_pos")

met_malignant$HOXD13_pos <- ifelse(met_malignant@assays$SCT@data["HOXD13", ] > 0, "positive", "negative")
Idents(met_malignant) <- 'HOXD13_pos'
pdf('Expression.pdf', width = 3, height = 5)
VlnPlot(met_malignant, features = "HOXD13", group.by = "HOXD13_pos", y.max = 1.3)  +
  stat_compare_means(label = "p.format", comparisons = list(c("positive", "negative")), method = "t.test") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1) + NoLegend()
VlnPlot(met_malignant, features = "CTL", group.by = "HOXD13_pos", y.max = 0.4)  +
  stat_compare_means(label = "p.format", comparisons = list(c("positive", "negative")), method = "t.test") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1) + NoLegend()
dev.off()

