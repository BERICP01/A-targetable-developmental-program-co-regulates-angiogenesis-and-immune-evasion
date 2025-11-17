library(Seurat)
library(AUCell)
library(GSEABase)
library(DESeq2)
library(dittoSeq)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(EnhancedVolcano)
library(Azimuth)
library(SeuratData)
library(patchwork)
options(future.globals.maxSize = 10000 * 1024^3)
setwd("B16_scRNA")
DirRes <- "B16_scRNA"

NCDneg <- Read10X(data.dir = 'pipseeker/NTC-CD45-NEG/raw_matrix')
NCDneg <- CreateSeuratObject(counts = NCDneg, min.cells = 10, min.features = 1000, project = "shNTC_CD45neg")
NCDneg@meta.data$sample <- "shNTC_CD45neg"
NCDneg@meta.data$group <- "shNTC"
NCDpos <- Read10X(data.dir = 'pipseeker/NTC-CD45-POS/raw_matrix')
NCDpos <- CreateSeuratObject(counts = NCDpos, min.cells = 10, min.features = 1000, project = "shNTC_CD45pos")
NCDpos@meta.data$sample <- "shNTC_CD45pos"
NCDpos@meta.data$group <- "shNTC"
HCDneg <- Read10X(data.dir = 'pipseeker/SH-HOX-CD45-NEG/raw_matrix')
HCDneg <- CreateSeuratObject(counts = HCDneg, min.cells = 10, min.features = 1000, project = "shHoxd13_CD45neg")
HCDneg@meta.data$sample <- "shHoxd13_CD45neg"
HCDneg@meta.data$group <- "shHoxd13"
HCDpos <- Read10X(data.dir = 'pipseeker/SH-HOX-CD45-POS/raw_matrix')
HCDpos <- CreateSeuratObject(counts = HCDpos, min.cells = 10, min.features = 1000, project = "shHoxd13_CD45pos")
HCDpos@meta.data$sample <- "shHoxd13_CD45pos"
HCDpos@meta.data$group <- "shHoxd13"
combined <- merge(NCDneg, c(NCDpos, HCDneg,HCDpos), merge.data = TRUE)
combined <- PercentageFeatureSet(combined, pattern = "^mt-", col.name = "percent.mt")
combined <- SCTransform(combined, vars.to.regress = "percent.mt", verbose = T, conserve.memory = T)
combined <- RunPCA(combined, verbose = T)
pdf(file.path(DirRes, "1.QC.pdf"), width = 14, height = 7)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = 'sample', pt.size = 0, raster=FALSE)
VlnPlot(combined, features = c("Hoxd13"), group.by = 'sample', pt.size = 0.1, raster=FALSE)
ElbowPlot(combined, ndims = 50, reduction = "pca")
dev.off()
combined <- RunUMAP(combined, dims = 1:30, reduction = "pca")
pdf(file.path(DirRes, "2.QC.pdf"), width = 14, height = 7)
DimPlot(combined, group.by = 'sample', pt.size = 0.02, label=FALSE, raster=FALSE)
dev.off()
saveRDS(combined, 'combined_1.rds')

pdf(file.path(DirRes, "3.QC.pdf"), width = 14, height = 7)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE, pt.size = 0, group.by = 'sample')
plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE, group.by = 'sample')
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE, group.by = 'sample')
CombinePlots(plots = list(plot1, plot2))
combined <- subset(combined, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 15)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE, pt.size = 0, group.by = 'sample')
dev.off()
# CELL CYCLE and clustering
CellCycle <- readRDS("mouse_cell_cycle_genes.rds")
s.genes <- CellCycle$s.genes
g2m.genes <- CellCycle$g2m.genes
combined  <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
combined <- SCTransform(combined, vars.to.regress = c('percent.mt', 'S.Score', 'G2M.Score'), verbose = T, conserve.memory = T)
combined  <- FindNeighbors(combined , dims = 1:30, reduction = "pca")
combined  <- FindClusters(combined , resolution = 0.9, reduction = "pca")
combined  <- RunUMAP(combined , dims=1:30, reduction = "pca")
pdf(file.path(DirRes,"4.UMAPs.pdf"), width = 10, height = 7)
DimPlot(combined, group.by = c('seurat_clusters'), pt.size = 0.02, label=T, raster=FALSE)
DimPlot(combined , reduction = "umap", group.by = 'group', raster=FALSE)
DimPlot(combined , reduction = "umap", group.by = 'Phase', raster=FALSE)
dittoBarPlot(combined, "seurat_clusters", group.by = 'group', scale = "percent", x.reorder = c(2,1))
dev.off()
# find markers for every cluster compared to all remaining cells, report only the positive ones
combined <- PrepSCTFindMarkers(combined, assay = "SCT", verbose = TRUE)
combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.5)
combined.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.csv(combined.markers, file = "5.combined.markers.csv")

# Subset
#Tirosh-Karras signature
genes<-c("Mia","Tyr","Slc45A2","Cdh19","Pmel","Slc24A5",
         "Magea6","Gjb1","Plp1","Prame","Capn3","Erbb3",
         "Gpm6b","S100b","Fxyd3","Pax3","S100a1","Mlana",
         "Slc26A2","Gpr143","Cspg4","Sox10","Mlph","Loxl4",
         "Plekhb1","Rab38","Qpct","Birc7","Mfi2","Sema3b",
         "Serpina3","Pir","Mitf","St6Galnac2","Ropn1B",
         "Cdh1","Abcb5","Qdpr","Serpine2","Atp1A1","St3Gal4",
         "Cdk2","Acsl3","Nt5Dc3","Igsf8","Mbp","Copg2","Cd59a",
         "Nceh1","Gjc3","Cers4","Sort1","Rapgef4","Kcnn4","Akr1b7",
         "Syngr1")
geneSets <- GeneSet(genes, setName="TK_malignant")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
TK_malignant<-getAUC(cells_AUC)
TK_malignant<-t(TK_malignant)
combined@meta.data<-cbind(combined@meta.data, TK_malignant)

#T-cells signature (Bagaev)
genes<-c("Cd28","Cd3d","Cd3e","Cd3g","Cd5","Itk","Tbx21","Trat1")

geneSets <- GeneSet(genes, setName="Tcells")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Tcells<-getAUC(cells_AUC)
Tcells<-t(Tcells)
combined@meta.data<-cbind(combined@meta.data, Tcells)

#NK-cells signature (Bagaev)
genes<-c("Cd160","Cd226","Cd244","Eomes","Fgfbp2","Gnly","Gzmb","Gzmh","Ifng","Kir2dl4","Klrc2","Klrf1","Klrk1","Ncr1","Ncr3","Nkg7","Sh2d1B")

geneSets <- GeneSet(genes, setName="NKcells")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
NKcells<-getAUC(cells_AUC)
NKcells<-t(NKcells)
combined@meta.data<-cbind(combined@meta.data, NKcells)

#B-cells signature (Tirosh)
genes<-c("Cd19","Cd79a","Cd79b","Blk","Ms4a1","Bank1","Fcrl1",
         "Pax5","Cd209c","Cd22","Bcl11a","Vpreb3","H2-Ob","Stap1",
         "Niban3","Tlr6","Ralgps2","Aff3","Pou2af1","Cxcr5",
         "Plcg2","Hvcn1","Ccr6","P2rx5","Blnk","Rubcnl","Pou2f2",
         "Irf8","Fcrla","Cd37")

geneSets <- GeneSet(genes, setName="Bcells")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Bcells<-getAUC(cells_AUC)
Bcells<-t(Bcells)
combined@meta.data<-cbind(combined@meta.data, Bcells)

#Macrophages signature (Tirosh)
genes<-c("Cd163","Cd14","Csf1r","C1qc","Vsig4","C1qa","Fcer1g",
         "F13a1","Tyrobp","Msr1","C1qb","Ms4a4a","Fpr1","S100a9",
         "Igsf6","Lilrb4b","Fpr2","Siglec1","Pirb","Lyz2","Hk3",
         "Slc11a1","Csf3r","Cd300e","Pilra","Aif1","Siglece",
         "Olr1","Tlr2","C5ar1","Fcgr1","Ms4a6d","C3ar1","Hck",
         "Il4i1","Lst1","Lilra5","Csta1","Ifi30","Cd68","Tbxas1",
         "Cxcl16","Ncf2","Rab20","Ms4a7","Nlrp3","Lrrc25","Adap2",
         "Spp1","Ccr1","Tnfsf13","Rassf4","Serpina1a","Mafb","Il18",
         "Fgl2","Sirpb1a","Clec4a2","Ifi204","Fcgr3","Clec7a",
         "Slamf8","Slc7a7","Itgax","Bcl2a1a","Plaur","Slco2b1",
         "Plbd1","Apoc1","Rnf144b","Slc31a2","Ptafr","Ninj1",
         "Itgam","Cpvl","Plin2","I830077J02Rik","Ftl1","Lipa",
         "Cd86","Glul","Fgr","Gk","Tymp","Gpx1","Npl","Acsl1")

geneSets <- GeneSet(genes, setName="Macro")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Macro<-getAUC(cells_AUC)
Macro<-t(Macro)
combined@meta.data<-cbind(combined@meta.data, Macro)

#Endothelial signature (Tirosh)
genes<-c("Pecam1","Vwf","Cdh5","Cldn5","Plvap","Ecscr","Slco2a1",
         "Ccl11","Mmrn1","Myct1","Kdr","Tm4sf5","Tie1","Erg",
         "Fabp4","Cavin2","Hyal2","Flt4","Egfl7","Esam","Tek",
         "Tspan18","Emcn","Mmrn2","Adgrl4","Pde2a","Nos3","Robo4",
         "Apold1","Ptprb","Rhoj","Ramp2","Adgrf5","F2rl3","Jup",
         "Ackr2","Gpr146","Rgs16","Tspan7","Ramp3","Pla2g4c","Tgm2",
         "Ldb2","Prcp","Id1","Smad1","Afap1l1","Elk3","Angpt2",
         "Lyve1","Arhgap29","Il3ra","Adcy4","Tfpi","Tnfaip1","Syt15",
         "Dysf","Podxl","Sema3a","Dock9","F8","Npdc1","Tspan15",
         "Cd34","Thbd","Itgb4","Rasa4","Col4a1","Ece1","Gfod2",
         "Efna1","Nectin2","Gng11","Mall","Ppm1f","Pkp4","Lims1",
         "Cd9","Rai14","Zfp521","Rgl2","Hspg2","Tgfbr2","Rbp1",
         "Fxyd6","Matn2","S1pr1","Piezo1","Pdgfa","Adam15","Hapln3","App")

geneSets <- GeneSet(genes, setName="Endo")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Endo<-getAUC(cells_AUC)
Endo<-t(Endo)
combined@meta.data<-cbind(combined@meta.data, Endo)

#Pericytes signature (Tirosh)
genes<-c("Rgs5","Acta2","Pdgfrb","Gm13889","Abcc9","Notch3","Myl9",
         "Epas1","Des","Kcnj8","Myh11","Ppp1r14a","Gucy1a1","Ndufa4l2",
         "Higd1b","Flt1","Gjc1","Heyl","Vtn","Adra2a","Aoc3","Cox4i2",
         "Adgrf5","Mustn1","Gucy1b1","Synpo2","Pla1a","Ano1","Tspan15",
         "Lmod1","Susd2","Trarg1","Trpc6","Tmem178","Esam","Irag1")

geneSets <- GeneSet(genes, setName="Peri")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Peri<-getAUC(cells_AUC)
Peri<-t(Peri)
combined@meta.data<-cbind(combined@meta.data, Peri)

#CAFs signature (Tirosh)
genes<-c("Fap","Col1a1","Col1a2","Col6a1","Col6a2","Cxcl14","Lum",
         "Col3a1","Dpt","Islr","Podn","Cd248","Fgf7","Mxra8","Pdgfrl",
         "Col14a1","Meg3","Sulf1","Aox1","Svep1","Lpar1","Pdgfrb","Tagln",
         "Igfbp6","Fbln1","Car12","Spock1","Tpm2","Thbs2","Fbln5","Tmem119",
         "Adam33","Prrx1","Pcolce","Igf2","Gfpt2","Pdgfra","Crispld2",
         "Cpe","F3","Mfap4","C1s2","Ptgis","Lox","Cyp1b1","Cldn11",
         "Serpinf1","Olfml3","Col5a2","Acta2","Msc","Vasn","Abi3bp",
         "Antxr1","Mgst1","C3","Palld","Fbn1","Cpxm1","Cybrd1","Igfbp5",
         "Prelp","Papss2","Mmp2","Ckap4","Ccdc80","Adamts2","Tpm1",
         "Pcsk5","Eln","Cxcl12","Olfml2b","Plac9","Rcn3","Ltbp2",
         "Nid2","Scara3","Amotl2","Tpst1","Mir100hg","Ccn2","Rarres2","Fhl2")

geneSets <- GeneSet(genes, setName="CAF")
geneSets
cells_rankings <- AUCell_buildRankings(combined@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
CAF<-getAUC(cells_AUC)
CAF<-t(CAF)
combined@meta.data<-cbind(combined@meta.data, CAF)
pdf(file.path(DirRes,"6a.signatures.pdf"), width = 20, height = 20)
FeaturePlot(combined, features = c("TK_malignant", "Ptprc", "Tcells", "NKcells",
                                   "Bcells", "Macro", "Endo", "Peri", "CAF"), order = T)
VlnPlot(combined, features = c("TK_malignant", "Ptprc", "Tcells", "NKcells",
                               "Bcells", "Macro", "Endo", "Peri", "CAF"), pt.size = 0)
dev.off()

cell_type <- c("Macrophages", "NKcells", "Macrophages", "Macrophages", "Malignant", "Tcells", "Macrophages", "DC", "Malignant", "Macrophages", "Malignant", "Tcells", "Malignant",
  "CAF", "Macrophages","Macrophages","Malignant","DC", "Tcells","Macrophages","Malignant","Pericytes","DC","Tcells","Macrophages","Tcells","Bcells","Endothelial","Macrophages","Monocytes","Macrophages","DC")
combined$cell_type <- factor(combined$seurat_clusters, 
                                 levels = 0:31, 
                                 labels = cell_type)
pdf(file.path(DirRes,"6b.UMAP_final_annotations.pdf"), width = 10, height = 7)
DimPlot(combined, group.by = c('cell_type'), pt.size = 0.02, label=T, raster=FALSE)
dittoBarPlot(combined, "cell_type", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Cell types')
dev.off()
freq_table <- dittoBarPlot(combined, "cell_type", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Cell types')
head(freq_table)
write.csv(freq_table$data, "celltype_percentages.csv", row.names = FALSE)

saveRDS(combined, "combined_2.rds")

#Refine malignant cells
Malignant <- subset(combined, idents = c('4','8','10','12','16', '20'))
Malignant <- subset(Malignant, TK_malignant > 0.05)
Malignant <- subset(Malignant, Ptprc == 0)
Malignant <- SCTransform(Malignant, vars.to.regress = c('percent.mt', 'S.Score', 'G2M.Score'), verbose = T, conserve.memory = T)
Malignant <- RunPCA(Malignant)
Malignant  <- FindNeighbors(Malignant, dims = 1:30, reduction = "pca")
Malignant  <- FindClusters(Malignant, resolution = 0.2, reduction = "pca")
Malignant  <- RunUMAP(Malignant, dims=1:30, reduction = "pca")
# Define an order of cluster identities
my_levels <- c('shNTC', 'shHoxd13')
Malignant$group <- factor(Malignant$group, levels = my_levels)

pdf(file.path(DirRes,"7.UMAPs_malignant.pdf"), width = 10, height = 7)
DimPlot(Malignant, group.by = c('seurat_clusters'), pt.size = 0.02, label=T, raster=FALSE)
DimPlot(Malignant, group.by = c('seurat_clusters'), split.by = 'group', pt.size = 0.02, label=F, raster=FALSE)
dittoBarPlot(Malignant, "seurat_clusters", group.by = 'group', scale = "percent", x.reorder = c(2,1))
dittoBarPlot(Malignant, "Phase", group.by = 'group', scale = "percent", x.reorder = c(2,1))
VlnPlot(Malignant, features = c("nFeature_RNA", "nCount_RNA"), pt.size = 0, group.by = 'group')
VlnPlot(Malignant, c('Hoxd13', 'Vegfa', 'Sema3a', 'Nt5e'), raster=FALSE, group.by = 'group', ncol = 2)
FeatureScatter(Malignant, 'Hoxd13', 'nFeature_RNA', split.by = 'group')
FeatureScatter(Malignant, 'Hoxd13', 'TK_malignant', split.by = 'group')
dev.off()
Idents(Malignant) <- 'group'
DE <- FindMarkers(Malignant, group.by = "group", ident.1 = "shHoxd13", ident.2 = "shNTC", 
                  assay = "SCT", logfc.threshold = 0.5)

p = VlnPlot(Malignant, features = c('Hoxd13', 'Vegfa', 'Sema3a', 'Nt5e'), 
            group.by = 'group', y.max = 4, ncol = 2) + NoLegend()
pdf(file.path(DirRes,"8.DE_malignant_volcano.pdf"))
EnhancedVolcano(DE, 
                rownames(DE),
                x ="avg_log2FC", 
                y ="p_val_adj", title = "Malignant shNTC vs shHoxd13")
p & stat_compare_means(
  comparisons = list(c("shNTC", "shHoxd13")),
  method = "t.test",
  label = "p.format"
)
dev.off()

# Malignant states from Karras
# Melanocytic
genes<-c("Mlph","Mitf","Mlana","Pmel","Slc45a2","Apoe","Dct","Gpnmb","Tyr","Cited1","Bace2","Cox17","Cox7a2","Ndufb2","Ndufb4","Uqcr10","Uqcrb")

geneSets <- GeneSet(genes, setName="Melanocytic")
geneSets
cells_rankings <- AUCell_buildRankings(Malignant@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Melanocytic<-getAUC(cells_AUC)
Melanocytic<-t(Melanocytic)
Malignant@meta.data<-cbind(Malignant@meta.data, Melanocytic)

# NeuralCrest
genes<-c("Tfap2b", "Prnp", "Mef2c", "Gfra1", "Cd200", "Syt11", "Thsd7a", "Cxxc4", 
         "Sema5a", "Tbx3", "Kif26b", "Efhd1", "Neto2", "Hmcn1", "Igsf10", "Olfml3", 
         "Rgs2", "Nt5e", "Morc4", "Aqp1", "Gfra2", "Wnt4", "Abcg2", "Elovl5", "Emilin1", "Fibin")

geneSets <- GeneSet(genes, setName="NeuralCrest")
geneSets
cells_rankings <- AUCell_buildRankings(Malignant@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
NeuralCrest<-getAUC(cells_AUC)
NeuralCrest<-t(NeuralCrest)
Malignant@meta.data<-cbind(Malignant@meta.data, NeuralCrest)

# Mesenchymal
genes<-c("Fap", "Lama2", "Prrx1", "Slit3", "Abi3bp", "Loxl1", "Cdh11", 
         "Col6a3", "Col6a2", "Loxl2", "Mfap5", "Fbn1", "Col4a2", "Pcolce", 
         "Lum", "Col4a1", "Col5a2", "Thy1", "Fbn2", "Bgn", "Tgfbi", "Pdgfrb", 
         "Sulf1", "Inhba", "Col1a1", "Col3a1", "Col1a2", "Ctsk", "Dcn", 
         "Serpinf1", "Sparc", "Fstl1")

geneSets <- GeneSet(genes, setName="Mesenchymal")
geneSets
cells_rankings <- AUCell_buildRankings(Malignant@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
Mesenchymal<-getAUC(cells_AUC)
Mesenchymal<-t(Mesenchymal)
Malignant@meta.data<-cbind(Malignant@meta.data, Mesenchymal)

# Antigen Presenting
genes<-c("Gbp4", "Gbp6", "Irf1", "Nlrc5", "Stat2", "Cd74", "Il12rb1", "B2m", 
         "Gbp2", "Ifit3", "Gbp3", "Psmb8", "Stat1", "Ifit1", "Isg15", "Psmb9", 
         "Gbp7", "Xaf1", "Tap1", "Gbp5", "Usp18", "Ciita", "Irf7", "Ube2l6", 
         "Tap2", "Psme1", "Ifitm3", "Psme2", "Tapbp", "Ifi35", "Eif2ak2", "Irf9", 
         "Ddx58", "Trim25", "Ifit2", "Irf8", "H2-DMa", "H2-Aa", "H2-DMb1", "H2-Eb1", 
         "H2-Ab1", "H2-T23")

geneSets <- GeneSet(genes, setName="AntigenPresenting")
geneSets
cells_rankings <- AUCell_buildRankings(Malignant@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
AntigenPresenting<-getAUC(cells_AUC)
AntigenPresenting<-t(AntigenPresenting)
Malignant@meta.data<-cbind(Malignant@meta.data, AntigenPresenting)

# Stress-like
genes<-c("Bnip3", "Tpi1", "Slc2a1", "Mif", "Vldlr", "Hk2", "Vegfa", 
         "Ldha", "Pfkl", "P4ha1", "Fam162a", "Bhlhe40", "Pgk1", 
         "Aldoa", "Pfkp", "Pdk1", "Hspa9", "Pgam1", "Pkm", "Atf4")

geneSets <- GeneSet(genes, setName="StressLike")
geneSets
cells_rankings <- AUCell_buildRankings(Malignant@assays[["SCT"]]$counts)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)
StressLike<-getAUC(cells_AUC)
StressLike<-t(StressLike)
Malignant@meta.data<-cbind(Malignant@meta.data, StressLike)

pdf(file.path(DirRes,"9a.Malignant_states.pdf"), width = 20, height = 20)
VlnPlot(Malignant, features = c("Mitf", "Sox10", "Dct", 
                                "Tyr", "Atf4", "Prrx1", "Aqp1", "Egfr", "Tcf4"), pt.size = 0, group.by = 'seurat_clusters')
FeaturePlot(Malignant, features = c("Melanocytic", "NeuralCrest", "Mesenchymal", 
                                   "AntigenPresenting", "StressLike"), order = T, cols = c("blue","lightgray", "red"))
VlnPlot(Malignant, features = c("Melanocytic", "NeuralCrest", "Mesenchymal", 
                                "AntigenPresenting", "StressLike"), pt.size = 0, group.by = 'seurat_clusters')
VlnPlot(Malignant, features = c("Melanocytic", "NeuralCrest", "Mesenchymal", 
                               "AntigenPresenting", "StressLike"), pt.size = 0, group.by = 'group', y.max = 1) & stat_compare_means(
                                 comparisons = list(c("shNTC", "shHoxd13")),
                                 method = "t.test",
                                 label = "p.format")
dev.off()
Idents(Malignant) <- 'seurat_clusters'
Malignant<- PrepSCTFindMarkers(Malignant, assay = "SCT", verbose = TRUE)
Malignant.markers <- FindAllMarkers(Malignant, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.5)
Malignant.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.csv(Malignant.markers, file = "9b.Malignant.markers.csv")
melanoma_states <- c("Stress-like", "Antigen_presenting", "Melanocytic", 
                     "Protein_processing", "RNA_processing")
Malignant$melanoma_state <- factor(Malignant$seurat_clusters, 
                                 levels = 0:4, 
                                 labels = melanoma_states)
pdf(file.path(DirRes,"9c.Malignant_UMAPs_final.pdf"), width = 10, height = 7)
DimPlot(Malignant, group.by = c('melanoma_state'), pt.size = 0.2, label=FALSE, raster=FALSE)
DimPlot(Malignant, group.by = c('melanoma_state'), split.by = 'group', pt.size = 0.2, label=FALSE, raster=FALSE)
dittoBarPlot(Malignant, "melanoma_state", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Malignant')
dev.off()
saveRDS(Malignant, "Malignant.rds")

pdf(file.path(DirRes,"9d.Malignant_some_genes.pdf"), width = 4, height = 4)
DotPlot(Malignant, 
        features = c("H2-K1","H2-D1","H2-M3","H2-Aa","H2-Ab1","H2-Eb1","H2-Eb2",
                     "H2-DMa","H2-DMb1","H2-DMb2","H2-Ob","Tap1","Tap2","Tapbp",
                     "Psmb8","Psmb9","Psmb10"), 
        group.by = "group") + 
  coord_flip() + RotatedAxis()
DotPlot(Malignant, 
        features = c("Vegfa", "Vegfb", "Vegfc", "Sema3a", "Sema3b", "Sema3c", "Sema3d", "Sema3f", "Sema3g", "Sema3e"), 
        group.by = "group") + 
  coord_flip() + RotatedAxis()
dev.off()

pdf(file.path(DirRes,"9e.Malignant_some_genes.pdf"), width = 4, height = 4)
VlnPlot(Malignant, 
        features = c("Hoxd13", "Sema3a"), 
        group.by = "group", ncol = 2, y.max = 2) + NoLegend() & stat_compare_means(
          comparisons = list(c("shNTC", "shHoxd13")),
          method = "t.test",
          label = "p.format")& stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1)
VlnPlot(Malignant, 
        features = c("Vegfa", "Vegfb"), 
        group.by = "group", ncol = 2, y.max = 4) + NoLegend() & stat_compare_means(
          comparisons = list(c("shNTC", "shHoxd13")),
          method = "t.test",
          label = "p.format")& stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1)
VlnPlot(Malignant, 
        features = c("Nt5e", "Tgfb1"), 
        group.by = "group", ncol = 2, y.max = 4) + NoLegend() & stat_compare_means(
          comparisons = list(c("shNTC", "shHoxd13")),
          method = "t.test",
          label = "p.format")& stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1)
dev.off()

# Extract cell metadata
meta <- Malignant@meta.data

# Calculate raw counts (number of cells per group and predicted cell type)
raw_counts <- meta %>%
  group_by(group, melanoma_state) %>%
  summarise(n = n(), .groups = 'drop')

# Calculate percentages within each group (same as scale = "percent" in dittoBarPlot)
raw_percent <- raw_counts %>%
  group_by(group) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

# Write to CSV
write.csv(raw_percent,
          file = file.path(DirRes, "melanoma_states_dittoBarPlot_values.csv"),
          row.names = FALSE)

#Refine TME
TME <- subset(combined, idents = c('0','1','2','3','5','6','7','9','11','13','14','15','17',
                         '18','19','21','22','23','24','25','26','27','28','29','30','31'))
TME <- SCTransform(TME, vars.to.regress = c('percent.mt', 'S.Score', 'G2M.Score'), verbose = T, conserve.memory = T)
TME <- RunPCA(TME)
TME  <- FindNeighbors(TME, dims = 1:30, reduction = "pca")
TME  <- FindClusters(TME, resolution = 0.9, reduction = "pca")
TME  <- RunUMAP(TME, dims=1:30, reduction = "pca")
my_levels <- c('shNTC', 'shHoxd13')
TME$group <- factor(TME$group, levels = my_levels)
pdf(file.path(DirRes,"10.UMAPs_TME.pdf"), width = 10, height = 7)
DimPlot(TME, group.by = c('seurat_clusters'), pt.size = 0.02, label=T, raster=FALSE)
FeaturePlot(TME, features = c("Ptprc", "Tcells", "NKcells", "Bcells", "Macro", "Endo", "Peri", "CAF"), 
            order = T)
VlnPlot(TME, features = c("Ptprc", "Tcells", "NKcells", "Bcells", "Macro", "Endo", "Peri", "CAF"), 
            pt.size = 0)
dev.off()
TME<- PrepSCTFindMarkers(TME, assay = "SCT", verbose = TRUE)
TME.markers <- FindAllMarkers(TME, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.5)
TME.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.csv(TME.markers, file = "11.TME.markers.csv")
saveRDS(TME, "TME.rds")
TME <- readRDS("/gpfs/data/HernandoLab/home/bericp01/scRNAseq/HOXD13/Analysis/B16/TME.rds")
#Pericytes, Endothelial and CAF
PEC <- subset(TME, idents = c('10', '21', '26'))
PEC <- RunPCA(PEC)
PEC  <- RunUMAP(PEC, dims=1:30, reduction = "pca")
identity <- c("CAF", "Pericytes", "Endothelial cells")
PEC$cell_identity <- factor(PEC$seurat_clusters, 
                                   levels = c(10,21,26), 
                                   labels = identity)
pdf(file.path(DirRes,"12.UMAPs_CAF_Peri_Endo.pdf"))
DimPlot(PEC, group.by = 'cell_identity', label=T)+ NoLegend()
DimPlot(PEC, group.by = 'group')
dev.off()
pdf(file.path(DirRes,"12.Peri_Endo_proportion.pdf"),height = 4, width = 3)
dittoBarPlot(PEC, "cell_identity", group.by = 'group', scale = "percent", x.reorder = c(2,1))
dev.off()
Endothelial <- subset(PEC, idents = c('26'), )
Idents(Endothelial) <- 'group'
DE.endo <- FindMarkers(Endothelial, group.by = "group", ident.1 = "shHoxd13", ident.2 = "shNTC", 
                  assay = "SCT", logfc.threshold = 0.5)
write.csv(DE.endo, "13.DE_Endothelial.csv")

Pericytes <- subset(PEC, idents = c('21'))
Idents(Pericytes) <- 'group'
DE.peri <- FindMarkers(Pericytes, group.by = "group", ident.1 = "shHoxd13", ident.2 = "shNTC", 
                  assay = "SCT", logfc.threshold = 0.5)
write.csv(DE.peri, "14.DE_Pericytes.csv")

CAF <- subset(PEC, idents = c('10'))
Idents(CAF) <- 'group'
DE.caf <- FindMarkers(CAF, group.by = "group", ident.1 = "shHoxd13", ident.2 = "shNTC", 
                  assay = "SCT", logfc.threshold = 0.5)
write.csv(DE.caf, "15.DE_CAF.csv")
pdf(file.path(DirRes,"16.DE_CAF_peri_endo.pdf"))
DimPlot(PEC, group.by = 'cell_identity', label=T)+ NoLegend()
EnhancedVolcano(DE.endo, 
                rownames(DE.endo),
                x ="avg_log2FC", 
                y ="p_val_adj", title = "Endothelial cells shNTC vs shHoxd13")
EnhancedVolcano(DE.peri, 
                rownames(DE.peri),
                x ="avg_log2FC", 
                y ="p_val_adj", title = "Pericytes shNTC vs shHoxd13")
EnhancedVolcano(DE.caf, 
                rownames(DE.caf),
                x ="avg_log2FC", 
                y ="p_val_adj", title = "CAFs shNTC vs shHoxd13")
VlnPlot(Endothelial, features = c('Pecam1','Angpt2','Tie1','Tie2','Vegfa'), group.by = 'group', y.max = 5) + NoLegend() + ggtitle('Pecam1', subtitle = 'Endothelial cells') & stat_compare_means(
  comparisons = list(c("shNTC", "shHoxd13")),
  method = "t.test",
  label = "p.format")
dittoBarPlot(Endothelial, "Phase", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Endothelial cells')
VlnPlot(Pericytes, features = 'Acta2', group.by = 'group', y.max = 6) + NoLegend() + ggtitle('Acta2', subtitle = 'Pericytes') & stat_compare_means(
  comparisons = list(c("shNTC", "shHoxd13")),
  method = "t.test",
  label = "p.format")
dittoBarPlot(Pericytes, "Phase", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Pericytes')
VlnPlot(CAF, features = 'Fap', group.by = 'group', y.max = 4) + NoLegend() + ggtitle('Fap', subtitle = 'CAFs') & stat_compare_means(
  comparisons = list(c("shNTC", "shHoxd13")),
  method = "t.test",
  label = "p.format")
dittoBarPlot(CAF, "Phase", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('CAFs')
dittoBarPlot(PEC, "cell_identity", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('CAFs')
dev.off()

pdf(file.path(DirRes,"16.Endo_markers.pdf"), height = 4, width = 6)
VlnPlot(Endothelial, features = c('Angpt2','Vegfa','Cxcl9'), group.by = 'group', y.max = 5) + NoLegend() & stat_compare_means(
  comparisons = list(c("shNTC", "shHoxd13")),
  method = "t.test",
  label = "p.format")& stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1)
dev.off()

pdf(file.path(DirRes,"16.Endo_proliferation.pdf"), height = 4, width = 3)
dittoBarPlot(Endothelial, "Phase", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Endothelial cells')
dev.off()

pdf(file.path(DirRes,"16.Peri_markers.pdf"), height = 5, width = 2)
VlnPlot(Pericytes, features = c('Acta2'), group.by = 'group', y.max = 6) + NoLegend() & stat_compare_means(
  comparisons = list(c("shNTC", "shHoxd13")),
  method = "t.test",
  label = "p.format")& stat_summary(fun = mean, geom = "crossbar", width = 0.3, color = "red", fatten = 1)
dev.off()

pdf(file.path(DirRes,"16.Peri_proliferation.pdf"), height = 4, width = 3)
dittoBarPlot(Pericytes, "Phase", group.by = 'group', scale = "percent", x.reorder = c(2,1)) + ggtitle('Pericytes')
dev.off()

##Characterization of immune compartment
VlnPlot(TME, 'Ptprc')
Immune <- subset(TME, idents = c('0','1','2','3','4','5','6','7','8','9','11','12','13','14','15','16','17','18','19','20','23','24','25','27','28','29'))

#InstallData("pbmcsca")
pbmcsca <- LoadData("pbmcsca")
DefaultAssay(Immune) <- 'RNA'
#Immune@assays$SCT <- NULL

pbmcsca <- RunAzimuth(
  query = Immune,
  reference = "pbmcref"
)

pdf(file.path(DirRes, "17.Immune_annotation.pdf"), width = 7, height = 7)
DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = F)+
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )
DimPlot(pbmcsca, group.by = "group")
dittoBarPlot(pbmcsca, "predicted.celltype.l2", group.by = 'group', scale = "percent", x.reorder = c(2,1))
dev.off()

# UMAP Tcells
Idents(pbmcsca) <- 'predicted.celltype.l2'
subsetT <- subset(pbmcsca, idents = c("CD4 Naive","CD4 Proliferating","CD4 TCM","CD4 TEM","CD8 Naive","CD8 Proliferating","CD8 TCM","CD8 TEM","Treg"))
DefaultAssay(subsetT) <- "RNA"
subsetT <- NormalizeData(subsetT)
subsetT  <- FindNeighbors(subsetT, dims = 1:30, reduction = "pca")
subsetT  <- FindClusters(subsetT , resolution = 0.6, graph.name = "SCT_snn")
# Get gene names (features) from RNA assay
genes <- rownames(subsetT@assays$RNA@features)

# Save to txt
write.table(genes,
            file = "genes_in_subsetT.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

pdf(file.path(DirRes, "18.UMAP_Tcells.pdf"), width = 10, height = 7)
DimPlot(subsetT, group.by = "predicted.celltype.l2", label = F)+
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )
DimPlot(subsetT, group.by = "predicted.celltype.l2", split.by = 'group', label = F)+
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )
DimPlot(subsetT, group.by = "seurat_clusters", label = F)+
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )
DimPlot(subsetT, group.by = "group")
#FeaturePlot(subsetT, features = '', order =T)
dev.off()
combined.markers <- FindAllMarkers(subsetT, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.5)
combined.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
write.csv(combined.markers, file = "18.subsetT.markers.csv")

# Extract immune cells percentage
meta <- pbmcsca@meta.data
df <- as.data.frame(table(meta$group, meta$predicted.celltype.l2))
colnames(df) <- c("group", "predicted.celltype.l2", "count")
df <- df %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(percent = count / sum(count) * 100)
write.csv(df, file.path(DirRes, "17.Immune_annotation_percent.csv"), row.names = FALSE)

Idents(pbmcsca) <- 'predicted.celltype.l2'
pbmcsca <- NormalizeData(pbmcsca)



subsetT <- subset(pbmcsca, idents = c("CD8 Proliferating", "CD8 TCM", "CD8 TEM", "CD8 Naive"))
Idents(subsetT) <- 'group'
DE.T <- FindMarkers(subsetT, group.by = "group", ident.1 = "shHoxd13", ident.2 = "shNTC", 
                       assay = "RNA", logfc.threshold = 0.5)
write.csv(DE.T, "18.DE_Tcells.csv")
