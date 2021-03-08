#Differentical Abundant Analysis Exploration  

library(Seurat)
library(MAST)
library(SingleCellExperiment)
library(monocle3)
library(scater)
library(DESeq2)
#reference: https://github.com/sansomlab/tenx/commit/d2e866d1c6e6ca23dcce9a481e6cbf79a1092773


setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/DEG/")
load("./../after 02 analysis/pbmc_and_mono.RData")

pbmc.sce <- as.SingleCellExperiment(pbmc.integrated)
p1 <- plotExpression(pbmc.sce, features = "MS4A1", x = "celltype") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1
p2 <- plotUMAP(pbmc.sce, colour_by = "celltype2")
p2

aheatmap(assay(scaRaw[1:1000,]), labRow='', annCol=as.data.frame(colData(scaRaw)[,c('condition', 'ourfilter')]), distfun='spearman')

mono.it@assays$SCT[1:10, 1:10]

DefaultAssay(mono.it) <- "RNA"
mono.it <- NormalizeData(mono.it,normalization.method = "LogNormalize",assay="RNA",scale.factor=10000)
mono.it <- ScaleData(mono.it)
Idents(mono.it) <- "celltype.pah"

test <- FindMarkers(mono.it, ident.1 = "X1_CD14Mono_CONT", ident.2 = "X1_CD14Mono_PAH", verbose = T, min.pct = 0,
                    test.use = "DESeq2", assay="RNA", pseudocount.use = T)

RNARowcounts <- as.data.frame(t(as.data.frame(mono.it@assays$RNA@counts)))
RNAdata <- as.data.frame(t(as.data.frame(mono.it@assays$RNA@data)))
SCTRowcounts <- as.data.frame(t(as.data.frame(mono.it@assays$SCT@counts)))
SCTdata <- as.data.frame(t(as.data.frame(mono.it@assays$SCT@data)))

integrateddata <- as.data.frame(t(as.data.frame(mono.it@assays$integrated@data)))

plot(rowSums(RNARowcounts), rowSums(integrateddata))

RNARowcounts[1:10, 1:10]
RNAdata[1:10, 1:10]

rowSums(RNARowcounts)
rowSums(RNAdata)

p <- ggplot(RNARowcounts, aes(x=rowSums(RNARowcounts), y=rowSums(RNAdata), color=mono.it@meta.data$condition)) + geom_point()
p

p <- ggplot(RNAdata, aes(x=rowSums(RNAdata), y=(RNAdata$CD14), color=mono.it@meta.data$condition)) + geom_point()

p <- ggplot(RNAdata, aes(x=RNAdata$FCGR3A, y=(RNAdata$CD14), color=mono.it@meta.data$condition)) + geom_point()
p <- p + facet_wrap(.~mono.it@meta.data$condition)
p

mono.it@assays$RNA@counts[1:10, 1:10]
mono.it@assays$RNA@data[1:10, 1:10]

mono.it@assays$SCT@counts[1:10, 1:10]
mono.it@assays$SCT@data[1:10, 1:10]
