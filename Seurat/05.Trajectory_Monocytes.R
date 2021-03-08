#Second round analysis May.10,2020
#PBMC analysis
#Monocle3 trajactory analysis
#work under new parameters

#Author: Xin Zhou, Ph.D. 
#Last Update:May, 21. 2020

#load necessary package
library("Seurat")
library("dplyr")
library("Matrix")
library("ggpubr")
library("DoubletFinder")
library("ggsci")
library("DESeq2")
library("ggrepel")
library("scales")
library("cowplot")
library("ggsci")
library("monocle3")
library("MAST")
options(stringsAsFactors = FALSE)

setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/")
load("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/after 02 analysis/pbmc_and_mono.RData")

DimPlot(mono.it)

Idents(mono.it) <- "celltype"
mono.nodc <- subset(mono.it, idents=c("X1_CD14Mono", "X12_CD14Mono_InterM", "X16_CD16MONO", "X24_CD14Mono_PPBP"))
DimPlot(mono.nodc)

####################################################################################################################################
#Only for control sample, subset 2000 sample for futher analysis
####################################################################################################################################
Idents(mono.nodc) <- "condition"
seurat  <- subset(mono.nodc, idents="CONT")
Idents(seurat) <- "celltype"
DimPlot(seurat)

seurat <- subset(seurat, cells= sample(x = colnames(seurat), size = 2104))
seurat

DimPlot(seurat)
seurat

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(seurat@assays[["integrated"]]@data@Dimnames[[2]], row.names = seurat@assays[["integrated"]]@data@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- seurat@assays[["integrated"]]@data
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


#list_cluster <- seurat@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
list_cluster <- seurat@meta.data[["celltype"]]
names(list_cluster) <- seurat@assays[["integrated"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

#cds_from_seurat@reduce_dim_aux@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]]<-seurat@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings


cds_from_seurat <- learn_graph(cds_from_seurat,use_partition = F)
cds_from_seurat <- order_cells(cds_from_seurat,reduction_method = "UMAP" ,root_pr_nodes = "Y_5")

p1 <- plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_roots=T,
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=4)

#cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 50)
#cds_from_seurat <- reduce_dimension(cds_from_seurat,preprocess_method="PCA",reduction_method = "UMAP")
#cds_from_seurat <- cluster_cells(cds_from_seurat, reduction_method ="UMAP",cluster_method = "louvain")

cluster_control <- plot_cells(cds_from_seurat, color_cells_by = "pseudotime")
cluster_control <- cluster_control + ggtitle("Control")
cluster_control

####################################################################################################################################
# Only for PAH sample: remove the outliers in UMAP
####################################################################################################################################
Idents(mono.nodc) <- "condition"
seurat  <- subset(mono.nodc, idents="PAH")
Idents(seurat) <- "celltype"

DimPlot(seurat)
seurat

umap <- as.data.frame(seurat@reductions[["umap"]]@cell.embeddings)
umap[umap$UMAP_2>0,]
pahcarcode  <- colnames(seurat)
subsettouse <- pahcarcode[pahcarcode!="ATACTTCGTTCTTCAT-8"]
seurat <- subset(seurat, cells = subsettouse)

DimPlot(seurat)
seurat

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(seurat@assays[["integrated"]]@data@Dimnames[[2]], row.names = seurat@assays[["integrated"]]@data@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- seurat@assays[["integrated"]]@data
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,cell_metadata = cell_metadata, gene_metadata = gene_annotation)

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- seurat@meta.data[["celltype"]]
names(list_cluster) <- seurat@assays[["integrated"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

#cds_from_seurat@reduce_dim_aux@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings

cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]]<-seurat@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings


cds_from_seurat <- learn_graph(cds_from_seurat,use_partition = F)
cds_from_seurat <- order_cells(cds_from_seurat,reduction_method = "UMAP",root_pr_nodes = "Y_16")
cds_from_seurat <- order_cells(cds_from_seurat,reduction_method = "UMAP")

p2 <- plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)

#cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 50)
#cds_from_seurat <- reduce_dimension(cds_from_seurat,preprocess_method="PCA",reduction_method = "UMAP")
#cds_from_seurat <- cluster_cells(cds_from_seurat, reduction_method ="UMAP",cluster_method = "louvain")

cluster_pah <- plot_cells(cds_from_seurat, color_cells_by = "pseudotime")
cluster_pah <- cluster_pah + ggtitle("PAH")
cluster_pah
####################################################################################################################################

pcluster <- DimPlot(mono.nodc, group.by = "celltype", split.by = "condition")
p <- cluster_control/cluster_pah/pcluster
p

ggsave(filename = "~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/trajactory/trajactory_final.pdf", p, width=8, height = 10, dpi=300)

DefaultAssay(mono.nodc) <- "integrated"
mono.nodc <- RunPCA(mono.nodc, dim=1:15)
mono.nodc <- RunUMAP(mono.nodc, dim=1:30)
mono.nodc <- FindNeighbors(mono.nodc, dims = 1:10)
mono.nodc <- FindClusters(mono.nodc,resolution = 0.8) 
DimPlot(mono.nodc)

FeaturePlot(mono.nodc, features =c("CD14", "FCGR3A","IFITM3","PPBP","HLA-DRA","STAT1"))

DefaultAssay(mono.nodc) <- "SCT"
#AHR also interesting
pmarker <- FeaturePlot(mono.nodc, features =c("FOS","JUN","JUNB", "DUSP1"), split.by = "condition", order = T)
pmarker
ggsave(filename = "~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/trajactory/fos.markers.pdf", pmarker, scale = 1)

FeaturePlot(mono.nodc, features =c("HLA-B"), split.by = "condition", order = T)

pseudotime <- plot_cells(cds_from_seurat, color_cells_by = 'pseudotime')

ciliated_genes <- c("CD14", "FCGR3A","IFITM3","PPBP","HLA-DRA","STAT1")

#heatmap
