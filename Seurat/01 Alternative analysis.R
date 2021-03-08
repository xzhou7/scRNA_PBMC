#Second round analysis Feb.23,2020
#PBMC analysis
#Set new cutoffs and compare parameters

#Author: Xin Zhou, Ph.D. 
#Last Update: March.26. 2020

library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggsci)
library(patchwork)

setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/")
load(file = "../PBMC_PAH.RData")

##add the house keeping gene list[Tirosh et al., 2016]
hkgene <- read.table("~/Box/XinZhouFiles/Projects/ZXP1_PAH/Single Cell Related/House_keeping_Gene.txt")
hkgenes <- as.vector(hkgene$V1)
hkgenes.found <- which(toupper(rownames(pbmc@assays$RNA)) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(pbmc@assays$RNA[hkgenes.found, ] > 0)
pbmc <- AddMetaData(object = pbmc, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

#add percentage of hkgenes
pbmc[["percent.RPS"]] <- PercentageFeatureSet(pbmc, pattern= "^RPS")
pbmc[["percent.RPL"]] <- PercentageFeatureSet(pbmc, pattern= "^RPL")
pbmc[["percent.PRPS"]] <- PercentageFeatureSet(pbmc, pattern= "^PRPS")
pbmc[["percent.RB"]] <- pbmc[["percent.RPS"]] + pbmc[["percent.RPL"]]+pbmc[["percent.PRPS"]]
pbmc[["percent.ACTB"]] <- PercentageFeatureSet(pbmc, pattern= "^ACTB")
pbmc[["percent.B2M"]] <- PercentageFeatureSet(pbmc, pattern= "^B2M")
pbmc[["percent.PSMB"]] <- PercentageFeatureSet(pbmc, pattern= "^PSMB")
pbmc[["percent.HNRPLL"]] <- PercentageFeatureSet(pbmc, pattern= "^HNRPLL")
pbmc[["percent.HPRT"]] <- PercentageFeatureSet(pbmc, pattern= "^HPRT")
pbmc[["percent.PPIA"]] <- PercentageFeatureSet(pbmc, pattern= "^PPIA")
pbmc[["percent.TRPS1"]] <- PercentageFeatureSet(pbmc, pattern= "^TRPS1")

pbmc[["percent.hkg"]] <- pbmc[["percent.RPS"]] + pbmc[["percent.ACTB"]] + pbmc[["percent.B2M"]] + pbmc[["percent.HNRPLL"]] + pbmc[["percent.PPIA"]] +
                          pbmc[["percent.RPL"]]+ pbmc[["percent.PRPS"]] + pbmc[["percent.PSMB"]] + pbmc[["percent.HPRT"]] + pbmc[["percent.TRPS1"]]
  
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"), ncol = 4,pt.size = 0.001)
VlnPlot(pbmc, features = c("percent.hkg", "percent.mt", "n.exp.hkgenes"), ncol = 4,pt.size = 0.001)

plot1 <- FeatureScatter(pbmc2, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size = 0.01, group = "chemistry") + facet_wrap(.~colors)
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.01,group = "chemistry") + facet_wrap(.~colors)
plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "n.exp.hkgenes",pt.size = 0.01,group = "chemistry") + facet_wrap(.~colors)
plot4 <- FeatureScatter(pbmc, feature1 = "nFeature_RNA", feature2 = "n.exp.hkgenes",pt.size = 0.01,group = "chemistry") + facet_wrap(.~colors)
plot5 <- FeatureScatter(pbmc, feature1 = "PPBP", feature2 = "n.exp.hkgenes",pt.size = 0.01,group = "chemistry") + facet_wrap(.~colors)
plot6 <- FeatureScatter(pbmc, feature1 = "percent.hkg", feature2 = "n.exp.hkgenes",pt.size = 0.01,group = "chemistry") + facet_wrap(.~colors)
plot1.zoom <- plot1 + xlim(0,5000)
plot1/plot1.zoom

#plot line for remove doublelets
p <- ggplot(pbmc[[]], aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(size=0.01, color="red")
p <- p + geom_abline(intercept = -170, slope = 0.22) + theme_cowplot()
p
#plot lines for remove by features/counts
p <- ggplot(pbmc[[]], aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point(size=0.01, color="red")
p <- p + geom_vline(xintercept=2000, color="blue") +geom_hline(yintercept=c(300,2100, 3500), color="blue") + theme_cowplot()
p

plot2

plot3

plot4

plot5

pbmc

#01 remove doublelets
pbmc <- subset(pbmc, subset =  nFeature_RNA > (0.22*nCount_RNA - 170))
pbmc
#02 remove high mt rna sample 
pbmc <- subset(pbmc, subset = (chemistry %in% "V3" & percent.mt < 20) | (chemistry %in% "V2" & percent.mt < 13))
pbmc
#03 remove high feature sample
pbmc <- subset(pbmc, subset = (chemistry %in% "V3" & nFeature_RNA < 3500) | (chemistry %in% "V2" & nFeature_RNA < 2100))
pbmc
#04 remove by count
pbmc <- subset(pbmc, subset = (chemistry %in% "V3" & nCount_RNA > 2100) | (chemistry %in% "V2" & nCount_RNA > 1100))
pbmc
#05 remove by hkgene
pbmc <- subset(pbmc, subset = n.exp.hkgenes > 50)
pbmc

#create list
pbmc.list <- SplitObject(pbmc, split.by = "chemistry")

#SCT normalization 
for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = T)
}

#merge sample together
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features,  verbose = T)

#the next step normally takes 25 minutes
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", anchor.features = pbmc.features, verbose = T)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", verbose = T)

#find varible features
pbmc.integrated <- FindVariableFeatures(pbmc.integrated, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(pbmc.integrated)
pbmc.integrated <- ScaleData(pbmc.integrated, features = all.genes)

pbmc.integrated <- RunPCA(pbmc.integrated, features = VariableFeatures(object = pbmc.integrated))
DimPlot(pbmc.integrated, reduction = "pca")
ElbowPlot(pbmc.integrated)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:15)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.3)

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:30)

pUMAP <- DimPlot(pbmc.integrated, reduction = "umap", label = T)#, group.by = "subject")
pUMAP

pbmc.integrated <- BuildClusterTree(pbmc.integrated,reorder.numeric = TRUE)
PlotClusterTree(pbmc.integrated)

table(pbmc.integrated$integrated_snn_res.1.3)

vplot <- VlnPlot(object = pbmc.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"),
                 ncol = 2, combine = T, pt.size=0.05,group.by = "integrated_snn_res.1.3")
vplot

vplot <- VlnPlot(object = pbmc.integrated, features = c("NKG7", "FOXP3", "RORC", "CD14"), ncol = 2, combine = T, pt.size=0.05)
vplot

vplot <- VlnPlot(object = pbmc.integrated, features = c("CD8A", "NKG7", "CD4"), ncol = 1, combine = T, pt.size=0.05)
vplot


#select 1.3 for downsteam analysis
Idents(pbmc.integrated) <- "integrated_snn_res.1.3"

plots1 <- DimPlot(pbmc.integrated, group.by ="integrated_snn_res.1.3",label=T)
plots1
ggsave2(filename = "./F1.UMAP.pdf", plots1, scale=0.77)

vplot <- VlnPlot(object = pbmc.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "n.exp.hkgenes"), ncol = 2, combine = T, pt.size=0.05)
vplot

grep(pattern = "^PTPRC", x = rownames(x = pbmc.integrated@assays$RNA), value = TRUE)

#plot single feature
FeaturePlot(pbmc.integrated, features ="FIT3",min.cutoff = 0, order = T)

#plot mutiple markers
markers = c("GZMB", "GNLY", "CD3D", "CCR7", "FCER1A", "LYZ", "HBB", 
            "IL7R","FCGR3A", "CST3","NKG7", "CD8A","IL3RA","SERPINF1", 
            "CD19", "CD79A", "SELL", "CD69","PPBP", "CD8B", "CD14", 
            "HSPB1", "MS4A7", "CD4", "FOXP3","GZMK", "GZMB")

FeaturePlot(pbmc.integrated, features=c(), min.cutoff = 0, order = T)

DefaultAssay(pbmc.integrated)

DefaultAssay(pbmc.integrated) <- "integrated"
for(gene in markers){
  print(gene)
  p=FeaturePlot(object = pbmc.integrated, features = c(gene), split.by="condition",cols=c("lightgrey", "#931F21"),pt.size=0.3, sort.cell=T,min.cutoff = 0)+
    theme(axis.title.x = element_text(vjust=0.5, size=24,face="bold"),axis.text.x  = element_text(vjust=0.5, size=24,face="bold"),
          axis.title.y = element_text(vjust=0.5, size=24,face="bold"),axis.text.y  = element_text(vjust=0.5, size=24,face="bold"),
          legend.text = element_text(size = 20, face = "bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave(p,filename = paste0("./integrated/features/",gene,".pdf"),height=7.5,width=17)
}

p.rnalist.1 <- FeaturePlot(pbmc.integrated, features = c("GZMB", "GNLY", "CD3D", "CCR7", "FCER1A", "LYZ", "HBB", "IL7R","FCGR3A", "CST3","NKG7", "CD8A"),
                           min.cutoff = 0, order = T)
p.rnalist.1
ggsave2("./integrated/plist1_In.png",p.rnalist.1,device = "png", dpi=300, width=25, height=15)

p.rnalist.2 <- FeaturePlot(pbmc.integrated, features = c("IL3RA","SERPINF1", "CD19", "CD79A", "SELL", "CD69","PPBP", "CD8B", "CD14", "HSPB1", "MS4A7", "CD4"), 
                           min.cutoff = 0, order = T)
p.rnalist.2
ggsave2("./integrated/plist2_In.png",p.rnalist.2,device = "png", dpi=300, width=25, height=15)


DefaultAssay(pbmc.integrated) <- "RNA"
for(gene in markers){
  print(gene)
  p=FeaturePlot(object = pbmc.integrated, features = c(gene), split.by="condition",cols=c("lightgrey", "#931F21"),pt.size=0.3, sort.cell=T)+
    theme(axis.title.x = element_text(vjust=0.5, size=24,face="bold"),axis.text.x  = element_text(vjust=0.5, size=24,face="bold"),
          axis.title.y = element_text(vjust=0.5, size=24,face="bold"),axis.text.y  = element_text(vjust=0.5, size=24,face="bold"),
          legend.text = element_text(size = 20, face = "bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave(p,filename = paste0("./RNA/features/",gene,".pdf"),height=7.5,width=17)
}

p.rnalist.1 <- FeaturePlot(pbmc.integrated, features = c("GZMB", "GNLY", "CD3D", "CCR7", "FCER1A", "LYZ", "HBB", "IL7R","FCGR3A", "CST3","NKG7", "CD8A"),
                           min.cutoff = 0, order = T)
p.rnalist.1
ggsave2("./RNA/plist1_In.png",p.rnalist.1,device = "png", dpi=300, width=25, height=15)

p.rnalist.2 <- FeaturePlot(pbmc.integrated, features = c("IL3RA","SERPINF1", "CD19", "CD79A", "SELL", "CD69","PPBP", "CD8B", "CD14", "HSPB1", "MS4A7", "CD4"), 
                           min.cutoff = 0, order = T)
p.rnalist.2
ggsave2("./RNA/plist2_In.png",p.rnalist.2,device = "png", dpi=300, width=25, height=15)

pbmc.integrated@meta.data$celltype <- "1"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==0] <- "X0_CD4_1"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==1] <- "X1_CD14Mono"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==2] <- "X2_CD4_2"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==3] <- "X3_CD4_3"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==4] <- "X4_NaiveCD8"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==5] <- "X5_CD4_4"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==6] <- "X6_CD4_5"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==7] <- "X7_D8_NKT?"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==8] <- "X8_B1"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3 ==9] <- "X9_NK1"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==10] <- "X10_CD8_1"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==11] <- "X11_B2"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==12] <- "X12_CD14Mono_InterM"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==13] <- "X13_CD4_FOXP3"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==14] <- "X14_CD4CD8"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==15] <- "X15_NK2"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==16] <- "X16_CD16MONO"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==17] <- "X17_CD8_RORC"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==18] <- "X18_CD4_UNK"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==19] <- "X19_cDC"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==20] <- "X20_NK2"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==21] <- "Megakaryocytes"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==22] <- "X22_pDC"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==23] <- "X23_CD4_UNK2"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==24] <- "X24_CD14Mono_PPBP"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==25] <- "X25_CD8_UK"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==26] <- "X26_B_UK"
pbmc.integrated$celltype[pbmc.integrated$integrated_snn_res.1.3==27] <- "X27_UK"

#identify cell type by color
pbmc.integrated@meta.data$celltype2 <- "2"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 0] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 1] <- "CD14Mono"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 2] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 3] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 4] <- "CD8_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 5] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 6] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 7] <- "CD8_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 8] <- "B"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3== 9] <- "NK"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==10] <- "CD8_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==11] <- "B"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==12] <- "CD14Mono"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==13] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==14] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==15] <- "NK"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==16] <- "CD16Mono"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==17] <- "CD8_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==18] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==19] <- "cDC"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==20] <- "NK"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==21] <- "Megakaryocytes"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==22] <- "pDC"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==23] <- "CD4_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==24] <- "CD14Mono"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==25] <- "CD8_T"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==26] <- "B"
pbmc.integrated$celltype2[pbmc.integrated$integrated_snn_res.1.3==27] <- "UK"

plots2 <- DimPlot(pbmc.integrated, group.by ="celltype2", label = T, label.size = 5) + scale_color_npg()
plots2

plots3 <- DimPlot(pbmc.integrated, group.by ="celltype2", label = F,split.by = "subject",ncol=4, label.size = 5) + scale_color_npg()
plots3

plots3 <- DimPlot(pbmc.integrated, group.by ="celltype2", label = F,split.by = "condition",ncol=4, label.size = 5) + scale_color_npg()
plots3

ggsave2(filename = "./F1.UMAP_2.pdf", plots2, scale=0.77)
ggsave2(filename = "./F1.UMAP_3.pdf", plots3, scale=0.77)

#calculate cell percentage
cellcount <- as.data.frame(table(pbmc.integrated[[]]$subject))
pbmc.integrated$count <- cellcount$Freq[match(pbmc.integrated[[]]$subject,cellcount$Var1)]
table(pbmc.integrated$subject)
count1 <- group_by(pbmc.integrated@meta.data, subject, celltype) %>% 
  filter( ! celltype %in%  c("Megakaryocytes", "X25_CD8_NK","X26_B_UK","X27_UK")) %>%
  filter(! subject %in% c("S6", "S7")) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(count1)[2] <- "group"

sum(count1$freq)

metatable1 <- unique(select(pbmc.integrated[[]], subject, celltype, condition, chemistry))
count01 <- merge(count1,metatable1,by.x = c("subject","group"), by.y = c("subject", "celltype"))

ptable1 <- ggplot(count01, aes(x=condition, y=freq, color=condition)) +
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75)
ptable1 <- ptable1 + scale_y_continuous(labels = percent_format())
ptable1 <- ptable1 + facet_wrap(.~group, scales="free_y") + theme_cowplot(10) + scale_color_aaas()
#p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ptable1 <- ptable1 + stat_compare_means(label =  "p.format") #hide.ns=T,method = "t.test")
ptable1
ggsave2(filename = "./celltypebypercent.pdf",ptable1)

#perform cell count analysis based on known cell type, but no sub-cell type
count2 <- group_by(pbmc.integrated@meta.data, subject, celltype2) %>% 
  filter( ! celltype2 %in%  c("UK", "Megakaryocytes")) %>%
  filter(! subject %in% c("S6", "S7")) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(count2)[2] <- "group"

sum(count2$freq)

metatable2 <- unique(select(pbmc.integrated[[]], subject, celltype2, condition, chemistry))
count02 <- merge(count2,metatable2,by.x =  c("subject","group"), by.y = c("subject", "celltype2"))

ptable2 <- ggplot(count02, aes(x=condition, y=freq, color=condition)) + #color=chemistry
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75)
ptable2 <- ptable2 + scale_y_continuous(labels = percent_format())
ptable2 <- ptable2 + facet_wrap(.~group, scales="free_y") + theme_cowplot(10) + scale_color_aaas()
#p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ptable2 <- ptable2 + stat_compare_means(label =  "p.format")#, method = "t.test")
ptable2

ggsave2(filename = "./celltype2bypercent.pdf",ptable2)

#save(pbmc.integrated, file="./PBMC.inte.final.RData")

xtable <- table(pbmc.integrated$subject, pbmc.integrated$celltype)
# write.table(xtable, file = "~/Desktop/table.csv", sep = ",")

pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.6)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.5)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.4)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.3)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.2)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.1)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 1.0)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.9)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.8)

pUMAP.16 <- DimPlot(pbmc.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.1.6")
pUMAP.15 <- DimPlot(pbmc.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.1.5")
pUMAP.14 <- DimPlot(pbmc.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.1.4")
pUMAP.13 <- DimPlot(pbmc.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.1.3")
pUMAP.12 <- DimPlot(pbmc.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.1.2")
pUMAP.11 <- DimPlot(pbmc.integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.1.1")

(pUMAP.16+pUMAP.15)/(pUMAP.14+pUMAP.13)
#