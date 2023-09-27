#Second round analysis Apr.03,2020
#PBMC analysis
#Set new cutoffs and compare parameters
#work under new parameters

#Author: Xin Zhou, Ph.D. 
#Last Update: March.28 2023

library(Seurat)
library(dplyr)
library(reshape2)
library(Matrix)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(ggsci)
library(patchwork)
library(scales)
library(DESeq2)
library(MAST)
library(gridExtra)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(robustSingleCell)
library(ReactomePA)

setwd("~/Library/CloudStorage/Box-Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/")

load("./PBMC.inte.final.RData")

#put nnormalized data in RNA and SCT slot
#note: normalizaton works differently for RNA assay and SCT assay, RNA assay has a scale factor, SCT assay doesn't work this way (direct log)

DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, normalization.method = "LogNormalize",scale.factor = 10000)
pbmc.integrated <- ScaleData(pbmc.integrated)

DefaultAssay(pbmc.integrated) <- "SCT"
pbmc.integrated <- NormalizeData(pbmc.integrated, normalization.method = "LogNormalize")
pbmc.integrated <- ScaleData(pbmc.integrated)

Idents(pbmc.integrated) <- "celltype2"
pbmc.integrated <- subset(pbmc.integrated, idents=c("CD4_T","CD8_T","cDC","B","CD14Mono","CD16Mono","pDC","NK"))

Idents(pbmc.integrated) <- "celltype"
DimPlot(pbmc.integrated, label = T)

table(pbmc.integrated$celltype2)

umap <- as.data.frame(pbmc.integrated@reductions[["umap"]]@cell.embeddings)
umap

rownames(umap[umap$UMAP_2< (-8),])

pahcarcode  <- colnames(pbmc.integrated)
subsettouse <- pahcarcode[pahcarcode!="ATCAGGTGTCTTGAGT-7"&
                            pahcarcode!="GGGACAAGTCGTTGCG-7"&
                            pahcarcode!="CGATGGCGTTATCTTC-8"]

pbmc.integrated <- subset(pbmc.integrated, cells = subsettouse)

pbmc.integrated$celltype3 <- pbmc.integrated$celltype2
pbmc.integrated$celltype3[pbmc.integrated$celltype3=="CD14Mono"] <- "Monocytes"
pbmc.integrated$celltype3[pbmc.integrated$celltype3=="CD16Mono"] <- "Monocytes"

pbmc.integrated$celltype3 <- factor(pbmc.integrated$celltype3, levels = c("CD4_T","CD8_T","cDC","B",
                                                                          "Monocytes","pDC","NK"))

Idents(pbmc.integrated) <- "celltype3"

pumap3 <- DimPlot(pbmc.integrated, label = T, pt.size = 0.005) + scale_color_d3()
pumap3

#ggsave2(filename = "./UMAP3.whole.pdf",pumap3, width = 6, height = 5, dpi=300)

table(Idents(pbmc.integrated))

Idents(pbmc.integrated) <- "celltype"
table(pbmc.integrated$celltype)

pumap1 <- DimPlot(pbmc.integrated, label = T, pt.size = 0.005) + NoLegend() + coord_fixed()
pumap1
ggsave2("./UMAP1.pdf", pumap1, width = 6, height = 5, dpi=300)

markers.to.plot <- c("CD3D","CD4","GPR183","FOXP3", "CD14","S100A9","FCGR3A", "HLA-DRA", "SELL","CD8A", "GNLY", "NKG7", "CCL5", 
                      "MS4A1", "CD79A", "VMO1", "S100A4", "HLA-DQA1","TSPAN13", "IL3RA", "PPBP", "GNG11")

DefaultAssay(pbmc.integrated) <- "integrated"

#subset(pbmc.integrated, idents = c( "X25_CD8_UK", "X26_B_UK", "X18_CD4_UNK"), invert = T) %>%

pbmc.integrated$celltype_order <- pbmc.integrated$celltype

pbmc.integrated$celltype_order <- factor(pbmc.integrated$celltype_order, levels = c("X0_CD4_1" , "X13_CD4_FOXP3", "X2_CD4_2", "X3_CD4_3", "X5_CD4_4", "X6_CD4_5",     
                                                                                           "X1_CD14Mono", "X12_CD14Mono_InterM","X16_CD16MONO"  ,"X19_cDC" ,  "X24_CD14Mono_PPBP",  
                                                                                           "X10_CD8_1", "X14_CD4CD8", "X17_CD8_RORC", "X4_NaiveCD8",             
                                                                                           "X11_B2", "X8_B1", 
                                                                                           "X15_NK2", "X20_NK2", "X7_D8_NKT?","X9_NK1", 
                                                                                           "X22_pDC",   
                                                                                           "X18_CD4_UNK", "X23_CD4_UNK2","X25_CD8_UK", "X26_B_UK"))
Idents(pbmc.integrated) <- "celltype_order"

pbmc.clean <- subset(pbmc.integrated, idents = c("X18_CD4_UNK", "X23_CD4_UNK2","X25_CD8_UK", "X26_B_UK"), invert = TRUE)
table(pbmc.clean$condition)

Dot_marker <- DotPlot(pbmc.clean,features = markers.to.plot) & theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Dot_marker

ggsave2(filename = "./Dot.Markers.pdf",Dot_marker, width = 8, height = 4.5, dpi = 300)

table(pbmc.integrated$condition)

#######################################################
table(pbmc.integrated$condition)

Idents(pbmc.integrated) <- "condition"

pbmc.integrated.control <- subset(pbmc.integrated, idents="CONT", downsample = 14030)

pbmc.even <- subset(pbmc.clean, cells=c(colnames(pbmc.integrated.control), colnames(subset(pbmc.integrated, idents="PAH"))))
pbmc.even$celltype3 <- factor(pbmc.even$celltype3, levels = c("CD4_T","CD8_T","cDC","B","Monocytes","pDC","NK"))
table(pbmc.even$condition)

pcp <- DimPlot(pbmc.even) + scale_color_npg() + coord_fixed()
pcp
#ggsave2(filename = "./UMAP4.cvsp.pdf",pcp, width = 6, height = 5, dpi=300)

p.percent <- table(pbmc.even$celltype_order, pbmc.even$condition) %>% data.frame() %>% filter(Freq != 0) %>% ggplot(aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(position="dodge", stat="identity")
p.percent <- p.percent +  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p.percent
ggsave("./celltype.by.conditions.pdf", p.percent, width = 7, height = 4, dpi = 300)

pumap.cp <- DimPlot(pbmc.even, label = T) #+ scale_color_d3() , group.by = "condition",split.by = "condition"
pumap.cp
#ggsave2(filename = "./UMAP5.cvsp.pdf",pumap.cp, width = 10, height = 5, dpi=300)

FeaturePlot(pbmc.even, features = "RORC", order = T, split.by = "condition")
RidgePlot(pbmc.even, features = "RORC")
VlnPlot(pbmc.even, features = "NEAT1", group.by = "")

DefaultAssay(pbmc.even)
Markers <- c("CD14", "CD3E","CD19","CD1C")
for (i in Markers){
  print(i)
  pi <- FeaturePlot(pbmc.even, features = i, order = T)
  ggsave2(filename = paste0("./UMAP.markers/",i, ".pdf"),pi, width = 3, height =2.5, dpi=300)
}

#######################################################

############################### (30minutes + )
DefaultAssay(pbmc.integrated) <- "integrated"
Idents(pbmc.integrated) <- "celltype"
MarkerforALL <- FindAllMarkers(pbmc.integrated,logfc.threshold = 0.25,test.use = "LR",only.pos = T)
write.csv(file="./All_Markers_by_Cluster_New.csv", MarkerforALL)
##############################

a <- as.data.frame(table(pbmc.integrated$subject,pbmc.integrated$celltype))
b <- as.data.frame(table(pbmc.integrated$subject, pbmc.integrated$condition)) %>% subset(Freq!=0)


colnames(a) <- c("Var1", "Celltype", "count") 

c <- merge(a, b, by= "Var1")
c$Percentage <- c$count/c$Freq

c$Percentage_Log <- log(c$Percentage)
c$Percentage_Log[c$Percentage_Log==Inf] <- "0"
c$Percentage_Log<- as.numeric(c$Percentage_Log) 

c_clean <- subset(c,  Var1 != "S6" & Var1 != "S7")

table(c$Var1)
table(c_clean$Var1)

p <- ggplot(c_clean, aes(x=Var2, y=Percentage,color=Var2, fill=Var2)) + geom_dotplot(binaxis='y', stackdir='center')
p <- p + facet_wrap(.~Celltype, scales = "free_y") + scale_color_npg() + theme_cowplot()+ scale_fill_npg()
p

p.indivi <- ggplot(c_clean, aes(x=Var2, y=Percentage,color=Var1)) + geom_dotplot(binaxis='y', stackdir='center')
p.indivi <- p.indivi + facet_wrap(.~Celltype, scales = "free_y")
p.indivi

p1 <- p + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5,alpha=0.5)
p1

p2 <- p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="grey", width=0.3,alpha=0.5) + stat_summary(fun.y=mean, geom="point", color="black", size=2, alpha=0.7)
p2

p3 <- p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="grey", alpha=0.5)
p3

subset(c, Celltype== "X26_B_UK")
subset(c, Celltype== "X25_CD8_UK")
subset(c, Celltype== "X24_CD14Mono_PPBP")
subset(c, Celltype== "X23_CD4_UNK2")

subset(c_clean, Celltype== "X15_NK2")
subset(c_clean, Celltype== "X20_NK2")

subset(c_clean, Celltype== "X17_CD8_RORC")

table(c_clean$Celltype)
target <- "X0_CD4_1"
t.testresult <- t.test(filter(c_clean, Var2  == "CONT" & Celltype==target)$Percentage,filter(c_clean, Var2  == "PAH" & Celltype==target)$Percentage)
t.testresult

log(t.testresult$estimate[2]/t.testresult$estimate[1])

wilcox.test(filter(c_clean, Var2  == "CONT" & Celltype==target)$Percentage,filter(c_clean, Var2  == "PAH" & Celltype==target)$Percentage)

save(pbmc.even, pbmc.integrated,MarkerforALL, file = "~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Monocytes/PBMC.even.and.Integreted.RData")

