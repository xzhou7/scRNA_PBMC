#Thrid round analysis Sep.21,2020
#PBMC analysis
#code clean for publications
#work under new parameters

#Author: Xin Zhou, Ph.D. 
#Last Update: Sep.21 2020

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

setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Monocytes/")
load("~/Desktop/PBMC related R objects/Monocytes.Final.RData")
#load("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Monocytes/PBMC.even.and.Integreted.RData")

#save(pbmc.even, file = "~/Desktop/PBMC.EVEN.RData")
#save(pbmc.integrated, file = "~/Desktop/PBMC.Clean.RData")

mono.nodc <- subset(pbmc.integrated, idents = c("X1_CD14Mono", "X16_CD16MONO", "X12_CD14Mono_InterM","X24_CD14Mono_PPBP"))

pahcarcode  <- colnames(mono.nodc)
subsettouse <- pahcarcode[pahcarcode!="ATACTTCGTTCTTCAT-8"]
mono.nodc <- subset(mono.nodc, cells = subsettouse)

DimPlot(mono.nodc)

table(mono.nodc$subject)

countmono <- group_by(mono.nodc@meta.data, subject, celltype) %>% 
  filter(! subject %in% c("S6","S7")) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(countmono)[2] <- "group"

sum(countmono$freq)

metatablemono <- unique(dplyr::select(mono.nodc[[]], subject, celltype, condition, chemistry))
count02 <- merge(countmono,metatablemono,by.x =  c("subject","group"), by.y = c("subject", "celltype"))

ptable2 <- ggplot(count02, aes(x=condition, y=freq, color=condition)) + #color=chemistry
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75) 
ptable2 <- ptable2 + scale_y_continuous(labels = percent_format())
ptable2 <- ptable2 + facet_wrap(.~group, scales="free_y") + theme_cowplot(10) + scale_color_npg()
#p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ptable2 <- ptable2 + stat_compare_means(label =  "p.format")#, method = "t.test")
ptable2

ggsave2(filename = "./MONObypercent.pdf",ptable2, width = 5, height = 5, dpi = 300)

count03 <- dcast(dplyr::select(count02,subject,group,condition,freq), subject + condition  ~ group, value.var= "freq")
ptable3 <- ggplot(count03, aes(x=X1_CD14Mono, y=X12_CD14Mono_InterM, color=condition)) + geom_point()
ptable3

count04 <- dcast(dplyr::select(count02,subject,group,condition,freq), subject + condition  ~ group, value.var= "freq")
ptable4 <- ggplot(count04, aes(x=condition, y=(X12_CD14Mono_InterM/X1_CD14Mono), color=condition)) + geom_point()
ptable4

countmono2 <- group_by(mono.nodc@meta.data, subject, celltype2) %>% 
  filter(! subject %in% c("S6", "S7")) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(countmono2)[2] <- "group"

sum(countmono2$freq)

metatablemono2 <- unique(dplyr::select(mono.nodc[[]], subject, celltype2, condition, chemistry))
count22 <- merge(countmono2,metatablemono2,by.x =  c("subject","group"), by.y = c("subject", "celltype2"))

ptable2 <- ggplot(count22, aes(x=condition, y=freq, color=condition)) + #color=chemistry
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75)
ptable2 <- ptable2 + scale_y_continuous(labels = percent_format())
ptable2 <- ptable2 + facet_wrap(.~group, scales="free_y") + theme_cowplot(10) + scale_color_npg()
#p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ptable2 <- ptable2 + stat_compare_means(label = "p.format")#, method = "t.test")
ptable2

ggsave2(filename = "./MONObypercent2.pdf", ptable2, width = 4, height = 3, dpi = 300)


DefaultAssay(mono.nodc) <- "integrated"

pvln <- VlnPlot(object = mono.nodc, features = c("CD14","FCGR3A", "HLA-DRA","PPBP"),cols=c("#00468BFF","#ED0000FF", "#42B540FF","#ADB6B6FF"),
                ncol=2, combine = T, pt.size=0.0, group.by = "celltype")
pvln


monocytes.marker <- c("CD14", "FCGR3A", "HLA-DRA", "PPBP")

for (i in monocytes.marker){
  print(i)
  p <- FeaturePlot(mono.nodc, features = i, pt.size = 0.01, order = T,col=c("grey","red"))
  ggsave(paste0("./CD14marker/", i, ".pdf"), p, width = 3, height = 2, dpi=300)
  }

#########################################Gate on Double Positive##############################################################
#assign doublepositive and doublenegative
mono.nodc$dp <- "1"
mono.nodc$dp[colnames(mono.nodc) %in% colnames(subset(mono.nodc, FCGR3A > 0 & CD14 > 0))] <- "doublepositive"
mono.nodc$dp[colnames(mono.nodc) %in% colnames(subset(mono.nodc, FCGR3A < 0 & CD14 < 0))] <- "doublenegative"
mono.nodc$dp[colnames(mono.nodc) %in% colnames(subset(mono.nodc, FCGR3A > 0 & CD14 < 0))] <- "CD16_only"
mono.nodc$dp[colnames(mono.nodc) %in% colnames(subset(mono.nodc, FCGR3A < 0 & CD14 > 0))] <- "CD14_only"

dndptable <- as.data.frame(table(mono.nodc$subject, mono.nodc$dp))
subjecttotal <- as.data.frame(table(mono.nodc$subject))
dndptable.meta <- merge(dndptable, subjecttotal, by = "Var1")
dndptable.meta$Freq <- dndptable.meta$Freq.x / dndptable.meta$Freq.y
dndptable.meta$condition <- "0"
dndptable.meta$condition[dndptable.meta$Var1 %in% c("S1", "S5", "S9", "S10", "S11", "S12")]<- "CONT"
dndptable.meta$condition[dndptable.meta$Var1 %in% c("S2","S3","S4","S6", "S7", "S8")]<- "PAH"

pdpdn <- ggplot(dndptable.meta, aes(x=condition, y = Freq)) + geom_point()
pdpdn <- pdpdn + facet_wrap(.~Var2)
pdpdn

table(mono.nodc$dp, mono.nodc$condition)
table( mono.nodc$condition)
monotable <- table(mono.nodc$dp, mono.nodc$condition)
monotable /c(3025,3025,3025,3025,2105 ,2105 ,2105 ,2105)

#plot CD16/CD14 scatter
pcd1416 <- FeatureScatter(mono.nodc, feature1 = "CD14", feature2 = "FCGR3A", pt.size = 0.5, 
                          slot="scale.data", cols=c("#42B540FF", "#00468BFF", "#AD002AFF","#0099B4FF"))
pcd1416 <- pcd1416 +  annotation_custom(tableGrob(round(table(mono.nodc$dp, mono.nodc$condition)/c(3025,3025,3025,3025,2105 ,2105 ,2105 ,2105 ),2)),
                                        xmin = 3,xmax = 8, ymin = 3, ymax=7)
pcd1416
#ggsave2(filename = "./CD14CD16scatter_in.pdf", pcd1416)

mono.dp <- subset(mono.nodc, FCGR3A > 0 & CD14 > 0)

FeatureScatter(mono.dp,feature1 = "CD14", feature2 = "FCGR3A") 
#################################################################################################################################
DefaultAssay(mono.nodc) <- "SCT"
table(mono.nodc$condition)

Idents(mono.nodc) <- "condition"
mono.control <- subset(mono.nodc, idents="CONT", downsample = 2104)
mono.even <- subset(mono.nodc, cells=c(colnames(nomo.control), colnames(subset(mono.nodc, idents="PAH"))))
mono.even$celltype <- factor(mono.even$celltype, levels = c("X1_CD14Mono","X12_CD14Mono_InterM","X24_CD14Mono_PPBP","X16_CD16MONO"))
table(mono.even$condition)

DefaultAssay(mono.even) <- "RNA"

p3 <- FeaturePlot(mono.even, features = c("TALDO1"), split.by = "condition", max.cutoff = 3, 
                  cols = c("grey", "red"), combine = F, order = T) 
#p3 <- lapply(X=p3, FUN = function(x) x + xlim(c(8, 15)) + ylim(c(-5, 8)))
p3 <- CombinePlots(p3, legend="right")
p3
ggsave2(filename = "./STAT1.pdf", p3, dpi=300)

mono.BMPR <- subset(mono.nodc, BMPR2 > 0)

table(mono.BMPR$celltype, mono.BMPR$subject)/table(mono.nodc$celltype, mono.nodc$subject)

p4 <- FeaturePlot(mono.even, features = c("PPBP", "S100A8", "HLA-DRA"), split.by = "condition", max.cutoff = 3, 
                  cols = c("grey", "red"), combine = F) 
#p4 <- lapply(X=p4, FUN = function(x) x + xlim(c(-15,-5)) + ylim(c(-5, 2.5)))
p4 <- CombinePlots(p4, legend="right")
p4
#ggsave2(filename = "./PPBPS100A8DRA.pdf", p4, dpi=300)

table(Idents(mono.even))

featureset1 <-  c("CD14", "FCGR3A", "HLA-DRA","CD36","STAT1","CXCR4")

#http://amigo.geneontology.org/amigo/term/GO:2000111#display-lineage-tab
featureset2 <- c("CDKN2A","SIRT1", "MEF2C", "CCR5", "CCL5", "BCL2", "SLC7A11", "PIK3CB")

plots <- VlnPlot(mono.even, features = "APAF1", split.by = "condition",  
                 pt.size = 0.01, combine = FALSE, split.plot = F) 
CombinePlots(plots = plots)

#DAG for different subpopulation
table(mono.it$celltype, mono.it$condition)
Idents(mono.it) <- "celltype"
DefaultAssay(mono.it) <- "integrated"
monomarker1.12 <- FindMarkers(mono.it, ident.1 = "X1_CD14Mono", ident.2 = "X12_CD14Mono_InterM", min.pct = 0.0, test.use="LR", logfc.threshold=0.3)
monomarker1.16 <- FindMarkers(mono.it, ident.1 = "X1_CD14Mono", ident.2 = "X16_CD16MONO", min.pct = 0.0, test.use="LR", logfc.threshold=0.3)
monomarker1.24 <- FindMarkers(mono.it, ident.1 = "X1_CD14Mono", ident.2 = "X24_CD14Mono_PPBP", min.pct = 0.0, test.use="LR", logfc.threshold=0.3)
monomarker12.16 <- FindMarkers(mono.it, ident.1 = "X12_CD14Mono_InterM", ident.2 = "X16_CD16MONO", min.pct = 0.0, test.use="LR", logfc.threshold=0.3)
monomarker12.24 <- FindMarkers(mono.it, ident.1 = "X12_CD14Mono_InterM", ident.2 = "X24_CD14Mono_3", min.pct = 0.0, test.use="LR",logfc.threshold=0.3)
monomarker16.24<- FindMarkers(mono.it, ident.1 = "X16_CD16MONO", ident.2 = "X24_CD14Mono_PPBP", min.pct = 0.0, test.use="LR",logfc.threshold=0.3)

write.csv(monomarker1.12, file = "./monomarker1.12.csv")
write.csv(monomarker1.16, file = "./monomarker1.16.csv")
write.csv(monomarker1.24, file = "./monomarker1.24.csv")
write.csv(monomarker12.16, file = "./monomarker12.16.csv")
write.csv(monomarker12.24, file = "./monomarker12.24.csv")
write.csv(monomarker16.24, file = "./monomarker16.24.csv")


cluster3.markers <- subset(monomarker1.24_LR, p_val<0.05)

cluster3.markers$genes <- row.names(cluster3.markers)

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on X1_X24 DEGs")
#p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.6 | avg_logFC >0.6 & log(p_val_adj) <25), aes(color="red"))
p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -100), aes(color="red"))
p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -100), aes(label=genes))
p


##############################################################################################################################
#DESEQ for PAH
##############################################################################################################################

Idents(mono.nodc) <- "condition"

Monocytes.DEG_MAST <- FindMarkers(mono.nodc,slot = "data", ident.1 = "PAH", ident.2 = "CONT", verbose = T,
                                  test.use = "DESeq2", assay = "SCT", min.pct = 0.02,logfc.threshold=0.3)

Monocytes.DEG_MAST$genes <- row.names(Monocytes.DEG_MAST)
cluster3.markers <- subset(Monocytes.DEG_MAST,  p_val_adj < 0.05)

#log(p_val_adj)
p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on Monocytes Cont/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0), aes(color="#4DBBD5FF"))
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < 0), aes(color="#E64B35FF"))
p <- p + geom_label_repel(data=subset(cluster3.markers, genes %in% X12.annotate), aes(label=genes))
#p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -20), aes(color="red"))
#p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -20), aes(label=genes))
p


mono.nodc$celltype.pah <- paste(mono.nodc$celltype, mono.nodc$condition, sep="_")
Idents(mono.nodc) <- "celltype.pah"

table(mono.nodc$celltype.pah)

#X1 DEG
DefaultAssay(mono.nodc) <- "SCT"
#DefaultAssay(mono.nodc) <- "RNA"
#mono.nodc <-  Seurat::NormalizeData(mono.it,normalization.method = "LogNormalize", scale.factor = 10000,verbose = TRUE,assay="RNA")
#mono.nodc <- ScaleData(mono.nodc, split.by = "condition")

#X1list  <- c(sample(x = colnames(subset(mono.it, idents="X1_CD14Mono_CONT")), size = 1000),
#                 sample(x = colnames(subset(mono.it, idents="X1_CD14Mono_PAH")), size = 1000))
#mono.X1 <- subset(mono.it, cells=X1list)

X1.annotate <- c("CD14", "MALAT1", "PLP2", "FOS", "ITGB2", "STAT1", "ANXA1", "JUNB", "HLA-B", "IFITM3", "DUSP1", "IRF1","PPDPF")
X12.annotate <- c("IFITM3", "CD14", "HLA-B", "PLP2", "FOS", "DUSP1", "JUNB", "STAT1", "ANXA1","PPDPF")
X16.annotate <- c("FCER1G", "PPDPF", "HLA-B", "STAT1", "HLA-DQB1", "IFITM3", "DUSP1", "JUNB")

totallist <- unique(c(X1.annotate, X12.annotate,X16.annotate))

DefaultAssay(mono.even) <- "integrated"

for (i in totallist){
  print(i)
  p3 <- FeaturePlot(mono.even, features = i, split.by = "condition", min.cutoff = 0,max.cutoff = 5, by.col=F, pt.size = 0.02,
                    cols = c("grey", "red"), order = T) 
  ggsave(paste0("./markers/MAX5", i, ".pdf"), p3, width=5, height = 4, dpi=300)
}

#table(mono.X1$celltype,mono.X1$condition)
#here we cannot use RNA but have to use SCT due to the merging of sequencing depth
mono.nodc
X1.DEG_DESEQ <- FindMarkers(mono.nodc,slot = "data", ident.2 = "X1_CD14Mono_CONT", ident.1 = "X1_CD14Mono_PAH", verbose = T,
                            test.use = "DESeq2", assay = "SCT", min.pct = 0.05,logfc.threshold=0.3, max.cells.per.ident = 1000)

X1.DEG_DESEQ$genes <- row.names(X1.DEG_DESEQ)
#write.csv(file = "./X1DEG.CD14.csv",X1.DEG_DESEQ)
cluster3.markers <- subset(X1.DEG_DESEQ, p_val_adj < 0.05)

#cluster3.markers[cluster3.markers$genes == "FTL",]
#VlnPlot(mono.it, features = "CXCR4", assay = "integrated",pt.size = 0, group.by = "celltype", split.by = "condition")

p <- ggplot(cluster3.markers, aes(x=avg_log2FC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on CD14 Monocyte Cont/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_log2FC > 0), aes(color="#4DBBD5FF"))
p <- p + geom_point(data=subset(cluster3.markers, avg_log2FC < 0), aes(color="#E64B35FF"))
p <- p + geom_label_repel(data=subset(cluster3.markers, genes %in% X1.annotate), aes(label=genes))
#p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -10), aes(label=genes))
p

ggsave(filename = "./deseq.X1.CD14_new.pdf", p, width=5, height = 6, dpi=300)
#write.csv(file = "./X1DEG-rev.csv", X1.DEG_DESEQ)

#X12 DEG
X12.DEG_DESEQ <- FindMarkers(mono.nodc, ident.2 = "X12_CD14Mono_InterM_CONT", ident.1 = "X12_CD14Mono_InterM_PAH", verbose = T,
                             test.use = "DESeq2",assay = "SCT",max.cells.per.ident=550, logfc.threshold=0.3, min.pct = 0.05)
X12.DEG_DESEQ$genes <- row.names(X12.DEG_DESEQ)
#write.csv(file = "./X12DEG.InterM.csv",X12.DEG_DESEQ)
cluster3.markers <- subset(X12.DEG_DESEQ,  p_val_adj < 0.05)

#log(p_val_adj)
p <- ggplot(cluster3.markers, aes(x=avg_log2FC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on X12_InterM Cont/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_log2FC > 0), aes(color="#4DBBD5FF"))
p <- p + geom_point(data=subset(cluster3.markers, avg_log2FC < 0), aes(color="#E64B35FF"))
p <- p + geom_label_repel(data=subset(cluster3.markers, genes %in% X12.annotate), aes(label=genes))
#p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -20), aes(color="red"))
#p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -20), aes(label=genes))
p

ggsave(filename = "./deseq.X12.Inter_new.pdf", p, width=5, height = 6, dpi=300)
write.csv(file = "./X12DEG-rev.csv", X12.DEG_DESEQ)

#X16 DEG
X16.DEG_DESEQ <- FindMarkers(mono.nodc, ident.2 = "X16_CD16MONO_CONT", ident.1 = "X16_CD16MONO_PAH", verbose = T,
                             test.use = "DESeq2",assay = "SCT",max.cells.per.ident=390, logfc.threshold=0.3, min.pct = 0.05)
X16.DEG_DESEQ$genes <- row.names(X16.DEG_DESEQ)
#write.csv(file = "./X16DEG.csv",X16.DEG_DESEQ)
cluster3.markers <- subset(X16.DEG_DESEQ, p_val_adj < 0.05)

p <- ggplot(cluster3.markers, aes(x=avg_log2FC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on X16_CD16 Mono Cont/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_log2FC > 0), aes(color="#4DBBD5FF"))
p <- p + geom_point(data=subset(cluster3.markers, avg_log2FC < 0), aes(color="#E64B35FF"))
p <- p + geom_label_repel(data=subset(cluster3.markers, genes %in% X16.annotate), aes(label=genes))
#p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -10), aes(color="red"))
#p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -10), aes(label=genes))
p

ggsave(filename = "./deseq.X16.Inter_new.pdf", p, width=5, height = 6, dpi=300)
write.csv(file = "./X16DEG-rev.csv", X16.DEG_DESEQ)

#X24 DEG
X24.DEG_DESEQ <- FindMarkers(mono.nodc, ident.2 = "X24_CD14Mono_PPBP_CONT", ident.1 = "X24_CD14Mono_PPBP_PAH", verbose = T,
                             test.use = "DESeq2",assay = "SCT",max.cells.per.ident=40, logfc.threshold=0.3, min.pct = 0.05)

X24.DEG_DESEQ$genes <- row.names(X24.DEG_DESEQ)
#write.csv(file = "./X24DEG.csv",X24.DEG_DESEQ)

cluster3.markers <- subset(X24.DEG_DESEQ, p_val_adj < 0.05)

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on X24_PPBP Mono Cont/PAH DEGs")
#p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -10), aes(color="red"))
#p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -10), aes(label=genes))
p

#cdc DEG
cDC.DEG_DESEQ <- FindMarkers(mono.nodc, ident.1 = "X19_cDC_CONT", ident.2 = "X19_cDC_PAH", verbose = T,
                             test.use = "DESeq2",assay = "SCT",max.cells.per.ident=90, logfc.threshold=0.2, min.pct = 0.05)
cDC.DEG_DESEQ$genes <- row.names(cDC.DEG_DESEQ)
#write.csv(file = "./cdcDEG.csv",cDC.DEG_DESEQ)
cluster3.markers <- subset(cDC.DEG_DESEQ, p_val_adj < 0.05)

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on cDC Cont/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -0.05), aes(color="red"))
p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -0.05), aes(label=genes))
p

X1list <- X1.DEG_DESEQ[X1.DEG_DESEQ$p_val_adj < 1e-7,]
X1list <- X1list[order(X1list$avg_logFC),"genes"]

X12list <- X12.DEG_DESEQ[X12.DEG_DESEQ$p_val_adj< 1e-7,]
X12list <- X12list[order(X12list$avg_logFC),"genes"]

X16list <- X16.DEG_DESEQ[X16.DEG_DESEQ$p_val_adj< 1e-7,"genes"]

DefaultAssay(mono.nodc) <- "RNA"
mono.nodc <- ScaleData(mono.nodc)
Idents(mono.nodc) <- "celltype"

heatobject1 <- subset(mono.nodc, idents="X1_CD14Mono")
p1 <- DoHeatmap(heatobject1, features = X1list, group.by= "condition")
p1

heatobject2 <- subset(mono.it, idents="X12_CD14Mono_InterM")
p12 <- DoHeatmap(heatobject2, features = X12list, group.by= "condition")
p12

listFC <- data.frame()
listFC.X1 <- X1.DEG_DESEQ[X1.DEG_DESEQ$p_val_adj < 0.05,]
listFC.X12 <- X12.DEG_DESEQ[X12.DEG_DESEQ$p_val_adj < 0.05,]

listFC <-  merge(listFC.X1, listFC.X12, by = "row.names", all.x=T, all.y=T)

p <- ggplot(listFC, aes(x=avg_logFC.x, y=avg_logFC.y)) + geom_point()
p <- p + geom_text_repel(aes(x=avg_logFC.x, y=avg_logFC.y,label=Row.names))
p

############################################################################################################
#cell cycle analysis
############################################################################################################
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

DimPlot(mono.nodc,split.by = "celltype")
DefaultAssay(mono.nodc) <- "integrated"

mono.nodc<- CellCycleScoring(mono.nodc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
p <- ggplot(mono.nodc[[]], aes(x=S.Score, y=G2M.Score, color=Phase)) + geom_point(size=0.2)
p <- p + facet_grid(condition~celltype)
p

table(mono.nodc$Phase, mono.nodc$condition)
table(mono.nodc$Phase, mono.nodc$celltype)

DimPlot(mono.nodc)

##################################################################
#pathway related analysis
##################################################################
uplist <- filter(X1.DEG_DESEQ, avg_logFC>0 & p_val_adj < 0.05)
downlist <- filter(X1.DEG_DESEQ, avg_logFC<0 & p_val_adj < 0.05)

dim(uplist)
dim(downlist)

upgeneList <- uplist[, 6]
names(upgeneList) <- as.character(uplist$genes)
upgeneList <- sort(upgeneList, decreasing = TRUE)
head(upgeneList)

downlistList <- downlist[, 6]
names(downlistList) <- as.character(downlist$genes)
downlistList <- sort(downlistList, decreasing = TRUE)
head(downlistList)

length(upgeneList)
length(downlistList)

upgene.df <- bitr(names(upgeneList), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Hs.eg.db)
downgene.df <- bitr(names(downlistList), fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"),OrgDb = org.Hs.eg.db)

gene <- upgene.df$ENTREZID

################################################################
#For gse analysis
################################################################
#cd14
geneList0 <- filter(Monocytes.DEG_MAST, avg_logFC < 0)
geneList1 <- geneList0[,2]
names(geneList1) <- as.character(geneList0$genes)
geneList1 <- sort(geneList1, decreasing = TRUE)

geneList2 <- bitr(names(geneList1), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)

names(geneList1) <- geneList2$ENTREZID[match(names(geneList1),geneList2$SYMBOL)]
geneList <- subset(geneList1, !is.na(names(geneList1)))

#Disease Ontology
DO14 <- gseDO(geneList, 
              nPerm         = 100000,
              pvalueCutoff  = 0.05,
              minGSSize    = 10,
              pAdjustMethod = "BH",
              seed = T,
              verbose       = T)

DO14 <- setReadable(DO14, 'org.Hs.eg.db')
DO14@result
cnetplot(DO14, foldChange=geneList, showCategory = 15)

write.csv(DO14@result, file = "./DiseaseOntology_cd14Mono.csv")

cnetplot(DO14, foldChange=geneList) 
heatplot(DO14, foldChange=geneList)
emapplot(DO14)

#GO
GO14 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              pvalueCutoff = 0.1,
              verbose      = T)
GO14 <- setReadable(GO14, 'org.Hs.eg.db')
GO14@result
cnetplot(GO14, foldChange=geneList,showCategory = 20) 

#KEGG
KEGG14 <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 3,
                  pvalueCutoff = 0.2,
                  use_internal_data=F,
                  verbose      = T)
KEGG14 <- setReadable(KEGG14, 'org.Hs.eg.db',keyType = "ENTREZID")
KEGG14@result
write.csv(KEGG14@result, file="./CD14mono_KEGG.csv")
cnetplot(KEGG14, foldChange=geneList) 

#RE
RE14 <- gsePathway(geneList  = geneList,
                   organism     = 'human',
                   nPerm        = 2000,
                   minGSSize    = 5,
                   pvalueCutoff = 0.1,
                   verbose      = T)
RE14 <- setReadable(RE14, 'org.Hs.eg.db',keyType = "ENTREZID")
RE14@result
cnetplot(RE14, foldChange=geneList) 
heatplot(RE14,foldChange=geneList)


###########
# Intermediate
###########
geneList0 <- filter(X12.DEG_DESEQ, p_val_adj < 0.05, abs(avg_logFC) > 0.1)
geneList1 <- geneList0[,2]
names(geneList1) <- as.character(geneList0$genes)
geneList1 <- sort(geneList1, decreasing = TRUE)

geneList2 <- bitr(names(geneList1), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)

names(geneList1) <- geneList2$ENTREZID[match(names(geneList1),geneList2$SYMBOL)]
geneList <- subset(geneList1, !is.na(names(geneList1)))

#DO
DO12 <- gseDO(geneList,
              nPerm         = 1000,
              minGSSize     = 5,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = T)

DO12 <- setReadable(DO12, 'org.Hs.eg.db')
cnetplot(DO12, foldChange=geneList) 

write.csv(DO12@result, file = "./DO_cd12Mono.csv")

cnetplot(DO12, foldChange=geneList) 
heatplot(DO12, foldChange=geneList)
emapplot(DO12)

#GO
GO12 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 10,
              pvalueCutoff = 0.2,
              by           = "DOSE",
              verbose      = T)

GO12 <- setReadable(GO12, 'org.Hs.eg.db')
GO12@result
cnetplot(GO12, foldChange=geneList) 

#KEGG
KEGG12 <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 3,
                  pvalueCutoff = 0.2,
                  use_internal_data=F,
                  verbose      = T)
KEGG12 <- setReadable(KEGG12, 'org.Hs.eg.db',keyType = "ENTREZID")
KEGG12@result
write.csv(KEGG12@result, file="./CD14CD16mono_KEGG.csv")
cnetplot(KEGG12, foldChange=geneList) 

#RE
RE12 <- gsePathway(geneList  = geneList,
                   organism     = 'human',
                   nPerm        = 2000,
                   minGSSize    = 5,
                   pvalueCutoff = 0.1,
                   verbose      = T)
RE12 <- setReadable(RE12, 'org.Hs.eg.db',keyType = "ENTREZID")
RE12@result
cnetplot(RE12, foldChange=geneList) 
heatplot(RE12,foldChange=geneList)

###########
# CD16
###########
geneList0 <- filter(X16.DEG_DESEQ, p_val_adj < 0.05, abs(avg_logFC) > 0.1)
geneList1 <- geneList0[,2]
names(geneList1) <- as.character(geneList0$genes)
geneList1 <- sort(geneList1, decreasing = TRUE)

geneList2 <- bitr(names(geneList1), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)

names(geneList1) <- geneList2$ENTREZID[match(names(geneList1),geneList2$SYMBOL)]
geneList <- subset(geneList1, !is.na(names(geneList1)))

#DO
DO16 <- gseDO(geneList,
              nPerm         = 1000,
              minGSSize     = 5,
              pvalueCutoff  = 0.2,
              pAdjustMethod = "BH",
              verbose       = T)

DO16 <- setReadable(DO16, 'org.Hs.eg.db')
cnetplot(DO16, foldChange=geneList) 

#write.csv(DO16@result, file = "./DO_cd16Mono.csv")

cnetplot(DO16, foldChange=geneList) 
heatplot(DO16, foldChange=geneList)
emapplot(DO16)

#GO
GO16 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 10,
              pvalueCutoff = 0.2,
              by           = "DOSE",
              verbose      = T)

GO16 <- setReadable(GO16, 'org.Hs.eg.db')
GO16@result
cnetplot(GO16, foldChange=geneList) 

#KEGG
KEGG16 <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 3,
                  pvalueCutoff = 0.2,
                  use_internal_data=F,
                  verbose      = T)
KEGG16 <- setReadable(KEGG16, 'org.Hs.eg.db',keyType = "ENTREZID")
KEGG16@result
write.csv(KEGG16@result, file="./CD16mono_KEGG.csv")
cnetplot(KEGG16, foldChange=geneList) 

#RE
RE16 <- gsePathway(geneList  = geneList,
                   organism     = 'human',
                   nPerm        = 2000,
                   minGSSize    = 5,
                   pvalueCutoff = 0.1,
                   verbose      = T)
RE16 <- setReadable(RE16, 'org.Hs.eg.db',keyType = "ENTREZID")
RE16@result
cnetplot(RE12, foldChange=geneList) 
heatplot(RE12,foldChange=geneList)


#####################################
#total CD14 Mono
#####################################
table(Idents(pbmc.integrated))
table((pbmc.integrated$celltype2))
CD14MONO <- subset(pbmc.integrated, idents = c("X1_CD14Mono", "X12_CD14Mono_InterM"))

table(CD14MONO$celltype)
table(CD14MONO$condition)
Idents(CD14MONO) <- "condition"
Idents(CD14MONO)

cd14markers <-FindMarkers(CD14MONO, ident.1 = "CONT", ident.2 = "PAH", verbose = T,
                          test.use = "DESeq2",assay = "SCT",max.cells.per.ident=1600, logfc.threshold=0.3, min.pct = 0.2)
cd14markers$genes <- row.names(cd14markers)

cluster3.markers <- subset(cd14markers, p_val_adj < 0.05)

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on CD14 Mono Cont/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, log(p_val_adj) < -10), aes(color="red"))
p <- p + geom_label_repel(data=subset(cluster3.markers,  log(p_val_adj) < -10), aes(label=genes))
p

#cd14
monomarker1.24_LR$genes <- row.names(monomarker1.24_LR)
geneList0 <- filter(cd14markers, p_val_adj < 0.00000005, abs(avg_logFC) > 0.1)
geneList1 <- geneList0[,2]
names(geneList1) <- as.character(geneList0$genes)
geneList1 <- sort(geneList1, decreasing = TRUE)

geneList2 <- bitr(names(geneList1), fromType = "SYMBOL", toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)

names(geneList1) <- geneList2$ENTREZID[match(names(geneList1),geneList2$SYMBOL)]
geneList <- subset(geneList1, !is.na(names(geneList1)))

#DO
DO14 <- gseDO(geneList,
              nPerm         = 1000,
              minGSSize     = 5,
              pvalueCutoff  = 0.4,
              pAdjustMethod = "BH",
              verbose       = T)

DO14 <- setReadable(DO14, 'org.Hs.eg.db')
DO14@result
cnetplot(DO14, foldChange=geneList) 

write.csv(DO14@result, file = "./DO_cd14Mono.csv")

cnetplot(DO14, foldChange=geneList) 
heatplot(DO14, foldChange=geneList)
emapplot(DO14)

#GO
GO14 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              nPerm        = 200,
              minGSSize    = 2,
              pvalueCutoff = 0.4,
              verbose      = T)
GO14 <- setReadable(GO14, 'org.Hs.eg.db')
GO14@result
cnetplot(GO14, foldChange=geneList) 

#KEGG
KEGG14 <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 500,
                  minGSSize    = 3,
                  pvalueCutoff = 0.4,
                  use_internal_data=F,
                  verbose      = T)
KEGG14 <- setReadable(KEGG14, 'org.Hs.eg.db',keyType = "ENTREZID")
KEGG14@result
write.csv(KEGG14@result, file="./CD14mono_KEGG.csv")
cnetplot(KEGG14, foldChange=geneList) 

#RE
RE14 <- gsePathway(geneList  = geneList,
                   organism     = 'human',
                   nPerm        = 500,
                   minGSSize    = 3,
                   pvalueCutoff = 0.2,
                   verbose      = T)
RE14 <- setReadable(RE14, 'org.Hs.eg.db',keyType = "ENTREZID")
RE14@result
cnetplot(RE14, foldChange=geneList) 
heatplot(RE14,foldChange=geneList)


upsetplot(y)
emapplot(y,foldChange=geneList)
heatplot(y,foldChange=geneList)
cnetplot(RE14, node_label="category",foldChange=geneList) 

dotplot(y,foldChange=geneList)
upsetplot(rc)

gseaplot2(rc, geneSetID = 1:3)
gseaplot(rc, geneSetID = 1, by = "runningScore", title = rc$Description[1])

#save(pbmc.integrated, mono.it, file = "./after 02 analysis/pbmc_and_mono.RData")



######################################################################
#GO analysis X1
######################################################################

#generate global background for the first time
geneList_glo <-bitr(rownames(mono.nodc@assays$RNA), fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

geneList0_u <- filter(X1.DEG_DESEQ, p_val_adj < 0.05, avg_logFC < (-0.1))
geneList0_d <- filter(X1.DEG_DESEQ, p_val_adj < 0.05, avg_logFC > 0.1)

#geneList_X <-bitr(X1.DEG_DESEQ$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

dim(geneList0_u)
dim(geneList0_d)

geneList2_u <- bitr(geneList0_u$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
geneList2_d <- bitr(geneList0_d$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

ego_1u <- enrichGO(gene       = (geneList2_u$ENTREZID),
                   universe      = (geneList_glo$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego_1u <- setReadable(ego_1u, 'org.Hs.eg.db',keyType = "ENTREZID")
head(summary(ego_1u))
GO1u <- ego_1u@result
up1 <- dotplot(ego_1u, showCategory=10)
up1
ggsave(filename = "./GO_X1_DowninPAH.pdf",up1, width = 9, height = 3, dpi=300)
up1_2 <- cnetplot(ego_1u, node_label="all",showCategory=10)
up1_2
ggsave(filename = "./GO_X1_DowninPAH_cnet.pdf",up1_2, width = 9, height = 6, dpi=300)
write.csv(file="./GOX1_DowninPAH.csv",GO1u)
#cnetplot(ego_1u, foldChange=geneList0_u$avg_logFC,showCategory=10)

ego_1d <- enrichGO(gene          = (geneList2_d$ENTREZID),
                   universe      = (geneList_glo$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego_1d <- setReadable(ego_1d, 'org.Hs.eg.db',keyType = "ENTREZID")
head(summary(ego_1d))
GO1d <- ego_1d@result
Dn1 <- dotplot(ego_1d, showCategory=10)
Dn1
ggsave(filename = "./GO_X1_UpinPAH.pdf",Dn1, width = 9, height = 3, dpi=300)
Dn1_2 <- cnetplot(ego_1d, node_label="all",showCategory=10)
Dn1_2
ggsave(filename = "./GO_X1_UpinPAH_cnet.pdf",Dn1_2, width = 9, height = 6, dpi=300)
write.csv(file="./GOX1_UpinPAH.csv",GO1d)

######################################################################
#GO analysis X12
######################################################################
geneList0_u <- filter(X12.DEG_DESEQ, p_val_adj < 0.05, avg_logFC < (-0.1))
geneList0_d <- filter(X12.DEG_DESEQ, p_val_adj < 0.05, avg_logFC > 0.1)

#geneList_X <-bitr(X1.DEG_DESEQ$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

dim(geneList0_u)
dim(geneList0_d)

geneList2_u <- bitr(geneList0_u$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
geneList2_d <- bitr(geneList0_d$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

ego_12u <- enrichGO(gene       = (geneList2_u$ENTREZID),
                   universe      = (geneList_glo$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego_12u <- setReadable(ego_12u, 'org.Hs.eg.db', keyType = "ENTREZID")
head(summary(ego_12u))
GO12u <- ego_12u@result
up12 <- dotplot(ego_12u, showCategory=20)
up12
ggsave(filename = "./GO_X12_DowninPAH.pdf",up12, width = 9, height = 3, dpi=300)
up12_2 <- cnetplot(ego_12u, node_label="all",showCategory=10)
up12_2
ggsave(filename = "./GO_X12_DowninPAH_cnet.pdf",up12_2, width = 9, height = 6, dpi=300)
write.csv(file="./GOX12_DowninPAH.csv",GO12u)
#cnetplot(ego_12u, foldChange=geneList0_u$avg_logFC,showCategory=10)

ego_12d <- enrichGO(gene          = (geneList2_d$ENTREZID),
                   universe      = (geneList_glo$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego_12d <- setReadable(ego_12d, 'org.Hs.eg.db',keyType = "ENTREZID")
head(summary(ego_12d))
GO12d <- ego_12d@result
Dn12 <- dotplot(ego_12d, showCategory=10)
Dn12
ggsave(filename = "./GO_X12_UpinPAH.pdf",Dn12, width = 10, height = 3, dpi=300)
Dn12_2 <- cnetplot(ego_12d, node_label="all",showCategory=10)
Dn12_2
ggsave(filename = "./GO_X12_UpinPAH_cnet.pdf",Dn12_2, width = 9, height = 6, dpi=300)
write.csv(file="./GOX12_UpinPAH.csv",GO12d)

######################################################################
#GO analysis X16
######################################################################

geneList0_u <- filter(X16.DEG_DESEQ, p_val_adj < 0.05, avg_logFC < (-0.1))
geneList0_d <- filter(X16.DEG_DESEQ, p_val_adj < 0.05, avg_logFC > 0.1)

dim(geneList0_u)
dim(geneList0_d)

geneList2_u <- bitr(geneList0_u$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
geneList2_d <- bitr(geneList0_d$genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

ego_16u <- enrichGO(gene       = (geneList2_u$ENTREZID),
                   universe      = (geneList_glo$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego_16u <- setReadable(ego_16u, 'org.Hs.eg.db',keyType = "ENTREZID")
head(summary(ego_16u))
GO16u <- ego_16u@result
up16 <- dotplot(ego_16u, showCategory=10)
up16
ggsave(filename = "./GO_X16_DowninPAH.pdf",up16, width = 9, height = 3, dpi=300)
up16_2 <- cnetplot(ego_16u, node_label="all",showCategory=10)
up16_2
ggsave(filename = "./GO_X16_DowninPAH_cnet.pdf",up16_2, width = 9, height = 6, dpi=300)
write.csv(file="./GOX16_DowninPAH.csv",GO16u)
#cnetplot(ego_1u, foldChange=geneList0_u$avg_logFC,showCategory=10)

ego_16d <- enrichGO(gene          = (geneList2_d$ENTREZID),
                   universe      = (geneList_glo$ENTREZID),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.1,
                   qvalueCutoff  = 0.1)

ego_16d <- setReadable(ego_16d, 'org.Hs.eg.db',keyType = "ENTREZID")
head(summary(ego_16d))
GO16d <- ego_16d@result
Dn16 <- dotplot(ego_16d, showCategory=10)
Dn16
ggsave(filename = "./GO_X16_UpinPAH.pdf",Dn16, width = 9, height = 3, dpi=300)
Dn16_2 <- cnetplot(ego_16d, node_label="all",showCategory=10)
Dn16_2
ggsave(filename = "./GO_X16_UpinPAH_cnet.pdf",Dn16_2, width = 9, height = 6, dpi=300)
write.csv(file="./GOX16_UpinPAH.csv",GO16d)

save(mono.nodc, file = "~/Desktop/Monocytes.Final.RData")

#############
#T and NK
############

T.NK <- subset(pbmc.integrated, idents = c("X6_CD4_5", "X13_CD4_FOXP3",
                                           "X4_NaiveCD8","X10_CD8_1",
                                           "X2_CD4_2","X17_CD8_RORC",
                                           "X3_CD4_3","X9_NK1",
                                           "X20_NK2","X0_CD4_1",
                                           "X5_CD4_4", "X15_NK2",
                                           "X14_CD4CD8", "X7_D8_NKT?",
                                           "X18_CD4_UNK","X23_CD4_UNK2",
                                           "X25_CD8_UK"))

mono.it <- subset(pbmc.integrated, idents = c("X1_CD14Mono", "X16_CD16MONO", "X12_CD14Mono_InterM", "X19_cDC","X24_CD14Mono_PPBP"))
mono.it$celltype <- factor(mono.it$celltype, levels = c("X1_CD14Mono", "X12_CD14Mono_InterM","X16_CD16MONO", "X24_CD14Mono_PPBP","X19_cDC"))
mono.it
Idents(mono.it) <- "celltype"

save(T.NK, file = "../T cell/t.nk.RData")


