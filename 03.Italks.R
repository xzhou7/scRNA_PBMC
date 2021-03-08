#PBMC with ligand receptor analysis

library(iTALK)
library(ggplot2)
library(Seurat)
library(patchwork)
library(dplyr)
library(monocle3)

#italk: https://www.biorxiv.org/content/10.1101/507871v1.full.pdf+html

#Set working directory
setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/italksPBMC+DC/")
load(file = "../after 02 analysis/pbmc_and_mono.RData")

#put essential normalization
DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated, normalization.method = "LogNormalize",scale.factor = 10000)
pbmc.integrated <- ScaleData(pbmc.integrated)

#DefaultAssay(pbmc.integrated) <- "SCT"
#pbmc.integrated <- NormalizeData(pbmc.integrated, normalization.method = "LogNormalize")
#pbmc.integrated <- ScaleData(pbmc.integrated)

#RNA_rawcounts <- as.data.frame(pbmc.integrated@assays$RNA@counts)
RNA_data <- as.data.frame(pbmc.integrated@assays$RNA@data)

#SCT_rawcounts <- as.data.frame(pbmc.integrated@assays$SCT@counts)
#SCT_data <- as.data.frame(pbmc.integrated@assays$SCT@data)

RNAsumhist <- as.data.frame(rowSums(RNA_data))
#RNAsumhist_raw <- as.data.frame(rowSums(RNA_rawcounts))
#SCTsumhist <- as.data.frame(rowSums(SCT_data))
#SCTsumhist_raw <- as.data.frame(rowSums(SCT_rawcounts))

colnames(RNAsumhist) <- "RNA_depth"
#colnames(RNAsumhist_raw) <- "RNA_depth_raw"
#colnames(SCTsumhist) <- "SCT_depth"
#colnames(SCTsumhist_raw) <- "SCT_depth_raw"

#p <- ggplot(RNAsumhist, aes(x=RNA_depth)) + geom_histogram()
#p <- p + scale_x_continuous(trans='log10') + ggtitle("RNA data after log normalize")
#p

#p1 <- ggplot(RNAsumhist_raw, aes(x=RNA_depth_raw)) + geom_histogram()
#p1 <- p1 + scale_x_continuous(trans='log10') + ggtitle("RNA data before log normalize")
#p1

#p + p1

#p2 <- ggplot(SCTsumhist, aes(x=SCT_depth)) + geom_histogram()
#p2 <- p2 + scale_x_continuous(trans='log10') + ggtitle("SCT data after log normalize")
#p2

#p3 <- ggplot(SCTsumhist_raw, aes(x=SCT_depth_raw)) + geom_histogram()
#p3 <- p3 + scale_x_continuous(trans='log10') + ggtitle("SCT data before log normalize")
#p3

#p2 + p3

#select data for iTalk analysis
#note: MAST request log normalized data, and DESEQ request raw RNA counts
tdata <- t(RNA_data)
tdata <- as.data.frame(tdata)

#check the distribution of gene, see if there are two peaks (because of batches)
p <- ggplot(tdata, aes(x=CD14)) + geom_histogram(bins = 30)
p <- p + scale_y_continuous(trans='log10')
p

#before adding cellID, needed to make sure that rowname of metadata is the same as rowname of tdata
identical(row.names(pbmc.integrated@meta.data),row.names(tdata))

#check the data you want to add
table(pbmc.integrated@meta.data$celltype)

#add metadata to the matrix
dim(tdata)
tdata$cell_type <- pbmc.integrated@meta.data$celltype
dim(tdata)

tdata$compare_group <- pbmc.integrated@meta.data$condition

tdata$compare_group[tdata$compare_group=="CONT"] <- 1
tdata$compare_group[tdata$compare_group=="PAH"] <- 2


#select data for ligand receptor analysis
data <-dplyr::filter(tdata, !cell_type %in% c("X27_UK", "Megakaryocytes"))

data <-dplyr::filter(tdata, !cell_type %in% c("X1_CD14Mono","X12_CD14Mono_InterM"))

data <-dplyr::filter(tdata, cell_type %in% c("X0_CD4_1", "X1_CD14Mono","X12_CD14Mono_InterM","X13_CD4_FOXP3","X14_CD4CD8"))

dim(data)
data[35288:35298,20738:20748]

#########################################################################################################
#control
#########################################################################################################
data_control <- filter(data, compare_group == 1)
dim(data_control)

# find top [] percent highly expressed genes
highly_exprs_genes.C <-rawParse(data_control,top_genes=40,stats='mean')

table(highly_exprs_genes.C$cell_type)

highly_exprs_genes.C[highly_exprs_genes.C$gene == "CD8A",]

highly_exprs_genes.C[highly_exprs_genes.C$gene == "CD14",]

highly_exprs_genes.C[highly_exprs_genes.C$gene == "CD19",]

highly_exprs_genes.C[highly_exprs_genes.C$gene == "CD3D",]

highly_exprs_genes.C[highly_exprs_genes.C$gene == "IL6",]

# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')

cell_col<-structure(c('#4a84ad','#b57b52', '#a6d608', '#5929f7','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a',
                      '#4a84ad','#b57b52', '#a6d608', '#5929f7','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a',
                      '#4a84ad','#b57b52', '#a6d608', '#5929f7','#4a1dc6','#e874bf'),names=unique(data$cell_type))

par(mfrow=c(1,2))
res.c<-NULL

for(comm_type in comm_list){
  res_cat<-FindLR(highly_exprs_genes.C,datatype='mean count',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  #plot by ligand category
  #overall network plot
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  #top 20 ligand-receptor pairs
  LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
  title(comm_type)
  res.c<-rbind(res.c,res_cat)
}

#########################################################################################################
#patient
#########################################################################################################
data_pah <- filter(data, compare_group == 2)
dim(data_pah)

# find top [] percent highly expressed genes
highly_exprs_genes.P <-rawParse(data_pah,top_genes=10,stats='mean')

table(highly_exprs_genes.P$cell_type)

highly_exprs_genes.P[highly_exprs_genes.P$gene == "CD8A",]

# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')

cell_col<-structure(c('#4a84ad','#b57b52', '#a6d608', '#5929f7','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a',
                      '#4a84ae','#b57b52', '#a6d608', '#5929f7','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a',
                  '#4a84af','#b57b52', '#a6d608', '#5929f7','#4a1dc6','#e874bf'),names=unique(data$cell_type))

c25 <- c("dodgerblue2", "#E31A1C", # red
  "green4", "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2", "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2","yellow",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown")

cell_col<-structure(c25,names=unique(data$cell_type))

par(mfrow=c(1,2))
res.p <- NULL
for(comm_type in comm_list){
  res_cat<-FindLR(highly_exprs_genes.P,datatype='mean count',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  #plot by ligand category
  #overall network plot
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  #top 20 ligand-receptor pairs
  LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
  title(comm_type)
  res.p<-rbind(res.p,res_cat)
}

unique(res$ligand)
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:80,]


#########################################################################################################
#plot specific ligand receptor pair
#########################################################################################################
ligand_S <- "APP"
res.c1 <- subset(res.c, ligand == ligand_S| receptor == ligand_S)
res.p1 <- subset(res.p, ligand == ligand_S| receptor == ligand_S)

par(mfrow=c(1,2))
NetView(res.c1,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
NetView(res.p1,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

pdf(file = "~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/italksPBMC+DC/iTALK.CD14.pdf", width = 16,height = 10)
par(mfrow=c(1,2))
LRPlot(res.c1[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res.c1$cell_from_mean_exprs[1:20],link.arr.width=res.c1$cell_to_mean_exprs[1:20],
       text.vjust="0.5cm", track.height_2=0.1,annotation.height_1 = 0.04,annotation.height_2 = 0.01)
LRPlot(res.p1[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res.p1$cell_from_mean_exprs[1:20],link.arr.width=res.p1$cell_to_mean_exprs[1:20],
       text.vjust="0.5cm", track.height_2=0.1,annotation.height_1 = 0.04,annotation.height_2 = 0.01)
title(paste(ligand_S, "Singnalling[Control-PAH]", sep=" "))
dev.off()

#########################################################################################################
#plot all ligand receptor pair
#########################################################################################################
plotlist<- unique(res.c$ligand)
plotlist <- plotlist[plotlist!="VEGFB"]
plotlist <- plotlist[plotlist!="PF4"]
plotlist <- plotlist[plotlist!="C1QA"]
plotlist <- plotlist[plotlist!="ADAM15"]
plotlist <- plotlist[plotlist!="NCAM1"]
plotlist <- plotlist[plotlist!="FASLG"]
plotlist <- plotlist[plotlist!="FLT3LG"]
plotlist <- plotlist[plotlist!="CCL4"]
plotlist <- plotlist[plotlist!="CXCL16"]
plotlist <- plotlist[plotlist!="CCL3"]
plotlist <- plotlist[plotlist!="CD24"]
plotlist <- plotlist[plotlist!="LGALS9"]

for (i in plotlist) {
  print(i)
  ligand_S <- i
  res.c1 <- subset(res.c, ligand == ligand_S| receptor == ligand_S)
  res.p1 <- subset(res.p, ligand == ligand_S| receptor == ligand_S)

  title <- paste(ligand_S, "iTALK.pdf", sep="_")

  pdf(file = paste("~/Box/XinZhouFiles/Projects/ZXP1_PAH/PBMC/Alternative Analysis (loose cutoff)/italksPBMC+DC/", title, sep=""), width = 16,height = 10)
  par(mfrow=c(1,2))
  LRPlot(res.c1[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res.c1$cell_from_mean_exprs[1:20],link.arr.width=res.c1$cell_to_mean_exprs[1:20],
       text.vjust="0.5cm", track.height_2=0.1,annotation.height_1 = 0.04,annotation.height_2 = 0.01)
  LRPlot(res.p1[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res.p1$cell_from_mean_exprs[1:20],link.arr.width=res.p1$cell_to_mean_exprs[1:20],
       text.vjust="0.5cm", track.height_2=0.1,annotation.height_1 = 0.04,annotation.height_2 = 0.01)
  title(paste(ligand_S, "Singnalling[Control-PAH]", sep=" "))
  dev.off()
}


#########################################################################################################
#two conditon
#########################################################################################################
tdata$compare_group <- pbmc.integrated@meta.data$condition

tdata$compare_group[tdata$compare_group=="CONT"] <- 1
tdata$compare_group[tdata$compare_group=="PAH"] <- 2

#MAST for RNA data is not accurate, use SCT instead
X1 <- tdata %>% filter(cell_type=='X1_CD14Mono')
X1_Clean <- X1[(colSums(X1[,1:20746]) != 0),]
dim(X1_Clean)

deg_monoX1<-DEG(tdata %>% filter(cell_type=='X1_CD14Mono'), method='Wilcox',contrast=c(2,1))

#deg_monoX1<-DEG(tdata %>% filter(cell_type=='X1_CD14Mono'), method='edgeR',contrast=c(2,1))

deg_monoX4<-DEG(tdata %>% filter(cell_type=='X17_CD8_RORC'), method='Wilcox',contrast=c(2,1))

deg_monoX1[deg_monoX1$gene=="TLR4",]

par(mfrow=c(1,2))
res<-NULL

for(comm_type in comm_list){
  res_cat<-FindLR(deg_monoX1,deg_monoX4,datatype='DEG',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
  #plot by ligand category
  if(nrow(res_cat)==0){
    next
  }else if(nrow(res_cat)>=20){
    LRPlot(res_cat[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC[1:20],link.arr.width=res_cat$cell_to_logFC[1:20])
  }else{
    LRPlot(res_cat,datatype='DEG',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_logFC,link.arr.width=res_cat$cell_to_logFC)
  }
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  title(comm_type)
  res<-rbind(res,res_cat)
}

if(is.null(res)){
  print('No significant pairs found')
}else if(nrow(res)>=20){
  res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:20,]
  NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  LRPlot(res[1:20,],datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC[1:20],link.arr.width=res$cell_to_logFC[1:20])
}else{
  NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  LRPlot(res,datatype='DEG',cell_col=cell_col,link.arr.lwd=res$cell_from_logFC,link.arr.width=res$cell_to_logFC)
}

