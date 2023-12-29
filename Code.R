####Univariate COX regression####
library(survival)
td=read.table("data.txt",header = T,row.names = 1,sep = "\t")
pFilter=0.05 
outResult=data.frame() 
sigGenes=c("status","time")
for(i in colnames(td[,3:ncol(td)])){ 
  tdcox <- coxph(Surv(time, status) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox) 
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"] 
  if(pvalue<pFilter){ 
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outResult,file="UniCoxSurvival.txt",sep="\t",row.names=F,quote=F)
UniCoxSurSigGeneExp=td[,sigGenes] 
UniCoxSurSigGeneExp=cbind(id=row.names(UniCoxSurSigGeneExp),UniCoxSurSigGeneExp)
write.table(UniCoxSurSigGeneExp,file="UniCoxSurSigGeneExp.txt",sep="\t",row.names=F,quote=F)
tducs <- read.table("UniCoxSurvival.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(tducs)
hr <- sprintf("%.3f",tducs$"HR")
hrLow  <- sprintf("%.3f",tducs$"L95CI")
hrHigh <- sprintf("%.3f",tducs$"H95CI")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pValue <- ifelse(tducs$pvalue<0.001, "<0.001", sprintf("%.3f", tducs$pvalue))

pdf(file="UniCoxSurForestPlot.pdf", width = 10,height = 20)
n <- nrow(tducs)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(2.5,2))

xlim = c(0,2.5)
par(mar=c(4,2.5,2,1))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.2-0.5*0.2,n:1,pValue,adj=1,cex=text.cex);text(1.2-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(2.5,n:1,Hazard.ratio,adj=1,cex=text.cex);text(2.5,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(0,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 1, col = boxcolor, cex=1.3)
axis(1)

dev.off()
####KM survival curve####
library(dplyr)
library(survival)
library(survminer)
score_2 <- read.csv("Clinical.csv")  
score_2 <- read.table("clinical.txt",header = T,sep = "\t")
Sscore <- score_2
sfit <- surv_fit(Surv(time,status==1)~score,data = Sscore)
sfit
sfit2 <- survfit(Surv(time,status==1)~score,data=Sscore)
sfit2
summary(sfit2,times = 1)
survdiff(Surv(time,status==1)~score ,data = Sscore)
coxph(Surv(time,status==1)~score,data = Sscore)
pairwise_survdiff(Surv(time,status==1)~score,data=Sscore)
pdf("TCGA-LUAD.pdf", height = 5.5,width = 6,onefile = FALSE)
ggsurvplot(sfit2,conf.int = F,pval = T,
           risk.table = T,surv.median.line = "hv",
           legend = c(0.8, 0.83),
           legend.labs=c("C1","C2"),break.time.by=500,
           legend.title="Group",palette = c('#1F78B4',"#D20A13"),
           title="TCGA-LUAD",ggtheme=theme_classic()+
           theme(plot.title = element_text(size=12,hjust=0.5)))+
           xlab("Time(Days)")#"#E7B800","#2E9FDF"
dev.off()
####pam clustering#### 
library(ConsensusClusterPlus)
title <-"D:/TRP/step2-Pam clustering"
rt=read.table("pam.txt", header=T, sep="\t", check.names=F, row.names=1)
dataset=as.matrix(rt)
mads <- apply(dataset, 1, mad)
dataset <- dataset[rev(order(mads))[1:14],]
dim(dataset)
results <- ConsensusClusterPlus(dataset, maxK = 6,
                                reps = 1000, pItem = 0.8,
                                pFeature = 1,  
                                clusterAlg = "km", 
                                distance = "pearson",
                                title = title,
                                plot = "png")#data为matrix
sample_cluster <- results[[2]]$consensusClass
sample_cluster_df <- data.frame(sample = names(sample_cluster),
                                cluster = sample_cluster)
head(sample_cluster_df)
write.table(sample_cluster_df, file="cluster2.txt", sep="\t", quote=F, col.names=T)
####GSEA enrichment####
library(limma)
library(dplyr)
df <- read.table("exp.txt", header = T, sep = "\t", row.names = 1, check.names = F)
list <- c(rep("cluster1",256),rep("cluster2",244)) %>% factor(., levels = c("cluster1", "cluster2"), ordered = F)
list <- model.matrix(~factor(list)+0)  
colnames(list) <- c("cluster1", "cluster2")
df.fit <- lmFit(df, list) 
df.matrix <- makeContrasts(cluster1 - cluster2 , levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,n = Inf, adjust = "fdr")
head(tempOutput)
nrDEG = na.omit(tempOutput)
diffsig <- nrDEG  
write.csv(diffsig, "all.limmaOut.csv")
d<-read.csv("all.limmaOut.csv") 
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
gene<-str_trim(d$symbol,"both") 
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=d$logFC, 
                      SYMBOL = d$symbol) 
gene_df <- merge(gene_df,gene,by="SYMBOL")

geneList<-gene_df $logFC 
names(geneList)=gene_df $ENTREZID 
geneList=sort(geneList,decreasing = T) 
kegmt<-read.gmt("h.all.v2023.2.Hs.entrez.gmt") 
KEGG<-GSEA(geneList,TERM2GENE = kegmt)
####boxplot and Violin diagrams####
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(ggsci)
df <- read_tsv("check.txt") %>% pivot_longer(-group)
p=ggplot(df,aes(name,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,outlier.fill  = 'black',outlier.size = 0.25,color = "black") + 
  labs(x = " ", y = "Immune infiltration") +#Gene expression, Immune infiltration
  theme(legend.position = "top") + 
  theme_bw() +
  theme(legend.position="top",
        panel.grid.major.x = element_line(size = 1),
        panel.grid.major.y = element_line(size = 1),
        panel.grid.minor.x = element_line(size = 0.8),
        panel.grid.minor.y = element_line(size = 0.8),
        
        axis.text.x = element_text(angle = 50, hjust = 1,vjust = 1, colour = "black")) + 
  theme(text = element_text(size=12,face = "bold")) + 
  scale_fill_manual(values= c("#7AC7E2","#E3716E"))+ 
  stat_compare_means(aes(group = group,label = ..p.signif..),method = "wilcox.test",)#kruskal.test
p
ggsave(p,file='ssGSEA.pdf',width = 10,height = 6.5)
fig2e<-ggboxplot(df, x='name', y='value', 
                 fill = "group", color = "black",
                 palette = ggsci::pal_aaas()(9), 
                 ylab="Gene Expression",xlab='',
                 add = "boxplot")+ 
  stat_compare_means(aes(group=group),method = 'wilcox.test',
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")+
  theme(axis.text.x = element_text(angle = 50,hjust = 1,face = "bold"))
fig2e
ggsave('Fig1a.pdf',fig2e,height = 5,width = 9)
df <- read_tsv("CIBERSORT-Results.txt") %>% pivot_longer(- group)
location <- df %>% group_by(name) %>% slice_max(value)
location$x <- seq(1,22,by=1)
head(location,3)
ggplot(df,aes(x = name, y = value,fill=group))+
      geom_violin(scale = "width",alpha=0.8,width=0.5,size=0.5)+
      scale_fill_manual(values = c("#7AC7E2","#E3716E"))+        
      stat_compare_means(aes(group=group),                       
                     method = "wilcox.test",
                     paired = F,                             
                     symnum.args = list(cutpoint=c(0,0.001,0.01,0.05,1),
                                        symbols=c("***","**","*","ns")),
                     label = "p.signif",
                     label.y = location$value+0.02,      
                     size=4.5)+ 
      geom_segment(data=location,                                
                  aes(x=x,y=value,xend=x+0.2,yend=value),size=1)+
      xlab("")+                                                  
      ylab("Fraction")+                                           
      theme_bw()+
      theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(size=1.2),                 
        axis.text.x = element_text(angle=50,size=10,vjust = 1,hjust =1,color = "black",face = "bold"),                                                   
        axis.text.y = element_text(size =10,color = "black",face = "bold"),
        legend.position = c(0.9,0.85) )+
  geom_line(position=position_dodge(0.5))+
  stat_summary(fun = "median", geom = "point", position=position_dodge(0.5),colour="white",size=3)
####GSVA enrichment####
library(GSVA)
gs = read.csv("cellmarker.csv", stringsAsFactors = FALSE, check.names = FALSE)
a = read.table("tcga_dat_T.txt", stringsAsFactors = FALSE, header = TRUE, row.names = 1, sep = "\t")
a = as.matrix(a)
gs = as.list(gs)
gs = lapply(gs, function(x) x[!is.na(x)])
ssgsea = gsva(a, gs, method = "gsva", kcdf='Gaussian',abs.ranking=TRUE)   # signature 'matrix,list'
write.csv(ssgsea, "ssgsea-kegg.csv")
####genomic alteration####
library(maftools)
library(RColorBrewer)
options(stringsAsFactors = F)
laml = read.maf(maf = 'combined_maf_value.txt')#合并后的体细胞突变数据矩阵
pdf(file="maf_summary.pdf",width =12,height=7)
plotmafSummary(maf = laml,addStat = 'median',color = col)
dev.off()
cols <- brewer.pal(8,"Set1")
names(cols) <- levels(laml@data$Variant_Classification)
cols = RColorBrewer::brewer.pal(n = 10, name = 'Paired')
names(cols) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Ins','In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Nonstop_Mutation','Translation_Start_Site','Multi_Hit')
genes <- c("OGDH","INMT","GCDH","IDO2","MAOA","ALDH2","ACAT1","MAOB","CAT","ACAT2","ASMT","KYNU","SLC7A5","SLC3A2")
pdf(file="TMB_cluster.pdf",width =8,height=8)
oncoplot(maf = laml, genes = genes,colors = cols,draw_titv = T)#top = 30
dev.off()
tmb_table=tmb(maf = laml)
tmb_table=tmb(maf = laml,logScale = F)
write.csv(tmb_table,"tmb_results.csv")
library(maftools)
g<- maftools::readGistic(
  gisticAllLesionsFile = "all_lesions.conf_95.txt",
  gisticAmpGenesFile = "amp_genes.conf_95.txt",
  gisticDelGenesFile = "del_genes.conf_95.txt",
  gisticScoresFile = "scores.gistic",
  isTCGA = TRUE,
  verbose = FALSE)
pdf("ChromPlot-low.pdf",width = 12,height = 6, onefile = TRUE)
maftools::gisticChromPlot(
  gistic = g,
  fdrCutOff = 0.25,
  txtSize = 0.8,
  cytobandTxtSize = 0.5,
  color = c("#D95F02","#1B9E77"),
  markBands = "all",
  ref.build ="hg38",
  y_lims = c(-0.5,0.5)
  )
dev.off()
####scRNA-seq data analysis ####
rm(list = ls())
options(stringsAsFactors = F)
setwd("C:\\Users\\86189\\Desktop\\代码\\55\\202200925.PAAD.LR\\01.scRNA")
dir.create('results',recursive = T)
library(devtools)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
dir.create('results')
file_id <- list.files("C:\\Users\\86189\\Desktop\\20230225.LUAD.TRP\\01.scRNA", full=FALSE)
datalist=list()
for (i in 1:length(file_id)){
  dir.10x = paste0("C:\\Users\\86189\\Desktop\\20230225.LUAD.TRP\\01.scRNA\\",file_id[i])
  my.data <- Read10X(data.dir = dir.10x) 
  colnames(my.data)=paste0(file_id[i],'_',colnames(my.data))
  datalist[[i]]=CreateSeuratObject(counts = my.data, 
                                   project = file_id[i],min.cells = 3, min.features = 250)
  datalist[[i]]$GSM=file_id[i]
}
rm(my.data)
names(datalist)=file_id

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
raw_meta=sce@meta.data
raw_count <- table(raw_meta$GSM)
raw_count
pearplot_befor<-VlnPlot(sce,group.by ='GSM', 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0, 
                        ncol = 4)
pearplot_befor
Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA',feature2 = 'nCount_RNA',group.by = 'GSM')
Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt',feature2 = 'nCount_RNA',group.by = 'GSM')
Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt',feature2 = 'nFeature_RNA',group.by = 'GSM')
Feature_ber1=Feature_ber1+theme(legend.position = 'none')
Feature_ber2=Feature_ber2+theme(legend.position = 'none')
Feature_ber<-ggarrange(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths = c(1,1,1.2))

datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 100 & 
              nFeature_RNA < 8000 & 
              percent.mt < 35 &
              nCount_RNA > 1000 )
})
#save(datalist,file = 'datalist.RData')
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_meta=sce@meta.data
clean_count <- table(clean_meta$GSM)
clean_count
pearplot_after <- VlnPlot(sce,group.by ='GSM', 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                          pt.size = 0, 
                          ncol = 4)
pearplot_after
gc()
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))

sce <- ScaleData(sce, features =  rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce)) 
dimplot1 <- DimPlot(sce, reduction = "pca",group.by = 'GSM') 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca
dev.off()
figs1abc=ggarrange(Feature_ber,pearplot_befor,pearplot_after,nrow = 3,labels = c('A','B','C'),ncol = 1)
figs1<-ggarrange(figs1abc,sc_pca,nrow = 2,ncol = 1,labels = c('','D'),heights = c(2,1))
figs1
ggsave(filename = 'results/FigS1.pdf',plot = figs1,he=20,width = 20)
dev.off()
sum(raw_count)
sum(clean_count)
Dims <- 40
sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)
#T cell：
VlnPlot(sce,features = c('CD2','CD3D','CD3E','CD3G'),pt.size = 0,group.by = 'seurat_clusters')
#CD8 T :
VlnPlot(sce,features = c('CD3D','CD8A','CD8B','GZMA'),pt.size = 0,group.by = 'seurat_clusters')
#B cells: 
VlnPlot(sce,features = c('CD19','CD79A','MS4A1'),pt.size = 0,group.by = 'seurat_clusters')
#Mast cell:
VlnPlot(sce,features = c('TPSAB1','CPA3'),pt.size = 0,group.by = 'seurat_clusters')
#Monocyte :
VlnPlot(sce,features = c('S100A12','VCAN','FCN1','S100A8'),pt.size = 0,group.by = 'seurat_clusters')
#Macrophage:
VlnPlot(sce,features = c('CD163','CD68','CD14'),pt.size = 0,group.by = 'seurat_clusters')
#Plasma cell:
VlnPlot(sce,features = c('CD79A','JSRP1'),pt.size = 0.1,group.by = 'seurat_clusters')
#epithelial cells (EPCAM):
VlnPlot(sce,features = 'EPCAM',pt.size = 0,group.by = 'seurat_clusters')
#proliferating cells (MKI67) 
VlnPlot(sce,features = 'MKI67',pt.size = 0,group.by = 'seurat_clusters')
#Fibroblast:
VlnPlot(sce,features = c('ACTA2','FAP','PDGFRB','NOTCH3'),pt.size = 0,group.by = 'seurat_clusters')
#endothelial cells:
VlnPlot(sce,features = c('PECAM1'),pt.size = 0,group.by = 'seurat_clusters')
#Mesenchymal stromal cell 8,10
VlnPlot(sce,features = c('ENG', 'CD44', 'NT5E', 'THY1'),pt.size = 0,group.by = 'seurat_clusters')
#Monocytes
VlnPlot(sce,features = c('CD14', 'CD68', 'CD74'),pt.size = 0,group.by = 'seurat_clusters')
library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
library(stringr)  
p_all_markers <- DotPlot(sce, features = genes_to_check,
                         assay='RNA'  )  + coord_flip()
p_all_markers
feat_gene<-c('PECAM1','TPSAB1','CPA3','CD3D','CD8A','CD8B','GZMA','CD79A','MS4A1','JSRP1','EPCAM','ACTA2','FAP','PDGFRB','NOTCH3','CD68','RUNX2', 'ALPL', 'IBSP')
cell_type=c('Epithelial','T_cell','Endothelial','Cancer','Macrophage','Mast_cell','Cancer','Macrophage','Fibroblasts','Epithelial','Cancer','Fibroblasts','Endothelial','B_cell','Epithelial','T_cell')
Idents(sce) <- sce@meta.data$seurat_clusters
names(cell_type) <- levels(sce)
sce<- RenameIdents(sce, cell_type)
sce@meta.data$cell_type <- Idents(sce)
Idents(sce)=sce@meta.data$cell_type
table(sce$cell_type)
length(feat_gene)
icd=read.delim('S1.txt',sep='\t',header = F)
feat_gene=icd$V1
pdf('results/FigS2-3.pdf',he=12,width =15)
FeaturePlot(sce,
            features = feat_gene,
            pt.size = 0.1,reduction = 'tsne',ncol = 5)
dev.off()
save(sce,file = 'results/sce.RData')
load("results/sce.RData")
library(scMetabolism)
library(ggplot2)
library(rsvd)
load(file = "sce.RData")
countexp.Seurat<-sc.metabolism.Seurat(obj = sce, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "tsne", dimention.reduction.run = F, size = 1)
BoxPlot.metabolism(obj = countexp.Seurat, pathway = input.pathway, phenotype = "ident", ncol = 1)
####Immune infiltration analysis####
rm(list = ls())
library(estimate)
in.file <- 'tcga_dat_T.txt' 
outfile2E <- 'ESTIMATE_input.gct' 
outputGCT(in.file, outfile2E) 
filterCommonGenes(input.f= in.file, output.f= outfile2E, id="GeneSymbol")
estimateScore("ESTIMATE_input.gct", "ESTIMATE_score.gct")
plotPurity(scores="ESTIMATE_score.gct", samples="s516")
ESTIMATE_score <- read.table("ESTIMATE_score.gct", skip = 2,
                             header = TRUE,row.names = 1) 
ESTIMATE_score <- ESTIMATE_score[,2:ncol(ESTIMATE_score)] 
ESTIMATE_score 
write.table(ESTIMATE_score,file = "ESTIMATE_score.txt",quote = F,sep = "\t")
remove(list = ls()) 
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
LM22.file <- "LM22.txt"
TCGA_exp.file<- "tcga_dat_T.txt"
TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 1000, QN = F)  
write.csv(TCGA_TME.results, "TCGA_CIBERSORT_Results_fromRcode.csv")










