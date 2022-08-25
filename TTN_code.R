##counts & TPM download
rm(list = ls())
library(tidyverse) # 数据清洗与读取的利器
library(stringr)
library(dplyr)
library(limma)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(tidyverse)
library(stringr)
library(dplyr)
counts=read.csv('counts.csv',row.names = 1,check.names = FALSE)
counts$gene_id=substr(counts$gene_id,1,15)
output=bitr(counts$gene_id,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = org.Hs.eg.db,drop = FALSE)
output$gene_id=output$ENSEMBL
counts_all <- counts %>% 
  inner_join(output,by="gene_id")
counts_all$gene_id=counts_all$SYMBOL
counts=counts_all[,-c(ncol(counts_all)-1,ncol(counts_all))]
eliminate_duplicated_lines=function(x){
  counts_df=x
  counts_df$median=apply(counts_df[,-1],1,median)
  counts_df=counts_df[order(counts_df$gene_id,counts_df$median,decreasing = T),]
  counts_df=counts_df[!duplicated(counts_df$gene_id),]
  counts_df=counts_df[,-c(ncol(counts_df))]
  counts_df=subset(counts_df,gene_id!= 'NA')
  x=counts_df
}
counts=eliminate_duplicated_lines(counts)

##TPM
gene_length<- read.table('All_hg19gene_len.txt',header = TRUE)
colnames(gene_length)=c('gene_id','Length')
merge<-inner_join(counts,gene_length,by="gene_id")#根据基因那列进行合并
merge <- as.data.frame(na.omit(merge))#删除错误值行
rownames(merge)=merge$gene_id
mycounts=merge[,-1]  #删除GENE列
mycounts$Length=as.numeric(mycounts$Length)
kb <-( mycounts$Length) / 1000
countdata <- mycounts[,1:(ncol(mycounts)-1)]
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
log_counts_tpm=log2(tpm+1)
normalized_log_counts_tpm=as.data.frame(normalizeBetweenArrays(log_counts_tpm))
rownames(normalized_log_counts_tpm)=rownames(log_counts_tpm)
sample=rownames(read.csv('stageIII_IV_sample_name.csv',row.names = 1,header = FALSE))
normalized_log_counts_tpm=normalized_log_counts_tpm[,substr(colnames(normalized_log_counts_tpm),1,12) %in% sample]
normalized_log_counts_tpm=normalized_log_counts_tpm[,-which(duplicated(substr(colnames(normalized_log_counts_tpm),1,12)))]

write.csv(counts,'counts.csv')
write.csv(normalized_log_counts_tpm,'counts_tpm.csv')

#maf
rm(list = ls())
library(maftools)
maf=read.maf('TCGA.SKCM.varscan.6c961926-7792-42fa-9a16-c62f60e2557b.DR-10.0.somatic.maf.gz',
             clinicalData = 'TCGA-SKCM-clinical.csv',isTCGA = TRUE)
sample=rownames(read.csv('stageIII_IV_sample.csv',row.names = 1))
maf=read.maf(maf = maf@data[maf@data$Tumor_Sample_Barcode %in% sample,],
             clinicalData = 'TCGA-SKCM-clinical_select2.csv',isTCGA = TRUE)
plotmafSummary(maf = maf,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE,titvRaw = FALSE)
hnsc=titv(maf,plot = TRUE,useSyn = TRUE)
lollipopPlot(maf = maf, gene = 'TTN', showMutationRate = TRUE,AACol = "HGVSp_Short",pointSize = 0.8,printCount = TRUE,
             collapsePosLabel = TRUE,showDomainLabel = FALSE,legendTxtSize =1.1,domainBorderCol = NA,domainAlpha =0.8,
             roundedRect = TRUE,defaultYaxis = TRUE)
vc_cols =c("#1F78B4", "#33A02C", "#83a78d",'#0eb0c9' ,"#FB9A99" ,"#E31A1C" ,"#FDBF6F" ,"#FF7F00")
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = maf, colors = vc_cols, top = 25, bgCol = 'white',
         clinicalFeatures = c("Gender",'Age','Stage','N','M') ,annoBorderCol ="white",anno_height = 3,
         logColBar = FALSE,legendFontSize =1)
#TMB
maf_TMB=tmb(maf=maf)
write.csv(maf_TMB,'skcm_tmb.csv')

####DEGs
rm(list=ls())
stringsAsFactors = FALSE
library(DESeq2)
library(ggplot2)
counts=read.csv('counts.csv',row.names = 2,check.names = FALSE)
counts=counts[,-1]
TTN_mution_sample=read.csv('TTN突变样本.csv',header = FALSE)$V1
TTN_no_mutation_sample=read.csv('TTN未突变样本.csv',header = FALSE)$V1
counts=counts[,as.numeric(substr(colnames(counts),14,15)) < 10]#保留肿瘤样本
mat1 <- counts[,substr(colnames(counts),1,12) %in% TTN_mution_sample]#TTN突变的表达矩阵
mat2 <- counts[,substr(colnames(counts),1,12) %in% TTN_no_mutation_sample]#TTN未突变的表达矩阵
mat2=mat2[,-which(duplicated(substr(colnames(mat2),1,12)))]
gene <- cbind(mat2,mat1)
gene$sum=apply(gene,1,sum)
gene=gene[gene$sum != 0,]
gene=gene[,-ncol(gene)]
gene_1=apply(gene, 1,as.integer)#所有格式统一为integer
rownames(gene_1)=colnames(gene)
gene=as.data.frame(t(gene_1))
rm(gene_1)
group = factor(c(rep('TTN_WT',ncol(mat2)),rep('TTN_MUT',ncol(mat1))))
coldata <- data.frame(row.names = colnames(gene),group = factor(c(rep('TTN_WT',ncol(mat2)),rep('TTN_MUT',ncol(mat1)))))
dds <- DESeqDataSetFromMatrix(countData = gene, colData = coldata, design = ~group)
dds <- DESeq(dds, parallel = FALSE)  
suppressMessages(dds)
res <- results(dds, contrast = c('group',  'TTN_WT','TTN_MUT'), pAdjustMethod = "BH", alpha = 0.05)
deseq_res <- as.data.frame(res[order(res$padj), ])
deseq_res$gene_id <- rownames(deseq_res)
write.table(deseq_res[c(7, 1:6)], 'DESeq2.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.csv(deseq_res[c(7, 1:6)], 'DEG_result.csv')

deseq_res <- read.delim('DESeq2.txt', sep = '\t')
deseq_res[which(deseq_res$log2FoldChange >= 1 & deseq_res$padj < 0.05),'sig'] <- 'log2FC >= 1'
deseq_res[which(deseq_res$log2FoldChange <= -1 & deseq_res$padj < 0.05),'sig'] <- 'log2FC <= -1'
deseq_res[which(abs(deseq_res$log2FoldChange) < 1 & deseq_res$padj < 0.05),'sig'] <- 'p-value < 0.05'
deseq_res[which(abs(deseq_res$padj) >= 0.05),'sig'] <- 'no difference'
deseq_res[which(deseq_res$padj %in% NA),'sig'] <- 'no difference'
volcano_p <- ggplot(deseq_res, aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("#2983bb","#ee3f4d" ,'gray30' ,'#1ba784')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.position = c(0.26, 0.92)) +
  theme(legend.title = element_blank(),legend.text = element_text(size = 14), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size=14))+
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value', color = NA,size=10) +
  xlim(-5, 5)
volcano_p



####GSEA
rm(list = ls())
data=read.csv('DEG_result.csv')
data$SYMBOL=data$gene_id
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)
gene <- data$gene_id
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
data_all <- data %>% 
  inner_join(gene,by="SYMBOL")
data_all_sort <- data_all %>% 
  arrange(desc(log2FoldChange))
head(data_all_sort)
geneList = data_all_sort$log2FoldChange 
names(geneList) <- data_all_sort$ENTREZID 
head(geneList)
gsemf <- gseGO(
  geneList, 
  ont = "BP", 
  OrgDb = org.Hs.eg.db, 
  keyType = "ENTREZID",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
)

gseaplot(gsemf,1)
write.csv(gsemf,"GSEAresult.csv")
library(enrichplot)
path_1=c("GO:0031341",
         "GO:0002920",
         "GO:0001909",
         "GO:0002703",
         "GO:0002440",
         "GO:0070661",
         "GO:0050863",
         "GO:0006959",
         "GO:0002697")
gseaplot2(gsemf,path_1,color = colorspace::rainbow_hcl(4),subplots=c(1,2),base_size = 15)

####WGCNA
rm(list = ls())
normalized_counts_tpm=read.csv('counts_tpm.csv',row.names = 1,check.names = FALSE)
DEGs_result=read.csv('DEG_result.csv')
DEGs=unique(DEGs_result[DEGs_result$padj < 0.05,]$X)
library(MASS)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(WGCNA)
library(stringr)
options(stringsAsFactors = F)
library(DESeq2)
enableWGCNAThreads()
datExprdataOne=normalized_counts_tpm
datExprdataOne=datExprdataOne[rownames(datExprdataOne) %in% DEGs,]
datExprdataOne=as.data.frame(t((datExprdataOne)))
gsg = goodSamplesGenes(datExprdataOne, verbose = 3)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExprdataOne)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExprdataOne)[!gsg$goodSamples], collapse = ",")));
  datExprdataOne = datExprdataOne[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExprdataOne), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers1 <- c(seq(1, 10, by=1), seq(12, 20, by=2))
sft <- pickSoftThreshold(datExprdataOne, powerVector = powers1)
RpowerTable <- pickSoftThreshold(datExprdataOne, powerVector = powers1)[[2]]
cex1 = 0.9
par(mfrow = c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], xlab = "soft threshold (power)", ylab = "scale free topology model fit, signes R^2", type = "n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels = powers1, cex = cex1, col = "red")
abline(h = 0.9, col = "red")
plot(RpowerTable[,1], RpowerTable[,5], xlab = "soft threshold (power)", ylab = "mean connectivity", type = "n")
text(RpowerTable[,1], RpowerTable[,5], labels = powers1, cex = cex1, col = "red")

betal = sft$powerEstimate+1
k.dataOne <- softConnectivity(datExprdataOne, power = betal) -1
par(mfrow=c(2,2))
scaleFreePlot(k.dataOne, main = paste("data set I, power=", betal), truncated = F)
cor <- WGCNA::cor
dataOne_net = blockwiseModules(datExprdataOne, power = betal, maxBlockSize =2000,
                               TOMType = "unsigned",  deepSplit = 4,minModuleSize =50,
                               reassignThreshold = 0, mergeCutHeight = 0.25,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs=TRUE, corType = "pearson", 
                               loadTOMs=TRUE,
                               saveTOMFileBase = "DataI.tom",
                               verbose = 3)
table(dataOne_net$colors)
moduleLabels = dataOne_net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(dataOne_net$dendrograms[[1]], moduleColors[dataOne_net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dataOne_MEs <- dataOne_net$MEs
plotEigengeneNetworks(dataOne_MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
MEs = dataOne_net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

timer=read.csv('infiltration_estimation_for_tcga.csv',row.names = 1)
timer=timer[substr(rownames(timer),1,12) %in% substr(rownames(timer),1,12),1:6]
timer=timer[as.numeric(substr(rownames(timer),14,15)) < 10,]#保留肿瘤样本
timer=timer[substr(rownames(timer) ,1,12) %in% substr(rownames(datExprdataOne),1,12),]
traitData = timer
str(traitData)
colnames(traitData)=c("B.cell"      ,           "T.cell.CD4"      ,      "T.cell.CD8"    ,       
                      "Neutrophil"   ,          "Macrophage"    ,         "Myeloid.dendritic.cell")  #  噫嘘唏
sampleName = substr(rownames(datExprdataOne),1,15)
traitData = traitData[match(sampleName, rownames(traitData)),]
write.csv(traitData,'skcm_timer.csv')

nSamples = nrow(datExprdataOne)
modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
rownames(modTraitP)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col),
               font.lab.x = 3,
               font.lab.y = 3,
               cex.lab.y =0.65, 
               cex.lab.x = 0.75,
               ySymbols =  c("brown" ,"red" ,  "blue" , "turquoise", "green", "pink","yellow","black","grey"), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.9, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
load(dataOne_net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
probes = colnames(datExprdataOne)
dimnames(TOM) <- list(probes, probes)
#hub gene 1
module = "blue"
probes = colnames(datExprdataOne)
inModule = (moduleColors==module)
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
corType='spearman'
if (corType=='spearman') {
  geneModuleMembership = as.data.frame(cor(datExprdataOne, MEs_col, use = "p",method = 'spearman'))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExprdataOne, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=='spearman') {
  geneTraitCor = as.data.frame(cor(datExprdataOne, traitData, use = "p",method = 'spearman'))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExprdataOne, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

module = "blue"
pheno = "Myeloid.dendritic.cell"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,col = '#2486b9')+
  abline(a=1,b=1,h=0.5,col='grey')+abline(a=1,b=1,v=0.7,col='grey')
module_gene_matrix=(geneModuleMembership[moduleGenes,])
target_module_gene_matrix=(module_gene_matrix[abs(geneModuleMembership[moduleGenes, module_column])>0.7,])
target_moldule_gene=rownames(target_module_gene_matrix)
trait_gene_matrix=geneTraitCor[moduleGenes, ]
target_trait_gene_matrix=trait_gene_matrix[abs(geneTraitCor[moduleGenes, pheno_column])>0.5,]
target_trait_gene=rownames(target_trait_gene_matrix)
hub_gene_1=intersect(target_moldule_gene,target_trait_gene)
write.csv(hub_gene_1,'WGCNA_hungene_DC_cell.csv')
#hub gene 2
module = "blue"
colnames(traitData)
pheno = "Neutrophil"   
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = '#41ae3c')+
  abline(a=1,b=1,h=0.5,col='gray')+abline(a=1,b=1,v=0.7,col='gray')
module_gene_matrix=(geneModuleMembership[moduleGenes,])
target_module_gene_matrix=(module_gene_matrix[abs(geneModuleMembership[moduleGenes, module_column])>0.7,])
target_moldule_gene=rownames(target_module_gene_matrix)

trait_gene_matrix=geneTraitCor[moduleGenes, ]
target_trait_gene_matrix=trait_gene_matrix[abs(geneTraitCor[moduleGenes, pheno_column])>0.5,]
target_trait_gene=rownames(target_trait_gene_matrix)
hub_gene_2=intersect(target_moldule_gene,target_trait_gene)
#hub gene
hub_gene=c(hub_gene_1,hub_gene_2)
wgcna_hub_gene=unique(hub_gene)
write.csv(wgcna_hub_gene,'WGCNA_GENE.csv')

#cox
rm(list = ls())
options(stringsAsFactors = F)
library(survival)
library(survminer)
library(plyr)
library(dplyr)
library(stringr) 
library(forestplot)
library(lars) 
library(glmnet) 
library(rms)
library(survivalROC)
library(VIM)
library(tidyverse)
library(ggtext)
library(ggsci)
library(ggcorrplot)
library(corrgram)
library(hrbrthemes)
library(wesanderson)
library(limma)
counts=read.csv('counts_tpm.csv',row.names = 1,check.names = FALSE)
clinical=read.csv('TCGA-SKCM-clinical.csv',row.names = 2)
targetgene_1=read.csv('GSEA_result_gene.csv',header = TRUE)
targetgene_1=targetgene_1$x
targetgene_2=read.csv('WGCNA_GENE.csv')
targetgene_2=targetgene_2$x
targetgene=c(targetgene_1,targetgene_2)
targetgene=unique(targetgene)                
exprSet=counts[rownames(counts) %in% targetgene,]
exprSet=t(exprSet) %>% na.omit() %>% as.data.frame()
exprSet$sample=substr(rownames(exprSet),1,12)
eliminate_duplicated_lines=function(x){
  counts_df=x
  counts_df$median=apply(counts_df[,-ncol(counts_df)],1,median)
  counts_df=counts_df[order(counts_df$sample,counts_df$median,decreasing = T),]
  counts_df=counts_df[!duplicated(counts_df$sample),]
  counts_df=counts_df[,-c(ncol(counts_df))]
  x=counts_df
}
exprSet=eliminate_duplicated_lines(exprSet)
rownames(exprSet)=exprSet$sample
exprSet=exprSet[,-ncol(exprSet)]

clinical=clinical[row.names(clinical) %in% rownames(exprSet),]#得到相关样本的临床信息
clinical$days_to_last_follow_up
clinical$days_to_death[is.na(clinical$days_to_death)]=0
clinical$days_to_last_follow_up[is.na(clinical$days_to_last_follow_up)]=0
clinical$days=as.numeric(clinical$days_to_death)+as.numeric(clinical$days_to_last_follow_up)
clinical=clinical %>% 
  dplyr::select(c(vital_status,race,age_at_index,gender,tumor_stage,days,ajcc_pathologic_t,ajcc_pathologic_n,ajcc_pathologic_m))
colnames(clinical)=c('event','race','age','gender','stage',"days",'T','N','M') 
clinical$event=ifelse(clinical$event=='Alive',0,1)   #死亡为1，生存为0
clinical$stage=str_split(clinical$stage,' ',simplify = T)[,2]
clinical$time=clinical$days/30
clinical$names=row.names(clinical)
exprSet$names=rownames(exprSet)

phe_1=inner_join(exprSet,clinical)
rownames(phe_1)=phe_1$names
phe_1=phe_1[,-(ncol(phe_1)-(ncol(clinical))+1)]
write.csv(phe_1,'merge_clinical.csv') 
phe_1=read.csv('merge_clinical.csv',row.names = 1,check.names = FALSE) 
phe_1=phe_1[!phe_1$time ==0,]
GO_cox=function(x){
  FML=as.formula(paste0('Surv(time,event)~phe_1$\'',x,'\''))
  m=coxph(FML,data = phe_1)
  sumcox=summary(m)
  HR=round(exp(coef(m)),2)
  PValue=round(sumcox$coefficients[,5],4)
  CI=paste0(round(sumcox$conf.int[,3:4],2),collapse = '-') 
  CI_MIN=round(sumcox$conf.int[,3],2)
  CI_MAX=round(sumcox$conf.int[,4],2)
  coxresult=data.frame('Characteristics'= x,
                       'Hazaed Ratio'= HR,
                       'CI95'=CI,
                       'P.value'=PValue,
                       'CI_MIN'= CI_MIN,
                       'CI_MMAX'= CI_MAX,check.names = FALSE)
  
}
phe2=colnames(phe_1[,1:(ncol(exprSet)-1)])
ss=lapply(phe2,GO_cox)
ss=ldply(ss)
coxresult=ss[(ss$P.value)<=0.05,]
dim(coxresult)
coxresult=coxresult[!coxresult$CI_MIN==1,]
coxresult=coxresult[!coxresult$CI_MMAX==1,]

#randomSurvivalForest
library(ggRandomForests)
library(randomForestSRC)
library(randomForest)
library(pROC)

rsf_data=read.csv('rsf_data.csv',row.names = 1)
v.out=rfsrc(Surv(time,event)~.,data =rsf_data ,ntree = 1000,forest = TRUE,importance = TRUE,proximity = TRUE,seed = 123)
print(v.out)
plot(v.out)
importance_otu <- v.out$importance
importance_otu <- importance_otu [order(importance_otu,decreasing = FALSE)]
importance_otu=as.data.frame(importance_otu)
importance_otu$gene=rownames(importance_otu)
importance_otu$color=ifelse(importance_otu$importance_otu >=0,"#ee3f4d","#2983bb")
ggplot(importance_otu, aes(x=gene, y=importance_otu)) +
  geom_segment( aes(x=gene, xend=gene, y=0, yend=importance_otu), color="grey") +
  geom_point( color=importance_otu$color, size=1.8) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 60,size =8,face = "bold"),
  )+
  xlab("") +
  ylab("Variable Importance")

vh.out_2=var.select(v.out)
gg_md <- gg_minimal_depth(vh.out_2, lbls = st.labs)
plot(gg_minimal_vimp(gg_md), ) +
  theme(legend.position=c(0.8, 0.2))+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
randomForest_gene=intersect(vh.out_2$topvars,rownames(importance_otu[importance_otu$importance_otu > 0,]))
phe_3=cbind(phe_1[,c("TRAV9-2","KLHDC7B",'CHGA')],phe_1$time,phe_1$event)
m=coxph(Surv(phe_3$`phe_1$time`,phe_3$`phe_1$event`)~.,data = phe_3)
sumcox=summary(m)
b=as.data.frame(sumcox$coefficients)
b[b$`Pr(>|z|)` < 0.05,]

#forest gram
m=summary(m)
p = ifelse(
  m$coefficients[, 5] < 0.001,
  "<0.001 ***",
  ifelse(
    m$coefficients[, 5] < 0.01,
    "<0.01  **",
    ifelse(
      m$coefficients[, 5] < 0.05,
      paste(round(m$coefficients[, 5], 3), " *"),
      round(m$coefficients[, 5], 3)
    )
  )
)
dat2 = as.data.frame(round(m$conf.int[, c(1, 3, 4)], 2))
dat2 = tibble::rownames_to_column(dat2, var = "Trait")
colnames(dat2)[2:4] = c("HR", "lower", "upper")
dat2$HR2 = paste0(dat2[, 2], "(", dat2[, 3], "-", dat2[, 4], ")")
dat2$p = p
ins = function(x) {
  c(x, rep(NA, ncol(dat2) - 1))
}

dat2 = rbind(
  c(NA, NA, NA, NA, NA,NA),
  dat2[1, ],
  dat2[2:nrow(dat2), ]
)
for(i in 2:4) {
  dat2[, i] = as.numeric(dat2[, i])
}
str(dat2)
aaa=dat2[,1]
bbb=c(NA     ,   "TRAV9-2", "KLHDC7B"  , "CHGA" )
dat2[,1]=bbb
forestplot(
  dat2[, c(1, 5, 6)],
  mean = dat2[, 2],
  lower = dat2[, 3],
  upper = dat2[, 4],
  zero = 1,
  boxsize = 0.2,
  col = fpColors(box = '#1075BB', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 2,
  xticks = F,
  is.summary = c(F, F, F),
  align = "l",
  hrzl_lines = list(
    "2" = gpar(lty=1)),
  colgap = unit(5, 'mm')
)


#k-m
Surphe_1=phe_1
Surphe_1$riskscore=
  Surphe_1$KLHDC7B       *   -0.3619641  +
  Surphe_1$CHGA  *   0.1964393
Surphe_1$riskscore_group=ifelse(Surphe_1$riskscore>quantile(Surphe_1$riskscore,0.5),'high','low')
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d","#2983bb"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)

dat =Surphe_1
s=Surv(time, event) ~ KLHDC7B+CHGA+dat$'TRAV9-2'
model <- coxph(s, data = dat )
RiskScore<-predict(model,type = "risk")
names(RiskScore) = rownames(dat)
fp <- RiskScore
phe<-dat
fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
fp_dat$riskgroup= ifelse(fp_dat$fp>=median(fp_dat$fp),'high','low')
sur_dat=data.frame(patientid=1:length(fp),time=phe[names(sort(fp)),'time'],event=phe[names(sort(fp )),'event']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=dat[names(sort(fp)),rownames(b[b$`Pr(>|z|)` < 0.05,])]
exp_dat$KLHDC7B=scale(exp_dat$KLHDC7B)
exp_dat$CHGA=scale(exp_dat$CHGA)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw()+labs(x=NULL,y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)+
  theme(axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 13))

p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  scale_x_discrete(expand = c(0, 0)) +
  labs(x=NULL,y="Survival time(month)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)+
  theme(axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 13))
library(ggplot2)
require(reshape2)
require(scales)
tmp=t(scale(exp_dat))
nba=tmp
colnames(nba)=1:length(fp)
nba.m <- melt(nba) 
colnames(nba.m)=c( "Var1"  ,  "patientid"   ,  "rescale")
nba.m1=nba.m[nba.m$Var1=='CHGA',]
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
p3 <- ggplot(nba.m, aes(x=patientid, y=Var1)) + geom_tile(aes(fill = rescale)) +
  scale_colour_manual(values = mycolors)+
  theme_bw()+
  theme(axis.text.y=element_text(size =13),legend.text = element_text(size = 13),legend.title = element_text(size = 13))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x=NULL,y=NULL)+
  scale_fill_gradient2(low ='#1772b4' , high = '#c04851',mid='white')
library(aplot)
plots=p1 %>% insert_top(p2) %>% insert_top(p3)
plots

#Survival curves of the IRPS in TTN-mutant group
ttn_mutation_samples=read.csv('TTN突变样本.csv',header = FALSE)
ttn_mutation_samples=ttn_mutation_samples$V1     
Surphe_2=phe_1[row.names(phe_1) %in% ttn_mutation_samples,]
Surphe_2$riskscore=
  Surphe_2$KLHDC7B       *   -0.3619641  +
  Surphe_2$CHGA  *   0.1964393
Surphe_2$riskscore_group=ifelse(Surphe_2$riskscore>quantile(Surphe_2$riskscore,0.5),'high','low')
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_2)
ggsurvplot(sfit,palette = c("#ee3f4d","#2983bb"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)
#Survival curves of the IRPS in TTN-wild group
Surphe_3=phe_1[-which(row.names(phe_1) %in% ttn_mutation_samples),]
Surphe_3$riskscore=
  Surphe_3$KLHDC7B       *   -0.3619641  +
  Surphe_3$CHGA  *   0.1964393  
Surphe_3$riskscore_group=ifelse(Surphe_3$riskscore>quantile(Surphe_3$riskscore,0.5),'high','low')
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_3)
ggsurvplot(sfit,palette = c("#ee3f4d","#2983bb"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)

group=read.csv('merge_clinical.csv',row.names = 1,check.names = FALSE)   
group=group[!group$time == 0,]
group$riskscore=
  group$KLHDC7B       *   -0.3619641  +
  group$CHGA  *   0.1964393  
group$riskscore_group=ifelse(group$riskscore>quantile(group$riskscore,0.5),'riskscore-high','riskscore-low')
group$sample=rownames(group)
colnames(group)
group=group[,(ncol(group)-2):ncol(group)]
write.csv(group,'group.csv')

other=read.csv('clinical_more_191sample_1.csv',row.names = 1)
other=other[!duplicated(other$bcr_patient_barcode),]
other=dplyr::rename(other,names=bcr_patient_barcode)
other$breslow_depth_value=factor(other$breslow_depth_value,levels = c( '≤1','1.01~4','＞4.01'),labels = c( '≤1','1.01~4','＞4.01'),order=TRUE)
other$melanoma_clark_level_value=as.factor(other$melanoma_clark_level_value)
other$melanoma_ulceration_indicator=as.factor(other$melanoma_ulceration_indicator)
phe_1=inner_join(exprSet[,c('KLHDC7B','CHGA','names')],clinical) %>% inner_join(other)#合并后的矩阵
rownames(phe_1)=phe_1$names
phe_1=phe_1 %>% select(-(names))
phe_1 $riskscore=
  phe_1 $KLHDC7B       *   -0.3619641  +
  phe_1 $CHGA  *   0.1964393  

nomog1=phe_1
nomog1$gender=factor(nomog1$gender)
nomog1$race=factor(nomog1$race,levels = c('asian','white'),labels = c('asian','white'))
nomog1[nomog1$stage=='iii','stage']='iii'
nomog1[nomog1$stage=='iiia','stage']='iii'
nomog1[nomog1$stage=='iiib','stage']='iii'
nomog1[nomog1$stage=='iiic','stage']='iii'
nomog1[nomog1$stage=='iv','stage']='iv'
nomog1$stage=factor(nomog1$stage,levels = c( 'iii','iv' ),labels = c( 'iii', 'iv' ),ordered = F)
nomog1$T=substr(nomog1$T,1,2)
nomog1[nomog1$T %in% c('TX','T0','T1','T2'),'T']='≤T2'
nomog1$T=factor(nomog1$T,levels = c('≤T2' , 'T3','T4'),labels = c('≤T2', 'T3','T4'),ordered = F)
nomog1$N=substr(nomog1$N,1,2)
nomog1[nomog1$N %in% c('N0','NX'),'N']='Nx~N0'
nomog1[nomog1$N %in% c('N1','N2'),'N']='N1~N2'
nomog1$N=factor(nomog1$N,levels = c('Nx~N0','N1~N2','N3'),labels = c('Nx~N0','N1~N2','N3'),ordered = F)
nomog1$M=substr(nomog1$M,1,2)
nomog1$M=factor(nomog1$M,levels = c('M0',  'M1'),labels = c('M0',  'M1'),ordered = F)
nomog1$riskscore_group=ifelse(nomog1$riskscore>quantile(nomog1$riskscore,0.5),'high','low')
nomog1$riskscore_group=factor(nomog1$riskscore_group,levels = c('low','high'),labels = c('low','high'),ordered = F)
nomog1$TTN=ifelse(rownames(nomog1) %in% ttn_mutation_samples,'TTN_mutation','TTN_wild')
nomog1$TTN=factor(nomog1$TTN,levels = c('TTN_mutation','TTN_wild'),labels = c('TTN_mutation','TTN_wild'))
BRAF_mutation_samples=read.csv('BRAF突变样本.csv')
BRAF_mutation_samples=BRAF_mutation_samples$Tumor_Sample_Barcode
nomog1$BRAF=ifelse(rownames(nomog1) %in% BRAF_mutation_samples,'BRAF_mutation','BRAF_wild')
nomog1$BRAF=factor(nomog1$BRAF,levels = c('BRAF_mutation','BRAF_wild'),labels = c('BRAF_mutation','BRAF_wild'))
write.csv(nomog1,'nomog1.csv')
nomog1=nomog1[!nomog1$time == 0,]
nomog2=nomog1
nomog2$time=as.integer(nomog2$time)
nomog2=data.frame(lapply(nomog2,as.numeric)) 
rownames(nomog2)=rownames(nomog1)
Surphe_1=nomog1
sfit <- survfit(Surv(time, event)~TTN, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d","#2983bb"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)

#Univariate Cox analysis
GO_cox=function(x){
  FML=as.formula(paste0('Surv(time,event)~',x))
  m=coxph(FML,data = nomog2)                           ###
  sumcox=summary(m)
  HR=round(exp(coef(m)),2)
  PValue=round(sumcox$coefficients[,5],3)
  lower=round(sumcox$conf.int[,3],2)
  upper=round(sumcox$conf.int[,4],2)
  CI=paste0(round(sumcox$conf.int[,3:4],2),collapse = '-') 
  coxresult=data.frame('Characteristics'= x,
                       'Hazaed Ratio'= HR,
                       'lower'=lower,
                       'upper'=upper,
                       'CI95'=CI,
                       'P value'=PValue)
  
}
phe_2=c("race","age","gender","stage","T","N","M","breslow_depth_value","melanoma_clark_level_value", 
        "melanoma_ulceration_indicator","riskscore","TTN","BRAF" )
ss=lapply(phe_2,GO_cox)
ss=ldply(ss)#COX结果
coxresult=ss[ss$P.value<0.05,]

#forest gram
p = ifelse(
  ss[, 6] < 0.001,
  "<0.001 ***",
  ifelse(
    ss[, 6]  < 0.01,
    "<0.01  **",
    ifelse(
      ss[, 6]  < 0.05,
      paste(round(ss[, 6] , 3), " *"),
      round(ss[, 6] , 3)
    )
  )
)
ss$P.value=p
ss$CI95 = paste0(ss$Hazaed.Ratio, "(", ss$lower, "-", ss$upper, ")")
dat2=ss
colnames(dat2)=c("Trait" ,"HR"  ,  "lower", "upper" ,"HR2"  , "p" )
ins = function(x) {
  c(x, rep(NA, ncol(dat2) - 1))
}
dat2 = rbind(
  c(NA, NA, NA, NA, NA, NA),
  dat2[1, ],
  dat2[2:nrow(dat2), ]
)
for(i in 2:4) {
  dat2[, i] = as.numeric(dat2[, i])
}
forestplot(
  dat2[, c(1, 5, 6)],
  mean = dat2[, 2],
  lower = dat2[, 3],
  upper = dat2[, 4],
  zero = 1,
  boxsize = 0.2,
  col = fpColors(box = '#1075BB', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 2,
  xticks = F,
  is.summary = c(F, F, F),
  align = "l",
  hrzl_lines = list(
    "2" = gpar(lty=1)),
  colgap = unit(5, 'mm')
)
m=coxph(Surv(time,event)~T+N+age+melanoma_clark_level_value+riskscore+breslow_depth_value,data = nomog2)
m=summary(m)
p = ifelse(
  m$coefficients[, 5] < 0.001,
  "<0.001 ***",
  ifelse(
    m$coefficients[, 5] < 0.01,
    "<0.01  **",
    ifelse(
      m$coefficients[, 5] < 0.05,
      paste(round(m$coefficients[, 5], 3), " *"),
      round(m$coefficients[, 5], 3)
    )
  )
)
dat2 = as.data.frame(round(m$conf.int[, c(1, 3, 4)], 3))
dat2 = tibble::rownames_to_column(dat2, var = "Trait")
colnames(dat2)[2:4] = c("HR", "lower", "upper")
dat2$HR2 = paste0(dat2[, 2], "(", dat2[, 3], "-", dat2[, 4], ")")
dat2$p = p
ins = function(x) {
  c(x, rep(NA, ncol(dat2) - 1))
}
dat2 = rbind(
  c(NA, NA, NA, NA, NA, NA),
  dat2[1, ],
  dat2[2:nrow(dat2), ]
)
for(i in 2:4) {
  dat2[, i] = as.numeric(dat2[, i])
}
str(dat2)
aaa=dat2[,1]
bbb=c( NA,"T","N","age","melanoma_clark_level_value","riskscore","breslow_depth_value")  
dat2[,1]=bbb
forestplot(
  dat2[, c(1, 5, 6)],
  mean = dat2[, 2],
  lower = dat2[, 3],
  upper = dat2[, 4],
  zero = 1,
  boxsize = 0.2,
  col = fpColors(box = '#1075BB', lines = 'black', zero = 'grey'),
  lty.ci = "solid",
  graph.pos = 2,
  xticks = F,
  is.summary = c(F, F, F),
  align = "l",
  hrzl_lines = list(
    "2" = gpar(lty=1)),
  colgap = unit(5, 'mm')
)




nomo<-datadist(nomog1)
options(datadist='nomo')
#multivariate analysis
nomo1 <- cph(Surv(time,event==1)~T+N+riskscore,
             x=T,y=T,
             data=nomog1,
             surv=T,
             time.inc = 12*5)
Cindex <- rcorrcens(Surv(as.numeric(nomog2$time),nomog2$event==1)~predict(nomo1))
Cindex
surv <- Survival(nomo1)
surv1 <- function(x)surv(12*1,lp=x)
surv2 <- function(x)surv(12*3,lp=x)
surv3 <- function(x)surv(12*5,lp=x)
nomo2<-nomogram(nomo1,
                fun=list(surv1,surv2,surv3),
                funlabel=c('1-year LR probability',
                           '3-year LR probability',
                           '5-year LR probability'),
                lp =T, 
                maxscale=100,
                fun.at=c("0.99",'0.8','0.7','0.6','0.5','0.25','0.1')
)
plot(nomo2)
#5years
nomo1 <- cph(Surv(time,event==1)~T+N+riskscore,
             x=T,y=T,
             data=nomog2,
             surv=T,
             time.inc = 12*5)
p<- calibrate(nomo1,
              cmethod='KM',
              method='boot',
              u=12*5,
              m=50, 
              B=100)
plot(p,
     add=F,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.5, 
     lwd=2,
     lty=1,
     errbar.col="#2983bb",
     xlim=c(0.0,1),
     ylim=c(0.0,1),
     xlab="Nomogram-Predicted Probability of 5-year Survival",
     ylab="Actual 5-years Survival",
     col="#ee3f4d", 
     cex.lab=1.3)
abline(0,1,lty=3,lwd=1,col="black")
#3years
nomo1 <- cph(Surv(time,event==1)~T+N+riskscore,
             x=T,y=T,
             data=nomog2,
             surv=T,
             time.inc = 12*3)
p<- calibrate(nomo1,
              cmethod='KM',
              method='boot',#
              u=12*3,
              m=50, 
              B=500)
plot(p,
     add=F,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.5, #
     lwd=2,
     lty=1,
     errbar.col="#2983bb",
     xlim=c(0.0,1),
     ylim=c(0.0,1),
     xlab="Nomogram-Predicted Probability of 3-year Survival",
     ylab="Actual 3-years Survival",
     col="#ee3f4d",
     cex.lab=1.3)
abline(0,1,lty=3,lwd=1,col="black")
#1years
nomo1 <- cph(Surv(time,event==1)~T+N+riskscore,
             x=T,y=T,
             data=nomog2,
             surv=T,
             time.inc = 12*1)
p<- calibrate(nomo1,
              cmethod='KM',
              method='boot',
              u=12*1,
              m=50, 
              B=500)
plot(p,
     add=F,
     conf.int=T,
     subtitles = F,
     cex.subtitles=0.5, 
     lwd=2,
     lty=1,
     errbar.col="#2983bb",
     xlim=c(0.0,1),
     ylim=c(0.0,1),
     xlab="Nomogram-Predicted Probability of 10-year Survival",
     ylab="Actual 10-years Survival",
     col="#ee3f4d",
     cex.lab=1.3)
abline(0,1,lty=3,lwd=1,col="black")

#ROC 
library(survival)
library(survminer)
library(timeROC)
colnames(nomog2)
rt=nomog2 %>% select(T,N,riskscore,stage,gender,age,time,event)
bioCol=rainbow(ncol(rt)-2,0.4)
rt$time=rt$time/12
#绘制
aucText=c()
ROC.T<-timeROC(T=rt$time,
               delta=rt$event,
               marker=rt$T,
               cause=1,weighting="marginal",
               times=quantile(rt$time,probs=seq(0.2,0.8,0.02)))

ROC.N<-timeROC(T=rt$time,
               delta=rt$event,
               marker=rt$N,
               cause=1,weighting="marginal",
               times=quantile(rt$time,probs=seq(0.2,0.8,0.02)))

ROC.risk<-timeROC(T=rt$time,
                  delta=rt$event,
                  marker=rt$riskscore,
                  cause=1,weighting="marginal",
                  times=quantile(rt$time,probs=seq(0.2,0.8,0.02)))

ROC.nomo<-timeROC(T=rt$time,
                  delta=rt$event,
                  marker=rt$riskscore * 0.751839 +rt$T * 1.150215+rt$N* 1.123845,
                  cause=1,weighting="marginal",
                  times=quantile(rt$time,probs=seq(0.2,0.8,0.02)))

ROC.age<-timeROC(T=rt$time,
                 delta=rt$event,
                 marker=rt$age,
                 cause=1,weighting="marginal",
                 times=quantile(rt$time,probs=seq(0.2,0.8,0.02)))

ROC.stage<-timeROC(T=rt$time,
                   delta=rt$event,
                   marker=rt$stage,
                   cause=1,weighting="marginal",
                   times=quantile(rt$time,probs=seq(0.2,0.8,0.02)))

# plot AUC curve for albumin and bilirunbin  with pointwise confidence interval
plotAUCcurve(ROC.T ,conf.int=TRUE,col="#FF9999" )
plotAUCcurve(ROC.N ,conf.int=TRUE,col="#FFF099",add=TRUE )
plotAUCcurve(ROC.risk,conf.int=TRUE,col="#B6FF99" ,add=TRUE)
plotAUCcurve(ROC.age,conf.int=TRUE,col="#99D3FF",add=TRUE)
plotAUCcurve(ROC.stage,conf.int=TRUE,col="#99FFD3",add=TRUE)
plotAUCcurve(ROC.nomo,conf.int=TRUE,col="#B699FF",add=TRUE)
bioCol
legend("bottomright",c("T",'N','riskscore','age','stage','nomogram'),
       col=c( "#FF9999" ,"#FFF099" ,"#B6FF99"  ,"#99D3FF","#99FFD3", "#B699FF"),
       lty=1,lwd=2)
library(survivalROC)
survivalROC_helper <- function(t) {
  survivalROC(Stime=rt$time, status=rt$event, 
              marker = rt$riskscore * 0.8421  +rt$T * 0.5981 +rt$N* 0.6424,
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(1,3,5)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest(cols = c(df_survivalROC)) %>%
  arrange(t, FP, TP)

survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.3f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color=AUC),size = 1,linejoin = 'round')+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +xlab('False positive rate')+ylab('True positive rate')+
  theme(legend.position = c(0.75,0.2),legend.text = element_text(size =14),
        axis.text = element_text(size = 15),axis.title = element_text(size =20))+
  scale_color_manual(values=c( '#20a162',"#ee3f4d", "#2983bb"))

#DCA
library(ggDCA)
T <- cph(Surv(time,event==1)~T,data = nomog1)
N <- cph(Surv(time,event==1)~N,data = nomog1)
IPRS<- cph(Surv(time,event==1)~riskscore,data = nomog1)
age <- cph(Surv(time,event==1)~age,data = nomog1)
stage <- cph(Surv(time,event==1)~stage,data = nomog1)
nomogram <- cph(Surv(time,event==1)~T+N+riskscore,data = nomog1)
data_dca=dca(nomogram,T,N,IPRS,age,stage,times=c(36,60))
ggplot(data_dca,
       color = c('#ed9db2','#2376b7', "#12aa9c",'#70887d','#e2e7bf','#eed045',"#ee3f4d", "#2983bb"))+
theme(axis.title = element_text(size = 22),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 13))


# infiltration estimation 
#estimate
rm(list = ls())
options(stringsAsFactors = F)
library(utils)
rforge <- "http://r-forge.r-project.org"
if(!require("estimate"))install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
library(do)

exprSet=read.csv('counts_tpm.csv',row.names = 1,check.names = FALSE)
colnames(exprSet)=substr(colnames(exprSet),1,12)
dat=log2(edgeR::cpm(exprSet)+1)
library(estimate)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina") 
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='SKCM'
scores=estimate(dat,pro)
estimate_score=as.data.frame(scores)
estimate_score$TumorPurity=cos(0.6049872018+0.0001467884 * scores[,3])
rownames(estimate_score)=Replace(rownames(estimate_score),from = '\\.',to='-')
write.csv(estimate_score,'skcm_estimate_result.csv')

library(tidyverse)
library(dplyr)
nomog2=read.csv('nomog1.csv',row.names = 1)
estimate_score$names=rownames(estimate_score)
nomog2$names=rownames(nomog2)
meta=inner_join(estimate_score,nomog2)
rownames(meta)=meta$names



library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)
scores_2=meta
shapiro.test(scores_2$StromalScore)
shapiro.test(scores_2$ImmuneScore)
shapiro.test(scores_2$ESTIMATEScore)
shapiro.test(scores_2$TumorPurity)

eltext=12
p1=ggplot(scores_2,aes(x=riskscore_group, y=StromalScore,color=riskscore_group)) +
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.2) +
  geom_jitter()+
  ylab('StromalScore') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size=23), 
        axis.title.y = element_text(size=22),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  stat_compare_means(aes(group =riskscore_group ),
                     method = "wilcox.test",label.x.npc = 'middle',size=8,label.y = 2500)

p2=ggplot(scores_2,aes(x=riskscore_group, y=ImmuneScore,color=riskscore_group)) +
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.2) +
  geom_jitter()+
  ylab('ImmuneScore') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size=23), 
        axis.title.y = element_text(size=22),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  stat_compare_means(aes(group =riskscore_group ),
                     method = "wilcox.test",label.x.npc = 'middle',size=8,label.y = 4500)

p3=ggplot(scores_2,aes(x=riskscore_group, y=ESTIMATEScore,color=riskscore_group)) +
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.2) +
  geom_jitter()+
  ylab('ESTIMATEScore') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size=23), 
        axis.title.y = element_text(size=22),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  stat_compare_means(aes(group =riskscore_group ),
                     method = "wilcox.test",label.x.npc = 'middle',size=8,label.y = 6500)

p4=ggplot(scores_2,aes(x=riskscore_group, y=TumorPurity,color=riskscore_group)) +
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.2) +
  geom_jitter()+
  ylab('TumorPurity') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.text = element_text(size = 18),
        axis.title.x = element_text(size=23), 
        axis.title.y = element_text(size=22),
        panel.border = element_blank(),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))+
  stat_compare_means(aes(group =riskscore_group ),
                     method = "wilcox.test",label.x.npc = 'middle',size=8,label.y = 1.2)

p5 <- cowplot::plot_grid(p1,
                         p2,
                         p3,
                         p4,
                         nrow = 2, ncol = 2,labels = c('StromalScore','ImmuneScore','ESTIMATEScore','TumorPurity'))

#cibersort
library(plyr)
library(dplyr)
library(stringr) 
library(pheatmap)
cibersort=read.csv('infiltration_estimation_for_tcga.csv',row.names = 1)[,7:28]
re <- cibersort[substr(rownames(cibersort),1,12) %in% rownames(nomog2),]
re=re[as.numeric(substr(rownames(re),14,15)) < 10,]
colnames(re)=c( "B.cell.naive"             ,        "B.cell.memory"     ,              
                "B.cell.plasma"              ,      "T.cell.CD8."        ,             
                "T.cell.CD4..naive"            ,    "T.cell.CD4..memory.resting"  ,    
                "T.cell.CD4..memory.activated" ,    "T.cell.follicular.helper"    ,    
                "T.cell.regulatory..Tregs."     ,   "T.cell.gamma.delta"    ,          
                "NK.cell.resting"             ,     "NK.cell.activated"  ,             
                "Monocyte"                  ,       "Macrophage.M0"   ,                
                "Macrophage.M1"               ,     "Macrophage.M2"    ,               
                "Myeloid.dendritic.cell.resting" ,  "Myeloid.dendritic.cell.activated",
                "Mast.cell.activated"           ,   "Mast.cell.resting"    ,           
                "Eosinophil"                   ,    "Neutrophil"           )
re$names=substr(rownames(re),1,12)
elimanatr_duplicated_lines=function(x){
  counts_df=x
  counts_df$median=apply(counts_df[,-ncol(counts_df)],1,median)    #按自己需求更改
  counts_df=counts_df[order(counts_df$names,counts_df$median,decreasing = T),]
  counts_df=counts_df[!duplicated(counts_df$names),]
  counts_df=counts_df[,-ncol(counts_df)]
  x=counts_df
}
re=elimanatr_duplicated_lines(re)
re_nomog2=inner_join(re,nomog2)
rownames(re_nomog2)=re_nomog2$names
re_nomog2=re_nomog2[order(re_nomog2$riskscore,decreasing = TRUE),]

annotation_df=re_nomog2 %>%  dplyr::select(riskscore,riskscore_group,T,M,N,TTN,stage)
heatmap_df=as.data.frame(t(re_nomog2[,1:22]))
pheatmap(heatmap_df,scale = "row",cluster_cols =FALSE,cluster_rows = FALSE,
            show_colnames = F,
            annotation_col =annotation_df,
            color = colorRampPalette(c('#1772b4', "white", '#c04851' ))(80))
#Histogram
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set2"))
dat <- re_nomog2[,1:22] %>%  as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
dat$Sample=factor(dat$Sample,levels =rownames(re_nomog2))

ggplot(dat,aes(x=Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity",width = 1) +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.0001,0)) +
  scale_fill_manual(values = mypalette(22))

dat$riskscore_group=ifelse(dat$Sample %in% re_nomog2[re_nomog2$riskscore_group== 'high',]$names,'high','low')

# boxplot
ggplot(dat,aes(Cell_type,Proportion,fill = riskscore_group)) +
  stat_boxplot(geom = 'errorbar',size=0.7,linetype='solid')+
  geom_boxplot(outlier.shape =21,color = "black") + 
  theme_bw() + 
  labs(x = "", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=70,vjust = 0.5,size =12,face = 'bold'),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12,face = 'bold'),
        legend.text = element_text(size = 13),legend.title = element_text(size =13))+
  scale_fill_manual(values =c('#c04851',  '#1772b4'))+ 
  stat_compare_means(aes(group = riskscore_group,label = ..p.signif..),method = "kruskal.test")

library(corrplot)
cor_df=re_nomog2[,1:22]
corr <- round(cor(cor_df,method ="spearman" ), 2)
cor(cor_df) %>% corrplot(method = "circle",order = "hclust",
                         type = "lower",tl.srt = 45,tl.col = "black")
write.csv(re,'skcm_cibersort.csv')

#IPRS associated DEGs
rm(list = ls())
options(stringsAsFactors = F)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(DOSE)
library(stringr)
library(dplyr)
library(org.Hs.eg.db)
library(enrichplot)

counts=read.csv('counts.csv',row.names = 2,check.names = FALSE)
counts_tpm=read.csv('counts_tpm.csv',row.names = 1,check.names = FALSE)
group=read.csv('group.csv',row.names = 1)

counts=counts[,-1]
counts=counts[,-which(duplicated(substr(colnames(counts),1,12)))]

risk_high_sample=rownames(group[group$riskscore_group=='riskscore-high',])
risk_low_sample=rownames(group[group$riskscore_group=='riskscore-low',])

colnames(counts)=substr(colnames(counts),1,12)
mat1 <- counts[,substr(colnames(counts),1,12) %in% risk_high_sample]#risk_score_high的表达矩阵
mat2 <- counts[,substr(colnames(counts),1,12) %in% risk_low_sample]#risk_low的表达矩阵

gene <- cbind(mat1,mat2)
gene$sum=apply(gene,1,sum)
gene=gene[gene$sum != 0,]
colnames(gene[,(ncol(gene)-1):ncol(gene)])
gene=gene[,-ncol(gene)]
gene_1=apply(gene, 1,as.integer)
rownames(gene_1)=colnames(gene)
gene=as.data.frame(t(gene_1))
rm(gene_1)
coldata <- data.frame(row.names = colnames(gene),group = factor(c(rep('riskscore_high',ncol(mat1)),rep('riskscore_low',ncol(mat2)))))

dds <- DESeqDataSetFromMatrix(countData = gene, colData = coldata, design = ~group)

dds <- DESeq(dds, parallel = FALSE)  
suppressMessages(dds)
res <- results(dds, contrast = c('group',  'riskscore_high','riskscore_low'), pAdjustMethod = 'fdr', alpha = 0.05)
deseq_res <- as.data.frame(res[order(res$padj), ])
deseq_res$gene_id <- rownames(deseq_res)
write.table(deseq_res[c(7, 1:6)], 'pingfen_DEGs_DESeq2.txt', row.names = FALSE, sep = '\t', quote = FALSE)
deseq_res <- read.delim('pingfen_DEGs_DESeq2.txt', sep = '\t')
DEGs=deseq_res[deseq_res$padj < 0.05,]
DEGs=DEGs[order(DEGs$log2FoldChange),]
DEG_analyse=DEGs[abs(DEGs$log2FoldChange) >= 2,] %>% na.omit()
pingfen_DEGs_genelist=unique(DEG_analyse$gene_id)
pingfen_DEGs_counts=counts_tpm[rownames(counts_tpm) %in% pingfen_DEGs_genelist,] %>% t() %>%  as.data.frame()
pingfen_DEGs_counts$names=substr(rownames(pingfen_DEGs_counts),1,12)
group$names=rownames(group)
group_pingfen_counts=inner_join(group,pingfen_DEGs_counts,by='names')
rownames(group_pingfen_counts)=group_pingfen_counts$names
counts_analyze=group_pingfen_counts[,c(1,(ncol(group)+1):ncol(group_pingfen_counts))]         #711个免疫相关基因  36
#Correlation coefficient calculation 
data=as.data.frame(t(counts_analyze))
gene_name1<-c()
gene_name2<-c()
cor_r<-c()
pvalue<-c()
for (i in 1:nrow(data)){
  r=1
  g1=rownames(data)[i]
  g2=rownames(data)[r]
  c_r=cor(as.numeric(data[i,]),as.numeric(data[r,]),method="spearman")
  p=cor.test(as.numeric(data[i,]),as.numeric(data[r,]),method ="spearman")[[3]]
  gene_name1=c(gene_name1,g1)
  gene_name2=c(gene_name2,g2)
  cor_r=c(cor_r,c_r)
  pvalue=c(pvalue,p)
  
}

data_cor<-data.frame(gene_name1,gene_name2,cor_r,pvalue)
head(data_cor)
data_cor=data_cor[data_cor$pvalue < 0.05,]
data_cor=data_cor[abs(data_cor$cor_r)>0.7,]
data_cor=data_cor[-1,]

genelist=unique(data_cor$gene_name1)
is.vector(genelist)
#GO analysis
deg_GO=enrichGO(genelist,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL'
                ,ont = 'BP',pvalueCutoff = 0.01)
write.csv(deg_GO@result[deg_GO@result$p.adjust < 0.05,],'pingfen_GO.csv')
pingfen_GO_degs=read.csv('pingfen_GO_genelist.csv',header = FALSE)
pingfen_GO_degs=pingfen_GO_degs$V1
genelist=genelist[genelist %in% pingfen_GO_degs]
deg_GO=enrichGO(genelist,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL'
                ,ont = 'BP',pvalueCutoff = 0.01)
write.csv(deg_GO@result[deg_GO@result$p.adjust < 0.05,],'pingfen_GO.csv')
library(ggplot2)
GO_select=read.csv('pingfen_GO_select.csv',row.names = 1)
GO_select=GO_select[1:12,]$Description
p1=dotplot(deg_GO, showCategory=GO_select) + scale_y_discrete(labels = function(x) str_wrap(x, width =40))
ggplot(p1$data,aes(x = GeneRatio,y = Description))+
  geom_point(aes(color =p.adjust,
                 size = Count))+
  scale_color_gradient(name='p.adjust',low = "red", high = "blue",limit=c(0,0.051),breaks=c(0,0.01,0.05))+
  xlab("")+ylab('')+
  theme_bw()+ 
  theme(axis.text = element_text(size = 15),
        legend.text = element_text(size =13.5),
        legend.title = element_text(size = 13.5))+
  scale_y_discrete(labels = function(x) str_wrap(x, width =40))

#kegg analysis
gene_kegg=bitr(genelist,
               fromType = 'SYMBOL',
               toType = c('ENTREZID'),
               OrgDb = 'org.Hs.eg.db')
kk<-enrichKEGG(gene = gene_kegg$ENTREZID ,#需要为vector格式
               organism = 'hsa',
               keyType = "kegg",          
               pAdjustMethod = "BH",
               pvalueCutoff = 0.01)
kk2=setReadable(kk,
                OrgDb = 'org.Hs.eg.db',
                keyType ='ENTREZID' )
write.csv(kk@result[kk@result$p.adjust < 0.05,],'pingfen_KEGG.csv')
kegg_select=read.csv('pingfen_KEGG_select.csv')
kegg_select=kegg_select$Description
p=enrichplot::cnetplot(kk2,circular=TRUE,showCategory = kegg_select,node_label = "gene",colorEdge = TRUE,node_label_size=0.5)
p+ggraph::scale_edge_color_discrete(labels = function(x) str_wrap(x, 40) )


#heatmap
gene_select=read.csv('pingfen_GO_select_gene.csv',header = FALSE)
genelist=gene_select$V1
nomog2=read.csv('nomog1.csv',row.names = 1)
genelist_counts=counts_tpm[rownames(counts_tpm) %in% genelist,] %>% t() %>% as.data.frame()
genelist_counts$names=substr(rownames(genelist_counts),1,12)
nomog2$names=rownames(nomog2)
merge_df= inner_join(nomog2,genelist_counts)
rownames(merge_df)=merge_df$names
merge_df=merge_df[order(merge_df$riskscore,decreasing = TRUE),]
annotation_df2=merge_df %>% dplyr::select(TTN,riskscore,riskscore_group)
heatmap_counts=merge_df[,(ncol(nomog2)+1):ncol(merge_df)] %>% t() %>% as.data.frame()
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
pheatmap(heatmap_counts,scale = "row",cluster_cols =FALSE,
         show_colnames = F,
         annotation_col =annotation_df2,
         color = colorRampPalette(c('#1772b4', "white", '#c04851' ))(80))


# Correlation analysis between IPRS and various variables
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(qqboxplot)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggstatsplot)
library(fmsb)

nomog2=read.csv('nomog1.csv',row.names = 1)
nomog2$names=rownames(nomog2)
estimate=read.csv('skcm_estimate_result.csv',row.names = 1)
estimate$names=rownames(estimate)
cibersort=read.csv('skcm_cibersort.csv',row.names = 1)
timer=read.csv('skcm_timer.csv',row.names = 1)
colnames(timer)=c("B.cell.timer" ,"T.cell.CD4.timer" ,"T.cell.CD8.timer" , 
                  "Neutrophil.timer","Macrophage.timer" ,"Myeloid.dendritic.cell.timer")
timer$names=substr(rownames(timer),1,12)
tmb=read.csv('skcm_tmb.csv',row.names = 2)
tmb$names=rownames(tmb)
counts=read.csv('counts_tpm.csv',row.names = 1,check.names = FALSE)
target_gene=c('CTLA4','PDCD1','CD274','LAG3' , 'HAVCR2'  ,'TIGIT' ,'VSIR','CD4','CD8A','CD8B','TTN') 
checkpoint_gene=c('CTLA4','PDCD1','CD274','LAG3' , 'HAVCR2'  ,'TIGIT' ,'VSIR','CD4','CD8A','CD8B') 
counts_target_gene=counts %>% t() %>% as.data.frame() %>% dplyr::select(all_of(target_gene))
counts_target_gene=dplyr::rename(counts_target_gene,'TTN_exp'='TTN')
counts_target_gene$names=substr(rownames(counts_target_gene),1,12)



Merge_func <- function(x,y){
  
  df <- merge(x, y, by = "names",all=T)
  rownames(df) <- df$names
  return(df)
}
meta <- Reduce(Merge_func,list(nomog2,estimate,cibersort,timer,tmb,counts_target_gene))

ggscatterstats(
  data = meta,
  title = "StromalScore",
  x = "riskscore",
  y ="StromalScore",
  type = "nonparametric",            
  results.subtitle = TRUE,
  point.args =(aes(color="#ee3f4d",size=2,alpha=0.6)),
  xlab = "riskscore",        
  ylab = "StromalScore", 
  line.color = "red", 
  ggtheme =theme_bw(base_size = 14)+theme(axis.title = element_text(size = 18),axis.text = element_text(size = 14)),
  ggstatsplot.layer = FALSE, 
  smooth.line.args = aes(color="#2983bb",alpha=0.6,method='lm'),
  marginal.type = "density", 
  xfill = "#ee3f4d", 
  yfill = "#2983bb",
  xalpha = 0.6, 
  yalpha = 0.6, 
  centrality.para = "median", 
  messages = FALSE 
)


ggscatterstats(
  data = meta,
  title = "ImmuneScore",
  x = "riskscore",
  y ="ImmuneScore",
  type = "nonparametric",            
  results.subtitle = TRUE,
  point.args =(aes(color="#ee3f4d",size=2,alpha=0.6)),
  xlab = "riskscore",          
  ylab = "ImmuneScore",            
  line.color = "red", 
  ggtheme =theme_bw(base_size = 14)+theme(axis.title = element_text(size = 18),axis.text = element_text(size = 14)),
  ggstatsplot.layer = FALSE, 
  smooth.line.args = aes(color="#2983bb",alpha=0.6,method='lm'),
  marginal.type = "density", 
  xfill = "#ee3f4d", 
  yfill = "#2983bb",
  xalpha = 0.6, 
  yalpha = 0.6, 
  centrality.para = "median", 
  messages = FALSE 
)

ggscatterstats(
  data = meta,
  title = "ESTIMATEScore",
  x = "riskscore",
  y ="ESTIMATEScore",
  type = "nonparametric",            
  results.subtitle = TRUE,
  point.args =(aes(color="#ee3f4d",size=2,alpha=0.6)),
  xlab = "riskscore",         
  ylab = "ESTIMATEScore",             
  line.color = "red", 
  ggtheme =theme_bw(base_size = 14)+theme(axis.title = element_text(size = 18),axis.text = element_text(size = 14)),
  ggstatsplot.layer = FALSE, 
  smooth.line.args = aes(color="#2983bb",alpha=0.6,method='lm'),
  marginal.type = "density", 
  xfill = "#ee3f4d", 
  yfill = "#2983bb", 
  xalpha = 0.6, 
  yalpha = 0.6, 
  centrality.para = "median", 
  messages = FALSE 
)

ggscatterstats(
  data = meta,
  title = 'TumorPurity',
  x = "riskscore",
  y ="TumorPurity",
  type = "nonparametric",   
  results.subtitle = TRUE,
  point.args =(aes(color="#ee3f4d",size=2,alpha=0.6)),
  xlab = "riskscore",         
  ylab = "TumorPurity",       
  line.color = "red", 
  ggtheme =theme_bw(base_size = 14)+theme(axis.title = element_text(size = 18),axis.text = element_text(size = 14)),
  ggstatsplot.layer = FALSE, 
  smooth.line.args = aes(color="#2983bb",alpha=0.6,method='lm'),
  marginal.type = "density",
  xfill = "#ee3f4d", 
  yfill = "#2983bb", 
  xalpha = 0.6,
  yalpha = 0.6, 
  centrality.para = "median",
  messages = FALSE
)

eltext <- 23
library(RColorBrewer)
shapiro.test(meta$total_perMB_log)
ggplot(meta,aes(x=as.factor(TTN), y=total_perMB_log,color=TTN)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('TMB') + xlab('')+
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.text.x  = element_text(size=eltext),
        axis.title.y = element_text(size=eltext),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(aes(group =TTN),method = "wilcox.test",label.y = 3,size=5)


#IPRS & T stage
b=meta[-which(is.na(meta$T)),]
mycompare=list(c( '≤T2', 'T3'),c('T3' ,'T4' ),c('≤T2' , 'T4'))
ggplot(b,aes(x=T, y=riskscore,color=T)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('riskscore') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 15),
        panel.border = element_blank(), panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(comparisons = mycompare,label.y = c(1.1,1.35,1.6),size=4)+
  stat_compare_means(label.y = 2.1,size=5)
library(ggplotify)
library(survival)
library(survminer)
Surphe_1=b[b$T=='≤T2',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for ≤T2 Stage ",
           ncensor.plot = FALSE)

Surphe_1=b[b$T=='T3',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8, 
           title="Survival curve for T3 Stage ",
           ncensor.plot = FALSE)

Surphe_1=b[b$T=='T4',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months",
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for T4 Stage ",
           ncensor.plot = FALSE)

#IPRS & N stage
b=meta[-which(is.na(meta$N)),]
table(b$N)
mycompare=list(c( 'Nx~N0', 'N1~N2'),c('N1~N2' ,'N3' ),c('Nx~N0' , 'N3'))
ggplot(b,aes(x=N, y=riskscore,color=N)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('riskscore') + #ggtitle('a) boxplot') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 15),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(comparisons = mycompare,label.y = c(1.1,1.5,1.8),size=4)+
  stat_compare_means(label.y = 2.5,size=5)

Surphe_1=b[b$N=='Nx~N0',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for Nx~N0 Stage ",
           ncensor.plot = FALSE)

Surphe_1=b[b$N=='N1~N2',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for N1~N2 Stage ",
           ncensor.plot = FALSE)

Surphe_1=b[b$N=='N3',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for N3 Stage ", 
           ncensor.plot = FALSE)

#IPRS & M stage
b=meta[-which(is.na(meta$M)),]
ggplot(b,aes(x=M, y=riskscore,color=M)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('riskscore') + #ggtitle('a) boxplot') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 15),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(aes(group = M),method = "wilcox.test",label.y = 1.2,size=5)

Surphe_1=b[b$M=='M0',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for M0 Stage ",
           ncensor.plot = FALSE)


Surphe_1=b[b$M=='M1',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,pval.coord=c(50,0.2),
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for M1 Stage ", 
           ncensor.plot = FALSE)

#IPRS & Age
colnames(meta)
b=meta[-which(is.na(meta$age)),]
b$age=ifelse(as.numeric(b$age)<=65 ,'≤65',">65")
ggplot(b,aes(x=age, y=riskscore,color=age)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('riskscore') + #ggtitle('a) boxplot') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 15),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(aes(group = age),method = "wilcox.test",label.y = 1.1,size=5)

Surphe_1=b[b$age=='>65',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for >65 subgroup ", 
           ncensor.plot = FALSE)


Surphe_1=b[b$age=='≤65',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for ≤65 subgroup ", 
           ncensor.plot = FALSE)

#IPRS & Gender
colnames(meta)
b=meta
ggplot(b,aes(x=gender, y=riskscore,color=gender)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('riskscore') + #ggtitle('a) boxplot') +
  xlab('Sex')+
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 15),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  stat_compare_means(aes(group = gender),method = "wilcox.test",label.y = 1.3,size=5)
#K-M
Surphe_1=b[b$gender=='female',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for female ", 
           ncensor.plot = FALSE)

Surphe_1=b[b$gender=='male',]
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d", "#2983bb"),
           risk.table =FALSE,pval =TRUE,pval.coord=c(1,0.15),
           conf.int =FALSE,xlab ="Time in months", 
           ggtheme =theme(plot.title=element_text(size=eltext-2, face="plain", hjust=0.5), 
                          axis.title.x = element_text(size=eltext-3), 
                          axis.title.y = element_text(size=eltext-3),
                          axis.text = element_text(size = 15),
                          panel.border = element_rect(fill = 'transparent',colour = 'grey55'),
                          panel.background = element_rect(fill="transparent"),
                          panel.grid.major = element_line(colour = "grey70"),
                          panel.grid.minor = element_line(colour="grey80"),
                          legend.text = element_text(size = 15)), pval.size=8,
           title="Survival curve for male ", 
           ncensor.plot = FALSE)

#TMB & TTN mutation
colnames(meta)
b=meta
b$TTN=ifelse(b$TTN=='TTN_mutation','TTN-MUT','TTN-WT')
ggplot(b,aes(x=TTN, y=total_perMB_log,color=TTN)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('log(TMB)')+xlab('') + #ggtitle('a) boxplot') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 19),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14))+
  stat_compare_means(aes(group = TTN),method = "wilcox.test",label.y = c(3.1),size=6)

#TTN mutation and TTN expression
colnames(meta)
b=meta
b$TTN=ifelse(b$TTN=='TTN_mutation','TTN-MUT','TTN-WT')
ggplot(b,aes(x=TTN, y=TTN_exp,color=TTN)) +scale_color_brewer(palette = 'Set1')+
  geom_violin(trim = FALSE)+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.4) +
  geom_jitter(width = 0.15)+
  ylab('TTN expression')+xlab('') + #ggtitle('a) boxplot') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = 19),
        panel.border = element_blank(), 
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 14))+
  stat_compare_means(aes(group = TTN),method = "wilcox.test",label.y = c(2.1),size=6)



#IPRS &immune gene
data=meta %>% dplyr::select(all_of(c('riskscore',checkpoint_gene)))
data=as.data.frame(t(data))
gene_name1<-c()
gene_name2<-c()
cor_r<-c()
pvalue<-c()
for (i in 1:nrow(data)){
  r=1
  g1=rownames(data)[i]
  g2=rownames(data)[r]
  c_r=cor(as.numeric(data[i,]),as.numeric(data[r,]),method="spearman")
  p=cor.test(as.numeric(data[i,]),as.numeric(data[r,]),method ="spearman")[[3]]
  gene_name1=c(gene_name1,g1)
  gene_name2=c(gene_name2,g2)
  cor_r=c(cor_r,c_r)
  pvalue=c(pvalue,p)
  
}

data_cor<-data.frame(gene_name1,gene_name2,cor_r,pvalue)
head(data_cor)
data_cor=data_cor[-1,]
rownames(data_cor)=data_cor$gene_name1
rownames(data_cor)
rownames(data_cor)=c("CTLA4","PD-1","PD-L1","LAG-3","TIM-3",
                     "TIGIT","VISTA","CD4","CD8A","CD8B")
#pd-1--PDCD1,PD-L1---CD274,LAG-3--LAG3,TIM-3--HAVCR2,VISTA--VSIR
rt=data_cor[,c(3,4)]
rt$name <- rownames(rt)
rt <- rt[order(rt$cor_r), ]
rt$name <-factor(rt$name, levels = rt$name)
rt$scale_p=(rt$pvalue-min(rt$pvalue))/(max(rt$pvalue)-min(rt$pvalue))
colnames(rt)=c("correlation"  , "pvalue",  "name"   , "scale_p")
ggplot(rt,
       aes(y = name, x = correlation))+
  
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#ee3f4d","#2983bb"),limit=c(0,0.05),breaks=c(0,0.01,0.05))+
  geom_segment(aes(x = 0, y = name, xend = correlation, yend =name), size = 1, color = "grey")+
  geom_point(stat = 'identity', aes(color=pvalue,size =abs(correlation)))+
  theme_classic(base_size = 20)+
  scale_y_discrete(position = "right")+
  ylab('') + xlab('correlation')+xlim(-0.85,0)

eltext <- 18

shapiro.test(meta$PDCD1)
ggplot(meta,aes(x=riskscore_group, y=PDCD1,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('PD-1 expression') +  
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,9)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.x.npc = 'middle',size=6)

shapiro.test(meta$CD274)
ggplot(meta,aes(x=riskscore_group, y=CD274,color=riskscore_group))+
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('PD-L1 expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,8)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.x.npc = 'middle',size=6)


shapiro.test(meta$CTLA4)
ggplot(meta,aes(x=riskscore_group, y=CTLA4,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('CTLA-4 expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,8)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.x.npc = 'middle',size=6)


shapiro.test(meta$VSIR)
ggplot(meta,aes(x=riskscore_group, y=VSIR,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('VISTA expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,8)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.x.npc = 'middle',size=6)



shapiro.test(meta$TIGIT)
ggplot(meta,aes(x=riskscore_group, y=TIGIT,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('TIGIT expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,8)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.x.npc = 'middle',size=6)


shapiro.test(meta$LAG3)
ggplot(meta,aes(x=riskscore_group, y=LAG3,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('LAG-3 expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,10)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.y = 10,label.x.npc = 'middle',size=6)

shapiro.test(meta$HAVCR2)
ggplot(meta,aes(x=riskscore_group, y=HAVCR2,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('TIM-3 expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,10)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.y = 10,label.x.npc = 'middle',size=6)


shapiro.test(meta$CD8B)
ggplot(meta,aes(x=riskscore_group, y=CD8B,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('CD8B expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,10)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.y = 10,label.x.npc = 'middle',size=6)


shapiro.test(meta$CD8B)
ggplot(meta,aes(x=riskscore_group, y=CD4,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('CD4 expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,10)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.x.npc = 'middle',label.y =10,size=6)

shapiro.test(meta$CD8B)
ggplot(meta,aes(x=riskscore_group, y=CD8A,color=riskscore_group)) +
  scale_color_brewer(palette = 'Set1')+
  geom_boxplot(varwidth = FALSE, notch =FALSE,size=1,width=0.3) +
  geom_jitter(width = 0.03)+
  ylab('CD8A expression') + 
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), 
        axis.title.x = element_text(size=eltext), 
        axis.title.y = element_text(size=eltext),
        axis.text = element_text(size = eltext-1),
        panel.border  = element_blank(),panel.background = element_rect(fill="white"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  ylim(0,10)+
  stat_compare_means(aes(group = riskscore_group),method = "wilcox.test",label.y = 10,label.x.npc = 'middle',size=6)


#Survival curves of 6 immune cell infiltration scores from the timer2.0 database
library(survival)
library(survminer)
colnames(meta)
#$B.cell     
meta$B.cell_TIMER_group=ifelse(meta$B.cell.timer>median(meta$B.cell.timer),'high','low')
sfit <- survfit(Surv(time, event)~meta$B.cell_TIMER_group, data=meta)
ggsurvplot(sfit,
           title='B cell',xlab ="Time in months",
           palette = c("#ee3f4d","#2983bb"),
           risk.table =FALSE,ncensor.plot = FALSE,
           pval =TRUE,pval.size=6.5,conf.int =TRUE,
           ggtheme =theme_light()+
             theme(axis.text = element_text(size = 14),
                   axis.title.y = element_text(size =16),
                   axis.title.x = element_text(size = 17),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size =12),
                   plot.title = element_text(size = 20,hjust = 0.5)))

#T.cell.CD4   
meta$T.cell.CD4._TIMER_group=ifelse(meta$T.cell.CD4.timer>median(meta$T.cell.CD4.timer),'high','low')
sfit <- survfit(Surv(time, event)~ meta$T.cell.CD4._TIMER_group, data=meta)
ggsurvplot(sfit,
           title='CD4 T cell',xlab ="Time in months",
           palette = c("#ee3f4d","#2983bb"),
           risk.table =FALSE,ncensor.plot = FALSE,
           pval =TRUE,pval.size=6.5,conf.int =TRUE,
           ggtheme =theme_light()+
             theme(axis.text = element_text(size = 14),
                   axis.title.y = element_text(size =16),
                   axis.title.x = element_text(size = 17),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size =12),
                   plot.title = element_text(size = 20,hjust = 0.5)))


#T.cell.CD8
meta$T.cell.CD8._TIMER_group=ifelse(meta$T.cell.CD8.timer>median(meta$T.cell.CD8.timer),'high','low')
sfit <- survfit(Surv(time, event)~ meta$T.cell.CD8._TIMER_group, data=meta)
ggsurvplot(sfit,
           title='CD8 T cell',xlab ="Time in months",
           palette = c("#ee3f4d","#2983bb"),
           risk.table =FALSE,ncensor.plot = FALSE,
           pval =TRUE,pval.size=6.5,conf.int =TRUE,
           ggtheme =theme_light()+
             theme(axis.text = element_text(size = 14),
                   axis.title.y = element_text(size =16),
                   axis.title.x = element_text(size = 17),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size =12),
                   plot.title = element_text(size = 20,hjust = 0.5)))


#Neutrophil  
meta$Neutrophil_TIMER_group=ifelse( meta$Neutrophil.timer>median( meta$Neutrophil.timer),'high','low')
sfit <- survfit(Surv(time, event)~  meta$Neutrophil_TIMER_group, data=meta)
ggsurvplot(sfit,
           title='Neutrophil',xlab ="Time in months",
           palette = c("#ee3f4d","#2983bb"),
           risk.table =FALSE,ncensor.plot = FALSE,
           pval =TRUE,pval.size=6.5,conf.int =TRUE,
           ggtheme =theme_light()+
             theme(axis.text = element_text(size = 14),
                   axis.title.y = element_text(size =16),
                   axis.title.x = element_text(size = 17),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size =12),
                   plot.title = element_text(size = 20,hjust = 0.5)))

#Macrophage     
meta$Macrophage_TIMER_group=ifelse( meta$Macrophage.timer>median( meta$Macrophage.timer),'high','low')
sfit <- survfit(Surv(time, event)~  meta$Macrophage_TIMER_group, data=meta)
ggsurvplot(sfit,
           title='Macrophage',xlab ="Time in months",
           palette = c("#ee3f4d","#2983bb"),
           risk.table =FALSE,ncensor.plot = FALSE,
           pval =TRUE,pval.size=6.5,conf.int =TRUE,
           ggtheme =theme_light()+
             theme(axis.text = element_text(size = 14),
                   axis.title.y = element_text(size =16),
                   axis.title.x = element_text(size = 17),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size =12),
                   plot.title = element_text(size = 20,hjust = 0.5)))


#DC   
meta$Myeloid.dendritic.cell_TIMER_group=ifelse( meta$Myeloid.dendritic.cell.timer >median( meta$Myeloid.dendritic.cell.timer),'high','low')
sfit <- survfit(Surv(time, event)~ meta$Myeloid.dendritic.cell_TIMER_group, data=meta)
ggsurvplot(sfit,
           title='Myeloid dendritic cell',xlab ="Time in months",
           palette = c("#ee3f4d","#2983bb"),
           risk.table =FALSE,ncensor.plot = FALSE,
           pval =TRUE,pval.size=6.5,conf.int =TRUE,
           ggtheme =theme_light()+
             theme(axis.text = element_text(size = 14),
                   axis.title.y = element_text(size =16),
                   axis.title.x = element_text(size = 17),
                   legend.title = element_text(size = 13),
                   legend.text = element_text(size =12),
                   plot.title = element_text(size = 20,hjust = 0.5)))

#IPS & IPRS
IPS=read.csv('IPS.csv',row.names = 1,check.names = FALSE)
IPS$names=rownames(IPS)
dat=IPS %>% dplyr::inner_join(meta,by='names') %>% 
  dplyr::select(riskscore_group,ips_ctla4_neg_pd1_neg,ips_ctla4_neg_pd1_pos,ips_ctla4_pos_pd1_neg,ips_ctla4_pos_pd1_pos)
library(reshape2)
library(ggbeeswarm)
colnames(dat)=c("riskscore_group", "CTLA4(-) PD1(-)", "CTLA4(-) PD1(+)" ,"CTLA4(+) PD1(-)","CTLA4(+) PD1(+)")
dat1=melt(dat,id.vars = c("riskscore_group" ))
ggplot(dat1, aes(x=variable, y=value, color=riskscore_group)) +
  stat_boxplot(geom = 'errorbar',size=1,width=0.5,linetype='solid')+
  geom_boxplot(varwidth = TRUE, notch = TRUE,size=1,width=0.5)+
  ylab('IPS') +xlab('')+ #ggtitle('a) boxplot') +
  theme(plot.title=element_text(size=eltext, face="plain", hjust=0.5), axis.title.x = element_text(size=12), axis.title.y = element_text(size=20),
        panel.border = element_blank(), panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour="grey80"),
        axis.line =element_line(colour = 'black'),
        legend.key.size = unit(0.4,'inches'),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text.x = element_text(size = 13,face = 'bold'),
        axis.text.y = element_text(size = 15))+
  scale_color_manual(values =c("#ee3f4d", "#2983bb"))+ 
  stat_compare_means(aes(group = riskscore_group,label = ..p.signif..),method = "kruskal.test")


#Validation of IPRS with samples in GEO
rm(list = ls())
options(stringsAsFactors = F)
library(survival)
library(survminer)
library(plyr)
library(dplyr)
library(stringr) 
library(rms)
library(survivalROC)
library(VIM)
library(tidyverse)
GEO_counts=read.csv('GEO_counts_log2_tpm_normalize.csv',row.names = 1,check.names = FALSE)
geo_clinical=read.csv('GEO_clinical.csv',row.names = 1)
GEO_counts=GEO_counts[,colnames(GEO_counts) %in% rownames(geo_clinical)]#获得84个GEO样本
counts=GEO_counts %>% t() %>% as.data.frame()
counts$names=rownames(counts)
geo_clinical$names=rownames(geo_clinical)
phe_1=merge(counts,geo_clinical,by='names',all = FALSE)
phe_1$riskscore=phe_1$KLHDC7B*-0.3619641+phe_1$CHGA * 0.1964393
phe_1=phe_1 %>% dplyr::select(KLHDC7B,CHGA,riskscore,time,event)

Surphe_1=phe_1
Surphe_1$riskscore_group=ifelse(Surphe_1$riskscore>quantile(Surphe_1$riskscore,0.5),'high','low')
sfit <- survfit(Surv(time, event)~riskscore_group, data=Surphe_1)
ggsurvplot(sfit,palette = c("#ee3f4d","#2983bb"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = TRUE)

dat =Surphe_1
s=Surv(time, event) ~ riskscore
model <- coxph(s, data = dat )
RiskScore<-predict(model,type = "risk")
names(RiskScore) = rownames(dat)
fp <- RiskScore
phe<-dat
fp_dat=data.frame(patientid=1:length(fp),fp=as.numeric(sort(fp)))
fp_dat$riskgroup= ifelse(fp_dat$fp>=median(fp_dat$fp),'high','low')
sur_dat=data.frame(patientid=1:length(fp),time=phe[names(sort(fp)),'time'],event=phe[names(sort(fp )),'event']) 
sur_dat$event=ifelse(sur_dat$event==0,'alive','death')
sur_dat$event=factor(sur_dat$event,levels = c("death","alive"))
exp_dat=dat[names(sort(fp)),c('KLHDC7B','CHGA')]
exp_dat$KLHDC7B=scale(exp_dat$KLHDC7B)
exp_dat$CHGA=scale(exp_dat$CHGA)
p1=ggplot(fp_dat,aes(x=patientid,y=fp))+geom_point(aes(color=riskgroup))+
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw()+labs(x=NULL,y="Risk score")+
  geom_hline(yintercept=median(fp_dat$fp),colour="black", linetype="dotted",size=0.8)+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)+
  theme(axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 13))

p2=ggplot(sur_dat,aes(x=patientid,y=time))+geom_point(aes(col=event))+theme_bw()+
  scale_colour_manual(values = c("#ee3f4d","#2983bb"))+
  scale_x_discrete(expand = c(0, 0)) +
  labs(x=NULL,y="Survival time(month)")+
  geom_vline(xintercept=sum(fp_dat$riskgroup=="low"),colour="black", linetype="dotted",size=0.8)+
  theme(axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.y = element_text(size = 13))
library(ggplot2)
require(reshape2)
require(scales)
tmp=t(scale(exp_dat))
nba=tmp
colnames(nba)=1:length(fp)
nba.m <- melt(nba) 
colnames(nba.m)=c( "Var1"  ,  "patientid"   ,  "rescale")
nba.m1=nba.m[nba.m$Var1=='CHGA',]
mycolors <- colorRampPalette(c("white", "green", "red"), bias = 1.2)(100)
p3 <- ggplot(nba.m, aes(x=patientid, y=Var1)) + geom_tile(aes(fill = rescale)) +
  scale_colour_manual(values = mycolors)+
  theme_bw()+
  theme(axis.text.y=element_text(size =13),legend.text = element_text(size = 13),legend.title = element_text(size = 13))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x=NULL,y=NULL)+
  scale_fill_gradient2(low ='#1772b4' , high = '#c04851',mid='white')
library(aplot)
plots=p1 %>% insert_top(p2) %>% insert_top(p3)
plots

library(survivalROC)
rt=phe_1
rt$time=rt$time/12
survivalROC_helper <- function(t) {
  survivalROC(Stime=rt$time, status=rt$event, 
              marker = rt$riskscore ,
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(1,3,5)) %>%
  mutate(survivalROC = map(t, survivalROC_helper),
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         df_survivalROC = map(survivalROC, function(obj) {
           as_data_frame(obj[c("cut.values","TP","FP")])
         })) %>%
  dplyr::select(-survivalROC) %>%
  unnest(cols = c(df_survivalROC)) %>%
  arrange(t, FP, TP)

survivalROC_data1 <- survivalROC_data %>% 
  mutate(auc =sprintf("%.3f",auc))%>% 
  unite(year, t,auc,sep = " year AUC: ")

AUC =factor(survivalROC_data1$year)

survivalROC_data1 %>%
  ggplot(mapping = aes(x = FP, y = TP)) +
  geom_path(aes(color=AUC),size = 1,linejoin = 'round')+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  theme_bw() +xlab('False positive rate')+ylab('True positive rate')+
  theme(legend.position = c(0.75,0.2),legend.text = element_text(size =14),
        axis.text = element_text(size = 15),axis.title = element_text(size =20))+
  scale_color_manual(values=c( '#20a162',"#ee3f4d", "#2983bb"))

#drug analysis
rm(list = ls())
options(stringsAsFactors = F)
library(reshape2)
library(ggpubr)
library(oncoPredict)
library(oncoPredict)
library(data.table)
library(gtools)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res)

testExpr<- GDSC2_Expr[,sample(1:ncol(GDSC2_Expr),10)]
colnames(testExpr)=paste0('test',colnames(testExpr))
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]

