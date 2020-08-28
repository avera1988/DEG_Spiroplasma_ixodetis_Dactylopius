#!/usr/bin/Rscript
library(DESeq2)
library(tximport)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(BiocParallel)

print("Starting DESeq Analysis")
dir <- "./"
system("ls -l |grep rsem|awk '{print $9}' > samples.id.txt")
run <- readLines("samples.id.txt")
files <- file.path(dir, run)
names(files) <- run
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
colnames(txi.rsem$counts)<- factor(sub(".spiro_rsem.dir.genes.results", "", colnames(txi.rsem$counts)))
colnames(txi.rsem$abundance)<- factor(sub(".spiro_rsem.dir.genes.results", "", colnames(txi.rsem$abundance)))
colnames(txi.rsem$length)<- factor(sub(".spiro_rsem.dir.genes.results", "", colnames(txi.rsem$length)))
head(txi.rsem$counts)
write.table(txi.rsem$counts,"spiro.transcripts.rawCounts.tab",quote=F,sep="\t")
txi.rsem$length[txi.rsem$length == 0 ] <- 1
sampleTable <- data.frame(condition=factor(sub("\\.","",sub("\\d", "",colnames(txi.rsem$counts)))))
rownames(sampleTable) <- colnames(txi.rsem$counts)
colnames(sampleTable) <- "treatment"
dds <- DESeqDataSetFromTximport(txi.rsem, colData=sampleTable,  design = ~ treatment)
dds <- DESeq(dds, fitType = "local")
conditionsName <- unique(as.vector(dds$treatment))
#Obtaining fpkm
FPKM <- fpkm(dds)
FpkmDF <- as.data.frame(FPKM)
save(FpkmDF, file="fpkm.RData")
#Obtainign Normalized counts DESEQ
table_counts_normalized <- counts(dds, normalized=TRUE)
save(table_counts_normalized, file="DDS.Normalized.counts.RData")
#COmparisons GUT vs Others
for (i in seq(2,3)){
  print (c(conditionsName[1],"vs",conditionsName[(i)]))
  res05 <- results(dds, contrast=c("treatment",conditionsName[1],conditionsName[(i)]),alpha=0.05)
  res05Ordered <- res05[order(res05$pvalue),]
  data05 <- na.omit(as.data.frame(res05Ordered))
  #writing all results 
  write.table(res05Ordered, paste0("DESEQ.",conditionsName[1],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F)
  print(summary(res05))
  #Getting differential expressed genes below p-adjust value 0.05
  keepAllDEgenes <- (data05$padj<=0.05)
  genesDE<-data05[keepAllDEgenes,]
  print (c("Number of all DE genes:", dim(genesDE)[1]))
  write.table(genesDE, paste0("All.DE.genes.",conditionsName[1],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F)
  #Obtainig genes with logFoldchage 1.5
  DF05 <- as.data.frame(res05)
  DF05$Gene_id <- row.names(res05)
  DF05 <- DF05[,c(7,1:6)]
  genesDELFC1.5 <- dplyr::filter(DF05, padj <= 0.05 & abs(log2FoldChange) >= 0.585)
  print (c("Number of DE genes with LFC_1.5:", dim(genesDELFC1.5)[1]))
  write.table(genesDELFC1.5, paste0("All.DE.genes.LFC1.5.",conditionsName[1],"_vs_",conditionsName[(i)],".tab"), sep="\t", quote=F,row.names=F)
  #Generating MA plots
  pdf(paste0("MA_plot.",conditionsName[1],"_vs_",conditionsName[(i)],".pdf"))	
  plotMA(results(dds, contrast=c("treatment",conditionsName[1],conditionsName[(i)]),alpha = 0.05), main=paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
  abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)	
  dev.off()
  #Generating Volcano plots with Enhanced Volcano (https://bioconductor.org/packages/devel/bioc/html/EnhancedVolcano.html)
  volcan <- EnhancedVolcano(res05,
			    lab = "",
			    x = "log2FoldChange",
			    y = "pvalue", 
			    FCcutoff = 0.5,
			    pCutoff=0.05,
			    transcriptPointSize = 1.5,
			    transcriptLabSize = 3.0,
                            title = paste0("DE Genes from:", conditionsName[1],"_vs_",conditionsName[(i)]))
  ggsave(volcan, file=paste0("volcano.plot.",conditionsName[1],"_vs_",conditionsName[(i)],".png"))
}

#Comparisson HEM vs OV
print (c(conditionsName[2],"vs",conditionsName[3]))
res05 <- results(dds, contrast=c("treatment",conditionsName[2],conditionsName[3]),alpha=0.05)
res05Ordered <- res05[order(res05$pvalue),]
data05 <- na.omit(as.data.frame(res05Ordered))
#writing all results 
write.table(res05Ordered, paste0("DESEQ.",conditionsName[2],"_vs_",conditionsName[3],".tab"), sep="\t", quote=F)
print(summary(res05))
#Getting differential expressed genes below p-adjust value 0.05
keepAllDEgenes <- (data05$padj<=0.05)
genesDE<-data05[keepAllDEgenes,]
print (c("Number of all DE genes:", dim(genesDE)[2]))
write.table(genesDE, paste0("All.DE.genes.",conditionsName[2],"_vs_",conditionsName[3],".tab"), sep="\t", quote=F)
#Obtainig genes with logFoldchage 1.5
DF05 <- as.data.frame(res05)
DF05$Gene_id <- row.names(res05)
DF05 <- DF05[,c(7,1:6)]
genesDELFC1.5 <- dplyr::filter(DF05, padj <= 0.05 & abs(log2FoldChange) >= 0.585)
print (c("Number of DE genes with LFC_1.5:", dim(genesDELFC1.5)[2]))
write.table(genesDELFC1.5, paste0("All.DE.genes.LFC1.5.",conditionsName[2],"_vs_",conditionsName[3],".tab"), sep="\t", quote=F,row.names=F)
#Generating MA plots
pdf(paste0("MA_plot.",conditionsName[2],"_vs_",conditionsName[3],".pdf"))	
plotMA(results(dds, contrast=c("treatment",conditionsName[2],conditionsName[3]),alpha = 0.05), main=paste0("DE Genes from:", conditionsName[2],"_vs_",conditionsName[3]))
abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)	
dev.off()
#Generating Volcano plots with Enhanced Volcano (https://bioconductor.org/packages/devel/bioc/html/EnhancedVolcano.html)
volcan <- EnhancedVolcano(res05,
			  lab = "",
			  x = "log2FoldChange",
			  y = "pvalue", 
			  FCcutoff = 0.5,
			  pCutoff = 0.05,
			  transcriptPointSize = 1.5,
			  transcriptLabSize = 3.0,
                          title = paste0("DE Genes from:", conditionsName[2],"_vs_",conditionsName[3]))
ggsave(volcan, file=paste0("volcano.plot.",conditionsName[2],"_vs_",conditionsName[3],".png"))

#for PCA using VST form dds object generated above
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color=treatment,shape=treatment, label=pcaData$name)) + 
  geom_point(size=3.5)+ 
  geom_text(size=3)+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste("PC2:",percentVar[2],"% variance"))+
  coord_fixed()+
  scale_y_continuous(limits=c(-5,10))+
  scale_x_continuous(limits=c(-10,10))+ 
  theme_bw()
ggsave(PCA_plot, file="PCA.pdf")
