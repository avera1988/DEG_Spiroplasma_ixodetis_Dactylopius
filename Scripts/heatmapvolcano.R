######################################################################################################
#This script produces heatmaps, volcanoplots and VennDiagrams of DEG from S. ixodetis DCF in Gut,
#  Hemolymph and ovary
#Author Arturo Vera
######################################################################################################
library(EnhancedVolcano)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(stringr)
library(pheatmap)
library(functional)
library(gplots)
library(RColorBrewer)
library(patchwork)
library(systemPipeR)
#Reading data
load("DESEQ.Spiro.total.RData")
Annot <- openxlsx::read.xlsx("DCF.annotation.tab.xlsx")
Annot <- Annot %>% 
  dplyr::rename("Gene_id"=locusID)
DEseqSpiroTotal <- DEseqSpiroTotal %>% 
  full_join(Annot,by="Gene_id")

#Obtain DEG from DESEQ2 tables
GUTHM <- dplyr::select(DEseqSpiroTotal,matches("Gene_id|GUTHM"))
GUTHMDEG <- GUTHM %>% dplyr::filter(abs(log2FoldChange.GUTHM) >= 0.5 & padj.GUTHM <= 0.05) %>% 
  arrange(desc(log2FoldChange.GUTHM))

GUTOV <- dplyr::select(DEseqSpiroTotal,matches("Gene_id|GUTOV"))
GUTOVDEG <- GUTOV %>% dplyr::filter(abs(log2FoldChange.GUTOV) >= 0.5 & padj.GUTOV <= 0.05) %>% 
  arrange(desc(log2FoldChange.GUTOV))

HMOV <- dplyr::select(DEseqSpiroTotal,matches("Gene_id|HMOV"))
HMOVDEG <- HMOV %>% dplyr::filter(abs(log2FoldChange.HMOV) >= 0.5 & padj.HMOV <= 0.05) %>% 
  arrange(desc(log2FoldChange.HMOV))

DEGTotLFC <- dplyr::full_join(GUTHMDEG,GUTOVDEG,by="Gene_id") %>%
  dplyr::full_join(HMOVDEG,by="Gene_id") %>%
  dplyr::select(matches("Gene_id|log2FoldChange|padj"))

#merge DESEQ and annotation
DEGTotLFCAnnot <- dplyr::inner_join(DEGTotLFC,Annot,by="Gene_id")

DEGTotLFCAnnot <- DEGTotLFCAnnot %>%
  mutate_if(is.factor,as.character) %>%
  mutate(annotation=if_else(grepl("hyp",Prokka_Annotation) & !grepl("^-",COG_Annotation),
                      COG_Annotation,Prokka_Annotation)) %>%
  mutate(annotation=if_else(grepl("hypot",annotation),"hyp",annotation)) %>%
  arrange(Gene_id)

#Create a matrix for pheatmap
DEGForMatrix <- tidyr::unite(DEGTotLFCAnnot,
                             "Gene_id_Annot",
                             c("Gene_id","annotation","GeneName"),sep=";")
DEGForMatrix <- DEGForMatrix %>%
  select(matches("Gene_id_Annot|log2"))
rownames(DEGForMatrix) <- DEGForMatrix$Gene_id_Annot
DEGForMatrix <- DEGForMatrix[,2:4]

ph_annot <- data.frame(colnames(DEGForMatrix))
names(ph_annot) <- "names"
ph_annot$Tissue <- ifelse(grepl("GUTHM",ph_annot$names), "GUTvsHM",
                          ifelse(grepl("GUTOV",ph_annot$names),"GUTvsOV", "HMvsOV"))
row.names(ph_annot) <- ph_annot$names
ph_annot <- dplyr::select(ph_annot,Tissue)
annot_colors=list(Tissue = c(GUTvsHM= wesanderson::wes_palettes$Darjeeling1[2],
                               GUTvsOV=wesanderson::wes_palettes$Darjeeling1[3],
                               HMvsOV=wesanderson::wes_palettes$Darjeeling1[5]))  

Colors <- brewer.pal(10, "PRGn")

PHmapBacteria <- pheatmap(DEGForMatrix,
                          show_rownames = T,
                          fontsize_row=10,
                          cellwidth=18,
                          color =Colors
                          ,treeheight_row = 0,
                          cluster_rows = F,
                          treeheight_col = 0,
                          annotation_col = ph_annot,
                          cluster_cols=F,
                          show_colnames = T,
                          na_col = "Gray",
                          annotation_colors = annot_colors,
                          breaks=rev(c(2.5,2,1.5,1,0.5,-0.5,-1,-1.5,-2,-2.5)),
                          legend_breaks =rev(c(2.5,2,1.5,1,0.5,-0.5,-1,-1.5,-2,-2.5)),
                          legend_labels=rev(c(2.5,2,1.5,1,0.5,-0.5,-1,-1.5,-2,-2.5)))
ggsave(PHmapBacteria,filename = "Spiro.DEG.heatmap.mod.pdf",width = 12,height = 12)

pdf("scale.pdf")
plot(NULL, xlim=c(0,length(Colors)), ylim=c(0,1), 
     xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(Colors)-1), 0, 1:length(Colors), 1, col=Colors)
dev.off()

DEGNoHyp <- DEGTotLFCAnnot %>%
  filter(!grepl("hyp",annotation))

DEGHyp <- DEGTotLFCAnnot %>%
  filter(grepl("hyp",annotation))


##Volcano Plots##


GUTHMAnnot <- GUTHM %>%
  inner_join(Annot) %>%
  mutate_if(is.factor,as.character) %>%
  mutate(annot=if_else(grepl("hyp",Prokka_Annotation) & !grepl("^-",COG_Annotation),
                            COG_Annotation,Prokka_Annotation)) %>%
  mutate(annot=if_else(grepl("hypot",annot),"hyp",annot)) %>%
  select(matches("Gene_id|log|padj|annot|GeneName",ignore.case = FALSE)) %>%
  mutate(GeneLabels=if_else(!grepl("^-",GeneName),
                            GeneName,annot))


GHgeneLab <- GUTHMAnnot %>% 
  dplyr::filter(abs(log2FoldChange.GUTHM) >= 0.5 & padj.GUTHM <= 0.05) %>% 
  filter(!grepl("hyp",GeneLabels)) 
  

GHvolc <- EnhancedVolcano(as.data.frame(GUTHMAnnot), 
                x='log2FoldChange.GUTHM',
                y='padj.GUTHM',
                lab=GUTHMAnnot$GeneLabels,
                selectLab = GHgeneLab$GeneLabels,
                pointSize = 2.0, 
                FCcutoff= 0.5, 
                pCutoff = 0.05,
               labSize = 3.0,
                xlim = c(-5,5),
                ylim = c(0,8),
               drawConnectors = TRUE,
                title = "GUTHM",
               subtitle = '',
               caption = paste0('Total of Differential Exp. Genes(FC < 1.5 padj < 0.05) = ', 
                                nrow(GUTHMDEG)))


GUTOVAnnot <- GUTOV %>%
  inner_join(Annot) %>%
  mutate_if(is.factor,as.character) %>%
  mutate(annot=if_else(grepl("hyp",Prokka_Annotation) & !grepl("^-",COG_Annotation),
                       COG_Annotation,Prokka_Annotation)) %>%
  mutate(annot=if_else(grepl("hypot",annot),"hyp",annot)) %>%
  select(matches("Gene_id|log|padj|annot|GeneName",ignore.case = FALSE)) %>%
  mutate(GeneLabels=if_else(!grepl("^-",GeneName),
                            GeneName,annot))


GOgeneLab <- GUTOVAnnot %>% 
  dplyr::filter(abs(log2FoldChange.GUTOV) >= 0.5 & padj.GUTOV <= 0.05) %>% 
  filter(!grepl("hyp",GeneLabels)) 



GOvolc <- EnhancedVolcano(as.data.frame(GUTOVAnnot), 
                          x='log2FoldChange.GUTOV',
                          y='padj.GUTOV',
                          lab=GUTOVAnnot$GeneLabels,
                          selectLab = GOgeneLab$GeneLabels,
                          pointSize = 2.0, 
                          FCcutoff= 0.5, 
                          pCutoff = 0.05,
                          labSize = 3.0,
                          xlim = c(-5,5),
                          ylim = c(0,12),
                          drawConnectors = TRUE,
                          title = "GUTOV",
                          subtitle = '',
                          caption = paste0('Total of Differential Exp. Genes(FC < 1.5 padj < 0.05) = ', 
                                           nrow(GUTOVDEG)))

HMOVAnnot <- HMOV %>%
  inner_join(Annot) %>%
  mutate_if(is.factor,as.character) %>%
  mutate(annot=if_else(grepl("hyp",Prokka_Annotation) & !grepl("^-",COG_Annotation),
                       COG_Annotation,Prokka_Annotation)) %>%
  mutate(annot=if_else(grepl("hypot",annot),"hyp",annot)) %>%
  select(matches("Gene_id|log|padj|annot|GeneName",ignore.case = FALSE)) %>%
  mutate(GeneLabels=if_else(!grepl("^-",GeneName),
                            GeneName,annot))
HOgeneLab <- HMOVAnnot %>% 
  dplyr::filter(abs(log2FoldChange.HMOV) >= 0.5 & padj.HMOV <= 0.05) %>% 
  filter(!grepl("hyp",GeneLabels)) 

HOvolc <- EnhancedVolcano(as.data.frame(HMOVAnnot), 
                          x='log2FoldChange.HMOV',
                          y='padj.HMOV',
                          lab=HMOVAnnot$GeneLabels,
                          selectLab = HOgeneLab$GeneLabels,
                          pointSize = 2.0, 
                          FCcutoff= 0.5, 
                          pCutoff = 0.05,
                          labSize = 3.0,
                          xlim = c(-5,5),
                          ylim = c(0,12),
                          drawConnectors = TRUE,
                          title = "HMOV",
                          subtitle = '',
                          caption = paste0('Total of Differential Exp. Genes(FC < 1.5 padj < 0.05) = ', 
                                           nrow(HMOVDEG)))

TotalVolcano <- GHvolc+GOvolc+HOvolc #This works pretty well for pdf but not for eps using patchwork lib.

TptV <- grid.arrange(GHvolc,GOvolc,HOvolc,ncol=3) #Works for eps to illustrator 

ggsave(TotalVolcano,file="DEG.volcanos.mod.pdf",width = 20,height = 12,device = cairo_pdf)
ggsave(TptV,file="DEG.volcanos.mod.eps",width = 18,height = 10,device = cairo_ps)

#Generating a Total table with DESEq values and annotations
load("../DDS.Normalized.counts.RData")
DDS <- as.data.frame(table_counts_normalized)
rm(table_counts_normalized)

DDS$Gene_id <- rownames(DDS)

DDS <- DDS %>%
  select(Gene_id,everything())

DEseqSpiroTotal <- DEseqSpiroTotal %>%
  arrange(Gene_id)
DEseqSpiroTotal <- full_join(DEseqSpiroTotal, DDS, by="Gene_id")

list.DEdata <- list(DEGTotLFCAnnot,DEseqSpiroTotal)
names(list.DEdata) <- c("DEG","Expr.Values.w.Annotation")

#Saving in excel format

openxlsx::write.xlsx(list.DEdata,file="DEG.Exp.val.Annot.xlsx")


##Venn Diagram

#obtain Gene Id total DEG and list it
GH <- GUTHMDEG %>%
  select(Gene_id)
GO <- GUTOVDEG %>%
  select(Gene_id)
HO <- HMOVDEG %>%
  select(Gene_id)
list.data <- list(GH,GO, HO)
names(list.data) <- c("GutvsHM","GutvsOV","HMvsOV")

#Use a vector list to overlap
listaDaVector <- lapply(list.data,function(x){as.vector(x$Gene_id)})

#GEnerate a VennDiagram

VEnnset <- overLapper(listaDaVector,type="vennsets")
#plot Venn Diagramm
VennDiagrama <- vennPlot(VEnnset,lcol=c(annot_colors$Tissue[1],
                                        annot_colors$Tissue[2],
                                        annot_colors$Tissue[3]),
                         lines=c(annot_colors$Tissue[1],
                                         annot_colors$Tissue[2],
                                         annot_colors$Tissue[3]))
pdf("VennDiagramDEG.pdf")
vennPlot(VEnnset,lcol=c(annot_colors$Tissue[1],
                        annot_colors$Tissue[2],
                        annot_colors$Tissue[3]),
         lines=c(annot_colors$Tissue[1],
                 annot_colors$Tissue[2],
                 annot_colors$Tissue[3]))
dev.off()

#to extract information from the Diagramm 
Names <- names(VEnnset@vennlist)
GetDataVenset <- function(NUM,COMP){
  W <- as.data.frame(VEnnset@vennlist[NUM])
 # write.table(W, paste0(COMP,".tab"),sep="\t",quote=FALSE,row.names=FALSE)
  W$Gene_id <- W[,1]
  W <- W %>%
    inner_join(DEGTotLFCAnnot,by="Gene_id") %>%
    select(matches("Gene_id|log2|annot",ignore.case = FALSE))
  return(W)
}

GetDataVenset(4,"GutHMvsGutHM") %>% filter(!grepl("hyp",annotation))

GUTHMDEGAnnot <- GUTHMDEG %>%
  inner_join(Annot) %>%
  arrange(Gene_id)
GUTOVDEGAnnot <- GUTOVDEG %>%
  inner_join(Annot) %>%
  arrange(Gene_id)
HMOVDEGAnnot <- HMOVDEG %>%
  inner_join(Annot) %>%
  arrange(Gene_id)
  
