library(dplyr)
library(tidyr)
library(tibble)
#Midgut DESEQ Table
DESEQGUTHM <- read.delim(file = "DESEQ.GUT_vs_HM.tab", 
                       fill= T, h = T, sep = "\t")
DESEQGUTOV <- read.delim(file = "DESEQ.GUT_vs_OV.tab", 
                       fill= T, h = T, sep = "\t")
DESEQHMOV <-read.delim(file = "DESEQ.HM_vs_OV.tab", 
                      fill= T, h = T, sep = "\t")
DESEQGUTHM$Gene_id <- row.names(DESEQGUTHM)
DESEQGUTOV$Gene_id <- row.names(DESEQGUTOV)
DESEQHMOV$Gene_id <- row.names(DESEQHMOV)
DEseqTotal <- dplyr::full_join(DESEQGUTHM[,c(2,5,6,7)],DESEQGUTOV[,c(2,5,6,7)], by = "Gene_id",
                           suffix = c(".GUTHM", ".GUTOV")) %>%
  dplyr::full_join(DESEQHMOV[,c(2,5,6,7)],by = "Gene_id") %>%
  dplyr::rename(padj.HMOV = padj) %>% 
  dplyr::rename(log2FoldChange.HMOV = log2FoldChange)%>%
  dplyr::rename(pvalue.HMOV=pvalue) %>%
  dplyr::mutate_at(vars(log2FoldChange.GUTHM,log2FoldChange.GUTOV,log2FoldChange.HMOV),
                   funs(round(., 2)))
DEseqSpiroTotal <- DEseqTotal[,c(4,1,2,3,5,6,7,8,9,10)]
save(DEseqSpiroTotal,file="DESEQ.Spiro.total.RData")
