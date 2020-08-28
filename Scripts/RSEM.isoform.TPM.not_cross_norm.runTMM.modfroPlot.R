library(edgeR)
library(tidyverse)

rnaseqMatrix = read.table("RSEM.isoform.TPM.not_cross_norm", header=T, row.names=1, com='', check.names=F)
rnaseqMatrix = as.matrix(rnaseqMatrix)
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study)
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
write.table(exp_study$samples, file="RSEM.isoform.TPM.not_cross_norm.TMM_info.txt", quote=F, sep="\t", row.names=F)

Samples <- exp_study$samples


Samples <- Samples %>%
  mutate(Tissue=ifelse(grepl("HM",group),"HM",
                       ifelse(grepl("GUT",group),"GUT",
                       "OV"))) 

DotPlot <- ggplot(Samples,aes(x=Tissue,y=lib.size,fill=Tissue))+
  geom_dotplot(binaxis='y', 
               stackdir='center')+
  stat_summary(fun.data=mean_se, 
               fun.args = list(mult=1), 
               geom="pointrange",
               color="black",
               size=0.7,
               show.legend = FALSE)+
  stat_summary(fun.data = mean_se,
               geom = "errorbar", 
               size=0.5,
               show.legend = FALSE)+
  theme_classic()+
  scale_fill_brewer(palette="Dark2")+
  scale_y_continuous(limits=c(999950,1000050),labels = function(x) format(x, scientific = TRUE))
ggsave(DotPlot,file="Lib.size.pdf")
ggsave(DotPlot,file="Lib.size.eps")

BoxPlot <- ggplot(Samples,aes(x=Tissue,y=eff.lib.size,fill=Tissue))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_brewer(palette="Dark2")
  
Stats <- Samples %>% 
  group_by(Tissue) %>% 
  summarise(avg=mean(lib.size),SD=sd(lib.size))


