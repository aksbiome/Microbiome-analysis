################ Deseq2 ##########################
library(ggplot2)
library(ggplot2)
library(dplyr)
library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(DESeq2)
library(ggpubr)
library(phyloseq)
library(tidyverse)
library(ggpubr)


###
setwd("C:/Users/Akshay/PS_object_with tree/Trial")
pseq3 <- readRDS("pseq3.RDS")
pseq3

#add metadata
metadata <- read.table("C:/Users/Akshay/PS_object_with tree/Trial/Glutathione_metadat_phyloseq.csv", sep = ",", header = T, row.names = 1, check.names = F)
metadata$sampleID <- rownames(metadata)

#linking phyloseq
sample_data(pseq3) <- metadata
pseq3

#pseq1 <- subset_samples(pseq5, Diagnosis != "H-Nonsmoker")
#pseq2 <- subset_samples(pseq1, Diagnosis!= "H-Smoker")
#pseq2 <- subset_samples(pseq3, Diagnosis != "DG_Alpha")
pseq4 <- subset_samples(pseq3, Diagnosis != "D_Gamma")
pseq <- prune_taxa(taxa_sums(pseq4) >1, pseq4)
pseq
sample_data(pseq)



#perform DESeq2
diagdds = phyloseq_to_deseq2(pseq, ~ Diagnosis)
diagdds = DESeq(diagdds, test="Wald", fitType="local")


res = results(diagdds, cooksCutoff = FALSE)

alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq)[rownames(sigtab), ], "matrix"))

sigtab$Diagnosis <- ifelse(sigtab$log2FoldChange>0,'D_Alpha','DG_Gamma')
head(sigtab)

dim(sigtab)

write.csv(sigtab, "sigtab_D_Alpha_DG_Gamma_20_5_22_asv.csv")

#plot

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Paired", ...) {
  scale_fill_brewer(palette = palname, ...)
}

x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))


#Plot_2
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color= Diagnosis)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold") ,strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))

#Plot_3
ggplot(sigtab, aes(y=Phylum, x=log2FoldChange, color= Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_boxplot(size=1) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+ theme(axis.text.x = element_text(angle = 90, size = 15, colour = "black", vjust = 0.2, hjust = 0.5, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold") ,strip.text = element_text(size = 18, colour = "black", face = "bold"), legend.title = element_text(size = 20, colour = "black", face = "bold"), legend.text = element_text(size = 22, face = "bold", colour = "black"), axis.text.y = element_text(size = 16, colour = "black", face = "bold"))


#write.csv(sigtab, "DG_A_DG_G_Time_point.csv")

