### tissue-spf modification site enrichment analysis
library(ggplot2)
library(tidyr)
library(org.Hs.eg.db)
library(ggsci)
library(clusterProfiler,lib.loc = "/home/zhaoying/R_lib/")
library(ggsci)
library(ggthemes,lib.loc = "/home/zhaoying/R_lib/")
library(enrichplot,lib.loc = "/home/zhaoying/R_lib/")
library(org.Hs.eg.db)
library(scales)
setwd("/home/zhaoying/workbai/NC/Figure3/tissue-spf/")
tissues <- c("brain","heart","liver","lung","muscle","kidney","spleen")
m6A <- read.table("./m6A_tissue_spf_site.tsv",sep="\t",header = T)
head(m6A)
for (tissue in tissues) {
  ego <- enrichGO(m6A[m6A$tissue==tissue,]$geneID,keyType = 'SYMBOL',OrgDb = org.Hs.eg.db,ont = "BP")
  #  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  p1 <- dotplot(ego) + scale_color_material(palette = c("blue"),reverse = TRUE)
  ggsave(p1,filename = paste0("m6A",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("m6A",tissue,"_ego.csv"))
}

m5C <- read.table("./m5C_tissue_spf_site.tsv",sep="\t",header = T)
head(m5C)
for (tissue in tissues) {
  ego <- enrichGO(m5C[m5C$tissue==tissue,]$geneID,keyType = 'SYMBOL',OrgDb = org.Hs.eg.db,ont = "BP")
  #  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  p1 <- dotplot(ego) + scale_color_material(palette = c("blue"),reverse = TRUE)
  ggsave(p1,filename = paste0("m5C",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("m5C",tissue,"_ego.csv"))
}

psU <- read.table("./psU_tissue_spf_site.tsv",sep="\t",header = T)
head(psU)
for (tissue in tissues) {
  ego <- enrichGO(psU[psU$tissue==tissue,]$geneID,keyType = 'SYMBOL',OrgDb = org.Hs.eg.db,ont = "BP")
  #  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  p1 <- dotplot(ego) + scale_color_material(palette = c("blue"),reverse = TRUE)
  ggsave(p1,filename = paste0("psU",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("psU",tissue,"_ego.csv"))
}

