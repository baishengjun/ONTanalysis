library(data.table,lib.loc = "/home/zhaoying/miniconda3/lib/R/library")
library(methylKit)
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
# library(methylKit)


setwd("/home/zhaoying/workbai/NC/07.m5C/m5Cnet_res/methylKit")

tissues <- c("brain","heart","liver","lung","muscle","kidney","spleen")
com_tissue <- combn(tissues,2)

com_tissue

# 加载必要的包
library(stringr)

# 获取当前目录下的文件列表
files <- list.files()

# 定义组织列表
tissues <- c("brain", "heart", "liver", "lung", "muscle", "kidney", "spleen")

# 生成两两不重复的组合
combinations <- combn(tissues, 2)

# 转置组合矩阵以便每行表示一个组合
combinations <- t(combinations)

# 初始化一个列表来存储结果
results <- list()

# 遍历每一对组合
for (i in 1:nrow(combinations)) {
  m <- combinations[i, 1]
  n <- combinations[i, 2]
  
  # 过滤以 m 结尾的文件
  f1 <- files[str_detect(files, paste0( m, "$"))]
  
  # 过滤以 n 结尾的文件
  f2 <- files[str_detect(files, paste0( n, "$"))]
  
  # 合并文件列表
  file_ls <- as.list(c(f1, f2))
  
  # 生成对比向量
  contrast <- c(rep(1, length(f1)), rep(0, length(f2)))
  myobj <- methRead(file_ls,sample.id=file_ls,treatment=contrast,mincov=5,assembly='hg38')
  meth=methylKit::unite(myobj, destrand=FALSE,min.per.group=1L)
  myDiff=calculateDiffMeth(meth,mc.cores=20,adjust = 'qvalue')
  write.table(myDiff,file=paste0(m,"_",n,".tsv"),row.names=FALSE,quote=F,sep="\t")
  # 将结果存储在列表中
  results[[i]] <- list(file_list = file_ls, contrast = contrast)
}

tissues[c(2:7)]
setwd("/home/zhaoying/workbai/NC/Figure6/NAR")
for (tissue in tissues[c(2:7)]) {
  d1 <- read.table(paste0(tissue,"_ol.tsv"),sep="\t",header = T)
  ego <- enrichGO(unique(d1$geneID),OrgDb = org.Hs.eg.db,ont = "BP",keyType = "SYMBOL")
  write.table(ego@result,file = paste0(tissue,"_ego.tsv"),
              sep="\t",row.names = F,quote = F)
  p2 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  ggsave(p2,filename = paste0("/home/zhaoying/workbai/NC/Figure6/NAR/",tissue,"_ego.pdf"))
  # break
}
head(d1)
# read.table(paste0())
#














file_ls <- list()
for (tissue in tissues) {
  fs <- files[str_detect(files, paste0( "brain", "$"))]
  # fs <- as.list(fs)
  file_ls <- append(file_ls,fs)
}
file_ls

data <- read.table("../tissue_diff.tsv",sep="\t",header = T)
head(data)
library(clusterProfiler)
# data <-  data$geneID
library(org.Hs.eg.db)
ego <- enrichGO(data$geneID,keyType = "SYMBOL",OrgDb = org.Hs.eg.db)
ego2 <- simplify(ego)
dotplot(ego2)





data <- read.table("../data.csv",header = T,sep=",")
head(data)
library(clusterProfiler)
# data <-  data$geneID
library(org.Hs.eg.db)
ego <- enrichGO(data$geneID,keyType = "SYMBOL",OrgDb = org.Hs.eg.db)
ego2 <- simplify(ego)
dotplot(ego2)


## enrichment analysis of tissue enriched expression gens
setwd("/home/zhaoying/BENAGEN/pydeseq2/")
data <- read.table("/home/zhaoying/BENAGEN/pydeseq2/spf_gene.csv",sep=",",header = T)
data$ID <- gsub("\\.\\d+","",data$ID)
head(data)

## tissue spf expression enrich
for (tissue in tissues) {
  ego <- enrichGO(data[data$tissue==tissue,]$ID,keyType = 'ENSEMBL',OrgDb = org.Hs.eg.db,ont = "BP")
  #  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  p1 <- dotplot(ego) + scale_color_material(palette = c("blue"),reverse = TRUE)
  ggsave(p1,filename = paste0("/home/zhaoying/BENAGEN/pydeseq2/GOanalysis/",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("/home/zhaoying/BENAGEN/pydeseq2/GOanalysis/",tissue,"_ego.csv"))
}

dotplot(ego) + scale_color_material(palette = c("blue"),reverse = TRUE)


### tissue spf m6A site hosted gene enrichment
setwd("/home/zhaoying/BENAGEN/methylKit/analysis/spf_GO")
data <- read.table("./tissue_spf.tsv",sep="\t",header = T)
for (tissue in tissues) {
  ego <- enrichGO(data[data$tissue==tissue,]$geneID,keyType = 'SYMBOL',OrgDb = org.Hs.eg.db,ont = "BP")
  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  ggsave(p1,filename = paste0("/home/zhaoying/BENAGEN/methylKit/analysis/spf_GO/",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("/home/zhaoying/BENAGEN/methylKit/analysis/spf_GO/",tissue,"_ego.csv"))
}
head(data)
#

### fetal spf m6A hosted gene enrichment
setwd("/home/zhaoying/BENAGEN/methylKit/analysis/NAR/")
data <- read.table("./fetal_spf_geneID.tsv",sep="\t",header=T)
head(data)
for (tissue in tissues) {
  ego <- enrichGO(data[data$tissue==tissue,]$geneID,keyType = 'SYMBOL',OrgDb = org.Hs.eg.db,ont = "BP")
  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  ggsave(p1,filename = paste0("/home/zhaoying/BENAGEN/methylKit/analysis/NAR/",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("/home/zhaoying/BENAGEN/methylKit/analysis/NAR/",tissue,"_ego.csv"))
}
#

ego_spf <- simplify(ego)
dotplot(ego_spf)

### tissue diff m6A host gene enrichment
setwd("/workd/workbai/NC/Figure4/m6A_methylKit/")
data <-read.table("./tissue_diff.csv",header = T,sep=",")
head(data)
data$index <- gsub("\\.\\d+","",data$index)
dim(data)
head(data)
ego <- enrichGO(unique(data$index) ,keyType = 'ENSEMBL',OrgDb = org.Hs.eg.db,ont = "BP")
p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
ego1 <- simplify(ego)
p1
ggsave(p1,filename =  "m6A_tissue_diff_ego.pdf", width = 6,height = 6)
write.csv(ego@result,file = "m6A_tissue_diff_ego.csv")
head(ego@result)

d1 <- read.table("/workd/workbai/NC/Figure4/m6A_methylKit/m6A_tissue_diff_ego2.csv",sep=",",header = T,row.names = 1)
# head(d1)
# ego@result <- d1
# dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
head(d1)

d1 <- head(d1,10)
d1 <- d1 %>%
  mutate(
    GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)) / 
      as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)),
    log10_qval = -log10(p.adjust)
  ) %>%
  arrange(desc(log10_qval)) %>%
  mutate(Description = factor(Description, levels = rev(Description)))
head(d1)
p2 <- ggplot(d1, aes(x = GeneRatio_num, y = Description, size = Count, color = log10_qval)) +
  geom_point(alpha = 1)  + scale_colour_material(palette = c("purple"),reverse = FALSE, name="-log10(P.adjust)") + theme_bw(15) + 
  labs(x="GeneRatio")  +scale_size_continuous(range = c(4, 10)) 
ggsave(p2,filename =  "m6A_tissue_diff_ego2.pdf", width = 12,height = 6)
  # scale_color_

#

### tissue diff m5C host gene enrichment
setwd("/workd/workbai/NC/Figure4/m5C_methylKit/")
data <-read.table("./tissue_diff.csv",header = T,sep=",")
head(data)
data$index <- gsub("\\.\\d+","",data$index)
dim(data)
ego <- enrichGO(unique(data$index) ,keyType = 'ENSEMBL',OrgDb = org.Hs.eg.db,ont = "BP")
p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)

p1
ggsave(p1,filename =  "m5C_tissue_diff_ego.pdf", width = 6,height = 6)
write.csv(ego@result,file = "m5C_tissue_diff_ego.csv")
#






### tissue specific m5C hosted gene enrichment analysis
setwd("/home/zhaoying/BENAGEN/CHEUI/tissue_spf/")
data <- read.table("./tissue_spf.tsv",sep="\t",header = T)
head(data)
data$ENSG <- gsub("\\.\\d+","",data$ENSG)
for (tissue in tissues) {
  ego <- enrichGO(data[data$tissue==tissue,]$ENSG ,keyType = 'ENSEMBL',OrgDb = org.Hs.eg.db,ont = "BP")
  p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
  ggsave(p1,filename = paste0("/home/zhaoying/BENAGEN/CHEUI/tissue_spf/",tissue,"_ego.pdf"),width = 6,height = 6)
  write.csv(ego@result,file = paste0("/home/zhaoying/BENAGEN/CHEUI/tissue_spf/",tissue,"_ego.csv"))
}


#### m5C analysis for tissue different
setwd("/home/zhaoying/BENAGEN/CHEUI/methylKit")
files <- list.files()

# 定义组织列表
tissues <- c("brain", "heart", "liver", "lung", "muscle", "kidney", "spleen")


combinations <- combn(tissues, 2)
combinations <- t(combinations)

results <- list()
for (i in 1:nrow(combinations)) {
  m <- combinations[i, 1]
  n <- combinations[i, 2]
  
  f1 <- files[str_detect(files, paste0( m, "$"))]
  
  f2 <- files[str_detect(files, paste0( n, "$"))]
  
  file_ls <- as.list(c(f1, f2))
  
  contrast <- c(rep(1, length(f1)), rep(0, length(f2)))
  myobj <- methRead(file_ls,sample.id=file_ls,treatment=contrast,mincov=10,assembly='hg38')
  meth=unite(myobj, destrand=FALSE,min.per.group=1L)
  myDiff=calculateDiffMeth(meth,mc.cores=2,adjust="hochberg")
  write.table(myDiff,file=paste0(m,"_",n,".tsv"),row.names=FALSE,quote=F,sep="\t")
  results[[i]] <- list(file_list = file_ls, contrast = contrast)
}



### tissue different m5C hosted gene enrichment analysis
setwd("/home/zhaoying/BENAGEN/CHEUI/methylKit")
data <- read.table("/home/zhaoying/BENAGEN/CHEUI/methylKit/tissue_diff.csv",sep=",",header = T)
head(data)
data$ENSG <- gsub("\\.\\d+","",data$ENSG)
dim(data)
ego <- enrichGO(data$ENSG ,keyType = 'ENSEMBL',OrgDb = org.Hs.eg.db,ont = "BP")
p1 <- dotplot(ego) + scale_colour_material(palette = c("purple"),reverse = TRUE)
# dim(data)
p1
ggsave(p1,filename =  "tissue_diff_ego.png", width = 6,height = 6)
write.csv(ego@result,file = "tissue_diff_ego.csv")





usethis::create_package("/home/zhaoying/worka/testR")




setwd("/home/zhaoying/BENAGEN/ISOQUANT/merged/suppa2/diffSplice")
data <- read.table("./AS_diffSplice_res.tsv",sep="\t",header = T)
head(data)
data$gene <- gsub("\\.\\d+","",data$gene)
data <- data[!duplicated(data$gene),]
ego <- enrichGO(data$gene ,keyType = 'ENSEMBL',OrgDb = org.Hs.eg.db,ont = "BP")
p1 <- dotplot(ego)  + scale_color_material(palette = c("blue"),reverse = TRUE)
p1
ggsave(p1,filename =  "AS_go.pdf", width = 6,height = 6)
write.csv(ego@result,file = "AS_go.csv")

barplot(ego)

