setwd("/home/zhaoying/workbai/NC/Figure3/metaplotR")
library(clusterProfiler,lib.loc = "/home/zhaoying/R_lib/")
library(ggsci)
library(ggthemes,lib.loc = "/home/zhaoying/R_lib/")
library(enrichplot,lib.loc = "/home/zhaoying/R_lib/")
library(org.Hs.eg.db)
library(scales)
library("ChIPseeker",lib.loc = "/home/zhaoying/R_lib/")
library(scales,lib.loc = "/home/zhaoying/R_lib/")
tissues <- c("brain","muscle","lung","heart","spleen","kidney","liver","concat","common","spf")
ls_res = list()
# tissues <- paste0("/home/zhaoying/BENAGEN/metaPlotR/annot_",tissues,".dist.measures.txt")
for (file in tissues) {
  for (mod in c("m6A","m5C","psU")) {
    
    filetxt <-  paste0("/home/zhaoying/workbai/NC/Figure3/metaplotR/",mod,"_",file,"_dist.measures.txt")
    m6a.dist <-  read.delim(filetxt,header=T)
    head(m6a.dist)
    # Determine longest length transcript for each gene
    trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
    temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
    colnames(temp) <- c("gene_name", "gid", "trx_len")
    temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
    temp.df <- temp[!duplicated(temp$gene_name),]
    
    # limit m6a data to one transcript per gene (longest)
    m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]
    
    # View size of our dataset (rows, columns)
    dim(m6a.dist)
    qplot(m6a.dist$rel_location, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()
    summary(data.frame(m6a.dist$utr5_size, m6a.dist$cds_size, m6a.dist$utr3_size))
    utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
    utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
    # assign the regions to new dataframes
    utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
    cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
    utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
    
    # rescale 5'UTR and 3'UTR
    # library("scales")
    utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
    utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
    # Combine and plot
    ## Histogram
    m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
    assign(paste0(mod,"_",file,".coord"),m6a.metagene.coord)
    qplot(m6a.metagene.coord, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()
    qplot(m6a.metagene.coord, geom="freqpoly") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
    p <- qplot(m6a.metagene.coord, geom="density") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
    write.table(m6a.metagene.coord,paste0("/home/zhaoying/workbai/NC/Figure3/metaplotR/result/",mod,"_",file,".txt"),row.names = FALSE)
    # ls_res <- append(ls_res,m6a.metagene.coord)
    print(p)
  }
}
# 
# metagene.cord <- c(brain.coord,liver.coord,heart.coord,lung.coord,spleen.coord,muscle.coord,kidney.coord)
# mod <- c(rep("brain",length(brain.coord)),
#          rep("liver",length(liver.coord)),
#          rep("heart",length(heart.coord)),
#          rep("lung",length(lung.coord)),
#          rep("spleen",length(spleen.coord)),
#          rep("muscle",length(muscle.coord)),
#          rep("kidney",length(kidney.coord)))
# df <- data.frame(metagene.cord, mod)
# df
# p1 <- ggplot(df) + geom_density(aes(x = metagene.cord, colour = mod),) + xlim(0.6, 3) +
#   theme_bw() + geom_vline(xintercept = 1:2, col = "grey",linetype="dashed") + scale_color_d3()
# p1
# ggsave(p1,filename = "/home/zhaoying/workbai/NC/Figure3/metaplotR/tissue_metagenePlot.pdf",height = 6,width = 8)
# 
# 
# 
# 
# tissues <- c("brain","muscle","lung","heart","spleen","kidney","liver","tissue_spf","tissue_common")
# ###tissue-specific metagene plot
# for (file in tissues) {
#   filetxt <-  paste0("/home/zhaoying/BENAGEN/metaPlotR/tissue_spf/annot_",file,".dist.measures.txt")
#   m6a.dist <-  read.delim(filetxt,header=T)
#   head(m6a.dist)
#   # Determine longest length transcript for each gene
#   trx_len <- m6a.dist$utr5_size + m6a.dist$cds_size + m6a.dist$utr3_size # Determine transcript length
#   temp <- data.frame(m6a.dist$gene_name, m6a.dist$refseqID, trx_len)
#   colnames(temp) <- c("gene_name", "gid", "trx_len")
#   temp.df <- temp[order(temp$gene_name, temp$gid, -temp$trx_len),]
#   temp.df <- temp[!duplicated(temp$gene_name),]
#   
#   # limit m6a data to one transcript per gene (longest)
#   m6a.dist <- m6a.dist[m6a.dist$refseqID %in% temp.df$gid,]
#   
#   # View size of our dataset (rows, columns)
#   dim(m6a.dist)
#   qplot(m6a.dist$rel_location, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()
#   summary(data.frame(m6a.dist$utr5_size, m6a.dist$cds_size, m6a.dist$utr3_size))
#   utr5.SF <- median(m6a.dist$utr5_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
#   utr3.SF <- median(m6a.dist$utr3_size, na.rm = T)/median(m6a.dist$cds_size, na.rm = T)
#   # assign the regions to new dataframes
#   utr5.m6a.dist <- m6a.dist[m6a.dist$rel_location < 1, ]
#   cds.m6a.dist <- m6a.dist [m6a.dist$rel_location < 2 & m6a.dist$rel_location >= 1, ]
#   utr3.m6a.dist <- m6a.dist[m6a.dist$rel_location >= 2, ]
#   
#   # rescale 5'UTR and 3'UTR
#   # library("scales")
#   utr5.m6a.dist$rel_location <- rescale(utr5.m6a.dist$rel_location, to = c(1-utr5.SF, 1), from = c(0,1))
#   utr3.m6a.dist$rel_location <- rescale(utr3.m6a.dist$rel_location, to = c(2, 2+utr3.SF), from = c(2,3))
#   # Combine and plot
#   ## Histogram
#   m6a.metagene.coord <- c(utr5.m6a.dist$rel_location, cds.m6a.dist$rel_location, utr3.m6a.dist$rel_location)
#   assign(paste0(file,".coord"),m6a.metagene.coord)
#   qplot(m6a.metagene.coord, geom="histogram") + geom_vline(xintercept = 1:2, col = "grey") + theme_bw()
#   qplot(m6a.metagene.coord, geom="freqpoly") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
#   p <- qplot(m6a.metagene.coord, geom="density") + geom_vline(xintercept = 1:2, col = "red") + theme_bw()
#   # ls_res <- append(ls_res,m6a.metagene.coord)
#   print(p)
# }
# metagene.cord <- c(brain.coord,liver.coord,heart.coord,lung.coord,spleen.coord,muscle.coord,kidney.coord)
# mod <- c(rep("brain",length(brain.coord)),
#          rep("liver",length(liver.coord)),
#          rep("heart",length(heart.coord)),
#          rep("lung",length(lung.coord)),
#          rep("spleen",length(spleen.coord)),
#          rep("muscle",length(muscle.coord)),
#          rep("kidney",length(kidney.coord)))
# df <- data.frame(metagene.cord, mod)
# df
# p1 <- ggplot(df) + geom_density(aes(x = metagene.cord, colour = mod),) + xlim(0.6, 3) +
#   theme_bw() + geom_vline(xintercept = 1:2, col = "grey",linetype="dashed") + scale_color_d3()
# p1
# ggsave(p1,filename = "/home/zhaoying/BENAGEN/metaPlotR/tissue_spf/tissue_spf_metagenePlot.pdf",height = 6,width = 8)
# 
# common_spf <- c(tissue_spf.coord,tissue_common.coord)
# mod <- c(rep("specific",length(tissue_spf.coord)),
#          rep("common",length(tissue_common.coord)))
# df2 <- data.frame(common_spf, mod)
# df2
# p2 <- ggplot(df2) + geom_density(aes(x = common_spf, colour = mod),) + xlim(0.6, 3) +
#   theme_bw() + geom_vline(xintercept = 1:2, col = "grey",linetype="dashed") + scale_color_d3()
# p2
# ggsave(p2,filename = "/home/zhaoying/BENAGEN/metaPlotR/tissue_spf/tissue_spf_common_metagenePlot.pdf",height = 6,width = 8)
