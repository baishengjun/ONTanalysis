library(devtools)
library(devtools)
devtools::install_github("karakulahg/regulaTER")

# devtools::install_github("karakulahg/TEffectR")
library(TEffectR)
library(stringr)
library(biomaRt)
library(biomartr)
library(dplyr)
library(Rsamtools)
library(edgeR)
library(rlist)
library(limma)
library(TEffectR)

repeatmasker.annotation <- TEffectR::rm_format(filepath = "~/reference/hg38/hg38.fa.out" )

exprs <- read.csv("/workd/workbai/NC/02.Quantification/DRS/DRS_estCount_samples.tsv",sep="\t", row.names = 1, header=T, stringsAsFactors = F)
head(exprs)
biomaRt::listEnsemblArchives()    
gene.annotation <- get_intervals(x = rownames(exprs), assembly="hg38", ID.type = "ensembl_gene_id", URL="http://dec2014.archive.ensembl.org" ) 
gene.annotation <- read.table("/workd/workbai/NC/08.TE/gene_annotation.tsv",sep="\t",header = T)
head(gene.annotation)

overlaps <- TEffectR::get_overlaps(g=gene.annotation, r=repeatmasker.annotation, strand = "strandness", distance = 5000, repeat_type = "LTR")
head(overlaps)
write.table(overlaps,sep="\t",file = "/workd/workbai/NC/08.TE/TE_geen.tsv",quote = F, row.names = F)
BAM.list <- c("/home/zhaoying/BENAGEN/BN220101NJ02S01N1-30RNAseq/22/sample6-brain.bam",
              "/home/zhaoying/BENAGEN/BN220101NJ02S01N1-30RNAseq/22/sample7-brain.bam")
TE.counts <- TEffectR::count_repeats(bamlist = BAM.list, 
                                     namelist = c("S6-brain","S7-brain"), ranges=overlaps)

TE_counts <- read.table("/home/zhaoying/workbai/NC/08.TE/TE_count.tsv",sep="\t",header = T,row.names = 1)
# head(TE_counts)

SumOfTEs <- TEffectR::summarize_repeat_counts(counts = TE_counts, namelist = colnames(TE_counts))
 
setwd("/workd/workbai/NC/08.TE/TEffectR-master/sampleInputs")
# read your gene annotation file
gene.annot<-read.table("gene.annotation.tsv", header= T, stringsAsFactors = F)

# read your gene expression file
gene.counts<-read.table("gene.counts.tsv", header= T, row.names=1, stringsAsFactors = F)

# read your summarised repeat annotation file
sum.repeat.counts<-read.table("sum.repeat.counts.tsv", header= T, stringsAsFactors = F)

# include covariates
covariates <- data.frame("TissueType" = c(rep("N",5), rep("T",6)) ) 

# OR

# in case the TE expression is only single predictor
covariates <- NULL

prefix<-"SampleRun"

# apply linear modeling function of TEffectR
lm<-TEffectR::apply_lm(gene.annotation = gene.annot,
                       gene.counts = gene.counts,
                       repeat.counts = sum.repeat.counts,
                       covariates = covariates,
                       prefix = prefix)





# devtools::install_github('robertamezquita/marge', ref = 'master')
library(regulaTER)
library(biomartr)
library(GenomicRanges)
library(dplyr)
library(biomaRt)
library(marge)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
raw.rmsk<-biomartr::read_rm("~/reference/hg38/hg38.fa.out")
# peak <-readPeakFile("/workd/workbai/NC/08.TE/roadmapData/E088_15_coreMarks_hg38lift_dense.active.bed")
peak <-readPeakFile("/workd/workbai/NC/08.TE/ATACseq/N8_peaks.narrowPeak")
peak
# names(mcols(peak))[8] <- "summit"
base_path <- "/workd/workbai/NC/08.TE/hg38_regions"
pathList <- list("Promoter" = paste0(base_path, "/hg38_promoter_complement.bed"),
                 "Exon" =  paste0(base_path, "/hg38_exons_complement.bed"),
                 "Intron" =  paste0(base_path, "/hg38_introns_complement.bed"),
                 "5UTR" =  paste0(base_path, "/hg38_5prime_complement.bed"),
                 "3UTR" =  paste0(base_path, "/hg38_3prime_complement.bed"),
                 "Downstream" =  paste0(base_path, "/hg38_downstream_complement.bed"),
                 "genomeSizePath" =  paste0(base_path, "/hg38.chrom.sizes"))
peakanno <- annotatePeak(peak,TxDb = txdb)
peakanno <- peakanno@anno
names(mcols(peakanno))[2] <- "summit"
head(peak)
enrich.peak <- TEAR(inputPeakFile = peakanno, pathList = pathList, numberOfShuffle = 2, repeatMaskerFile = raw.rmsk, 
                    format="narrow", minoverlap=0L, alternative = "greater", minobserved = 10)
write.table(enrich.peak$RepeatName, file = "/home/zhaoying/workbai/NC/08.TE/regulaTER/res/brain_enrich_peak_spf.tsv",sep="\t",quote = F,row.names = F)
genes <- read.table("/workd/workbai/NC/08.TE/regulaTER/brain_spf.tsv",sep="\t")
# head(genes)
head(genes)
genes <- genes[,c(2:7)]
colnames(genes) <- c("geneID","geneName","seqnames","start","end","strand")
genes$geneID <- sub("\\.\\d+$", "", genes$geneID)
head(genes)
IdDEGRepeats <- DATE(enrichTEARResult = enrich.peak, peaks = peakanno, rmsk = raw.rmsk, genes = genes, 
                     alternative = "greater" ,numberOfShuffle = 2, minobserved = 10, distance = 100000)
write.table(IdDEGRepeats, file = paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/res/","brain","_IdDEGRepeats.tsv"),
            sep="\t",quote = F,row.names = F)

options('homer_path' = '/home/zhaoying/')
check_homer()
FindMotifs(df = IdDEGRepeats, repeatMaskerFile = raw.rmsk, peak = peakanno, distance = 100000, genes = genes, genome = "hg38", 
           outDir = "/workd/workbai/NC/08.TE/regulaTER/res/res/brain_spf", type = "linkedRepeats",
           topRepeats = T,homerPath = '/home/zhaoying/' )
# FindMotifs(df = enrich.peak$RepeatName, repeatMaskerFile = raw.rmsk, peak = peak, genome = "hg38",
#            outDir = "../test/", homerPath = '/home/zhaoying/', type = "enrichPeak", topRepeats = T)


setwd("/home/zhaoying/workbai/NC/08.TE/regulaTER/res")
tissues <- c("liver","muscle","kidney","heart","spleen","lung")
for (i in tissues) {
  peak <-readPeakFile(paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/",i,"_atac.bed"))
  peak <- annotatePeak(peak,TxDb = txdb)
  peak <- peak@anno
  names(mcols(peak))[2] <- "summit"
  enrich.peak <- TEAR(inputPeakFile = peak, pathList = pathList, numberOfShuffle = 2, repeatMaskerFile = raw.rmsk, 
                      format="narrow", minoverlap=0L, alternative = "greater", minobserved = 10)
  write.table(enrich.peak$RepeatName, file = paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/res/",i,"_enrich_peak.tsv"),
              sep="\t",quote = F,row.names = F)
  genes <- read.table(paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/",i,"_exp.tsv"),sep="\t")
  genes <- genes[,c(2:7)]
  colnames(genes) <- c("geneID","geneName","seqnames","start","end","strand")
  genes$geneID <- sub("\\.\\d+$", "", genes$geneID)
  IdDEGRepeats <- DATE(enrichTEARResult = enrich.peak, peaks = peak, rmsk = raw.rmsk, genes = genes, 
                       alternative = "greater" ,numberOfShuffle = 2, minobserved = 10, distance = 100000)
  write.table(IdDEGRepeats, file = paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/res/",i,"_IdDEGRepeats.tsv"),
              sep="\t",quote = F,row.names = F)
  FindMotifs(df = IdDEGRepeats, repeatMaskerFile = raw.rmsk, peak = peak, distance = 100000, genes = genes, genome = "hg38", 
             outDir = i, type = "linkedRepeats", topRepeats = T,homerPath = '/home/zhaoying/' )
  
}


setwd("/home/zhaoying/workbai/NC/08.TE/regulaTER/res")
tissues <- c("liver","muscle","kidney","heart","spleen","lung")
tissues <- c("muscle","kidney","heart","spleen","lung")

for (i in tissues) {
  peak <-readPeakFile(paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/ATAC_",i,".bed"))
  peakanno <- annotatePeak(peak,TxDb = txdb)
  peakanno <- peakanno@anno
  names(mcols(peakanno))[2] <- "summit"
  enrich.peak <- TEAR(inputPeakFile = peakanno, pathList = pathList, numberOfShuffle = 2, repeatMaskerFile = raw.rmsk, 
                      format="narrow", minoverlap=0L, alternative = "greater", minobserved = 10)
  write.table(enrich.peak$RepeatName, file = paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/res/",i,"_enrich_peak_spf.tsv"),
              sep="\t",quote = F,row.names = F)
  genes <- read.table(paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/",i,"_spf.tsv"),sep="\t")
  genes <- genes[,c(2:7)]
  colnames(genes) <- c("geneID","geneName","seqnames","start","end","strand")
  genes$geneID <- sub("\\.\\d+$", "", genes$geneID)
  IdDEGRepeats <- DATE(enrichTEARResult = enrich.peak, peaks = peakanno, rmsk = raw.rmsk, genes = genes,
                       alternative = "greater" ,numberOfShuffle = 2, minobserved = 10, distance = 100000)
  write.table(IdDEGRepeats, file = paste0("/home/zhaoying/workbai/NC/08.TE/regulaTER/res/",i,"_IdDEGRepeats_spf.tsv"),
              sep="\t",quote = F,row.names = F)
  FindMotifs(df = IdDEGRepeats, repeatMaskerFile = raw.rmsk, peak = peakanno, distance = 100000, genes = genes, genome = "hg38",
             outDir = paste0("res/",i,"_spf"), type = "linkedRepeats", topRepeats = T,homerPath = '/home/zhaoying/' )
  # FindMotifs(df = enrich.peak$RepeatName, repeatMaskerFile = raw.rmsk, peak = peakanno, genome = "hg38",
  #            outDir = i, homerPath = '/home/zhaoying/', type = "enrichPeak", topRepeats = T)
}













