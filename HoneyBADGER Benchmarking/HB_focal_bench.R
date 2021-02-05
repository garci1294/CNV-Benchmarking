
# Jesus Garcia Garcia
# Benchmarking Single-cell RNA-seq Copy Number Variation (CNV) Computational Biology Tools
# Updated - 09/29/2020
# garci624@umn.edu
# Run "simulateData" first

library(bench)
library(ggplot2)
library(Seurat)
require(biomaRt) ## for gene coordinates
library(HoneyBADGER)

##################################################################
final_add_mat <- read.csv("/Users/garci624/Desktop/CNV/final_add_mat.csv", row.names=1,sep=",")
#final_del_mat <- read.csv("/Users/garci624/Desktop/CNV/final_del_mat.csv", row.names=1,sep=",")

mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl') ## current version

HB = function(dat){
  
  gma <-NormalizeData(dat)
  # Extract MGH36, MGH53 and MGH54 tumors
  coln = colnames(gma)
  
  Tumor_MGH_vals = coln[grepl( 'Tumor_' , coln )]
  
  MGH_vals = coln[grepl( '^MGH' , coln )]
  
  # combined_ALL = c(Tumor_MGH_vals,MGH_vals)
  
  # Storing all tumors data in three_mgh
  # We wil be using the three_mgh data to create 'infercnv' object
  gexp <- gma[, Tumor_MGH_vals]
  ref <- gma[, MGH_vals]
  
  #print(mart.obj)
  
  hb <- new('HoneyBADGER', name='MGH')
  hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=TRUE, verbose=TRUE)
  #minMeanBoth=0, 
  #minMeanTest=mean(gexp[gexp!=0]),
  #minMeanRef=mean(gexp[gexpt!=0]))
  
  ##################################################################
  png("/Users/garci624/Desktop/CNV/HoneyBADGER/HB_MB_AMP_8.png",width = 1000,height = 400)
  hb$plotGexpProfile() ## initial visualization
  dev.off()
}

##################################################################
HoneyBADGER <- bench::mark(iterations = 100, HB(final_add_mat))
summary(HoneyBADGER)


