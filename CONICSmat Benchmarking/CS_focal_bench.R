
# Jesus Garcia Garcia
# Benchmarking Single-cell RNA-seq Copy Number Variation (CNV) Computational Biology Tools
# Updated - 09/29/2020
# garci624@umn.edu
# Run "simulateData" first

library(CONICSmat)
library(Seurat)
library(bench)

##################################################################
final_add_mat <- read.csv("/Users/garci624/Desktop/CNV/final_add_mat.csv", row.names=1,sep=",")
#final_del_mat <- read.csv("/Users/garci624/Desktop/CNV/final_del_mat.csv", row.names=1,sep=",")

# return(final_mat)

# CONICSmat function being benchmarked
CS = function(dat){
  a=NormalizeData(dat)
  
  meta=read.table("/Users/garci624/Desktop/CNV/final_annotation.txt",header=F,row.names=1,check.names = F,sep="\t")
  meta=meta[colnames(a),]
  unique(meta)
  
  #hi = function(pat){
  patient="MGH"
  mat=a[,which(meta=="Tumor")]
  n=a[,which(meta=="Microglia/Macrophage"  | meta=="Oligodendrocytes (non-malignant)")]
  mat=cbind(n,mat)
  
  c <- colnames(mat)
  
  normal=1:ncol(n)
  tumor=(ncol(n)+1):ncol(mat)
  
  regions=read.table("/Users/garci624/Desktop/CNV/CONICSmat/chromosome_full_positions_grch38.txt",sep="\t",row.names = 1,header = T)
  gene_pos=getGenePositions(rownames(mat))
  mat=filterMatrix(mat,gene_pos[,"hgnc_symbol"],minCells=10)
  normFactor=calcNormFactors(mat)
  
  p <- plotChromosomeHeatmap(mat,normal = normal, plotcells = c[grepl( '^Tumor' , c )] ,
                             gene_pos = gene_pos,chr=T,windowsize = 121, expThresh=0.2, thresh = 1)
}

##################################################################
CONICSmat <- bench::mark(iterations = 100, CS(final_add_mat)) # Run for 100 iterations
summary(CONICSmat) # Outputs benchmarking stats


