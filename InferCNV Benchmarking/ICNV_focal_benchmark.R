
# Jesus Garcia Garcia
# Benchmarking Single-cell RNA-seq Copy Number Variation (CNV) Computational Biology Tools
# Updated - 09/29/2020
# garci624@umn.edu
# Run "simulateData" first

##################################################################
final_add_mat <- read.csv("/Users/garci624/Desktop/CNV/final_add_mat.csv", row.names=1,sep=",")
#final_del_mat <- read.csv("final_del_mat.csv", row.names=1,sep=",")

library(bench)
library(ggplot2)
library(Seurat)
require(biomaRt) ## for gene coordinates
library(infercnv)

# InferCNV function being benchmark
##################################################################
gma <- final_add_mat

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=gma,
                                    annotations_file="/Users/garci624/Desktop/CNV/final_annotation.txt",
                                    delim="\t",
                                    gene_order_file="/Users/garci624/Desktop/CNV/gencode_v19_gene_pos.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))


IC = function(dat){
  Infer <- infercnv::run(dat,
                         cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                         
                         ##################################################################
                         out_dir="/Users/garci624/Desktop/CNV/InferCNV/ICNV_MB_AMP_8",  # dir is auto-created for storing outputs
                         cluster_by_groups=T,   # cluster
                         denoise=T,
                         HMM=F
  )
}

InferCNV <- bench::mark(iterations = 100, IC(infercnv_obj)) # Runs InferCNV for 100 iterations
summary(InferCNV) # Benhcmark summary



