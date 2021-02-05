# Jesus Garcia Garcia
# Benchmarking Single-cell RNA-seq Copy Number Variation (CNV) Computational Biology Tools
# Updated - 09/29/2020
# garci624@umn.edu
# Run "simulateData" first

library(bench)
library(ggplot2)
library(Seurat)

location <- getwd()
setwd(location)

gene_pos_filename <- "gencode_v19_gene_pos.txt"
count_matrix <- "oligodendroglioma_expression_downsampled.counts.matrix"
anno_filename <- "oligodendroglioma_annotations_downsampled.txt"
selected_chr = c('chr1', 'chr2', 'chr3',
                 'chr4', 'chr5', 'chr6',
                 'chr7', 'chr8', 'chr9',
                 'chr10', 'chr11', 'chr12',
                 'chr13', 'chr14', 'chr15',
                 'chr16', 'chr17', 'chr18',
                 'chr19', 'chr20', 'chr21',
                 'chr22')

mbps = 9000000
selected_tumor_init = c("MGH36", "MGH54", "MGH53")
#add_or_red = "RED" # AMP | RED
val_add = 100
val_red = (1/0)

# import gene positions
gene_pos <- read.delim(gene_pos_filename, header = FALSE)

mat <- read.delim(count_matrix)

# import original annotations
ann <- read.delim(anno_filename, header=FALSE)
# original annotation
orig_ann <- ann[ which(ann$V2=='Microglia/Macrophage' | ann$V2=='Oligodendrocytes (non-malignant)'),1]

# original annotation's'
orig_anns <- ann[ which(ann$V2=='Microglia/Macrophage' | ann$V2=='Oligodendrocytes (non-malignant)'),]

# original matrix with only 'Microglia/Macrophage' and 'Oligodendrocytes (non-malignant).'
orig_mat <- mat[ ,orig_ann]

# start of simulated data
sim_add_mat <- orig_mat # copy of simulated data
sim_del_mat <- orig_mat
coln = colnames(sim_add_mat)

for (curr_chr in selected_chr){
  chr <- subset(gene_pos, V2 %in% curr_chr & V1 %in% row.names(sim_add_mat))
  total_chr = nrow(chr)
  
  chr_start = chr$V3[1]
  chr_end = chr$V4[total_chr]
  
  if (chr_end - chr_start <= mbps){
    reg_start_idx = 1
    reg_end_idx = total_chr
    
  }else{
    reg_start_idx = as.integer(runif(1, min = 0, max = total_chr/4))
    reg_start = chr$V3[reg_start_idx]
    reg_end_lim = reg_start + mbps
    
    match_list = chr$V3[chr$V3 >= reg_end_lim]
    reg_end_idx = which(chr$V3 == match_list[1])
  }
  
  chr_genes <- chr$V1[reg_start_idx: reg_end_idx]
  
  combined = c()
  for (val in selected_tumor_init){
    MGH_vals = coln[grepl( val , coln )]
    combined = c(combined, MGH_vals)
  }
  
  #print(curr_chr)
  for (gene in chr_genes){
    print(gene)
    #if ((gene %in% row.names(sim_mat))) {
    sim_add_mat[gene,combined]=(sim_add_mat[gene,combined]+1)*val_add
    sim_del_mat[gene,combined]=sim_del_mat[gene,combined]/val_red
    #}
  }
}

# final matrix

# generating new annotation file
# insert 'Tumor' in simulated data
colnames(sim_add_mat) <- paste("Tumor", colnames(sim_add_mat), sep = "_")
colnames(sim_del_mat) <- paste("Tumor", colnames(sim_del_mat), sep = "_")

# original matrix with only 'Microglia/Macrophage' and 'Oligodendrocytes (non-malignant).' + simulated tumor data
final_add_mat <- cbind(orig_mat,sim_add_mat)
write.csv(final_add_mat,"final_add_mat.csv")

final_del_mat <- cbind(orig_mat,sim_del_mat)
write.csv(final_del_mat,"final_del_mat.csv")

# 'Tumor' + 'Normal' annotations
tumor_ann <- data.frame("V1" = colnames(sim_add_mat), "V2" = "Tumor")
# tumor_ann <- data.frame("V1" = colnames(sim_mat), "V2" = "Tumor")
final_ann <- rbind(orig_anns, tumor_ann) # final annotations

# 'final TXT' annotation file
write.table(final_ann, "final_annotation.txt", sep = "\t",row.names = FALSE, col.names = FALSE)

