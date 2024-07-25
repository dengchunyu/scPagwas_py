#change the R data to files
#load panthway data
load("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas/data/Genes_by_pathway_kegg.RData")
gene_data <- data.frame(
  Gene_Name = names(Genes_by_pathway_kegg),
  Genes = unlist(sapply(Genes_by_pathway_kegg, paste, collapse = ","))
)

# Write the data frame to a CSV file
write.csv(gene_data, "/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/Genes_by_pathway_kegg.csv", row.names = FALSE)

load("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas/data/block_annotation.RData")
write.csv(block_annotation, "/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/block_annotation.csv", row.names = FALSE)


for(i in 1:22){
   a<-readRDS(paste0("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/LD/",i,".Rds"))
   file_name <- paste0("/Users/chunyudeng/Library/CloudStorage/OneDrive-共享的库-Onedrive/RPakage/scPagwas_py/data/LD/",i,".csv")
   write.csv(a,file_name,row.names = FALSE)
 }