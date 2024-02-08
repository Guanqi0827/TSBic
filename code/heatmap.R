
##看细胞在基因的表达情况(小规模数据集)
library(pheatmap)
library(RColorBrewer)
colorName <- c(brewer.pal(9,"OrRd")[6], brewer.pal(9,"Set1")[2])
colorName <- c('#d86967', '#58539f')
colorName <- c(rgb(243,104,176, maxColorValue = 255), rgb(92,219,253, maxColorValue = 255))
colorName <- c(rgb(146,197,222, maxColorValue = 255), rgb(217,83,79, maxColorValue = 255))
barplot(rep(1,2), 1, col = colorName)
display.brewer.pal(9,"OrRd")
Data_3_simulation <- Data_3[Simulation_GENE_index,]
pheatmap(Data_3_simulation, cluster_rows = F, cluster_cols = F, color = colorName, legend = F, cex = 0.9)
#画细胞在基因表达情况热图
drawBiclusterExpressed <- function(bicluster = B){
  library(pheatmap)
  colorName <- c(rgb(146,197,222, maxColorValue = 255), rgb(217,83,79, maxColorValue = 255))
  gene_index <- bicluster[[1]]
  Cell_names <- bicluster$`cell names`
  bicluster_expressed <- Data_3_simulation[sort(gene_index),which(colnames(Data_3) %in% Cell_names)]
  if(sum(bicluster_expressed) < length(gene_index) * length(Cell_names)){
    pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName, legend = F)
  }
  else{
    pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1)
  }
}

B <- GAbiclustering_fate_result[[1]]
B <- GAbiclustering_all_result[[1]]
B <- All_B[[1]]
drawBiclusterExpressed(B)
B <- All_B[[2]]
drawBiclusterExpressed(B)
B <- All_B[[3]]
drawBiclusterExpressed(B)
B <- All_B[[4]]
drawBiclusterExpressed(B)
B <- All_B[[5]]
drawBiclusterExpressed(B)
B <- All_B[[6]]
drawBiclusterExpressed(B)
B <- All_B[[7]]
drawBiclusterExpressed(B)
B <- All_B[[8]]
drawBiclusterExpressed(B)
B <- All_B[[9]]
drawBiclusterExpressed(B)
B <- All_B[[10]]
drawBiclusterExpressed(B)









                                                 