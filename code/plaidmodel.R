#计算基因表达速率的矩阵53*145
#CELL_fate_blot_velocity_174是只取正增量，这里是取全部增量均值

CELL_fate_blot_velocity_53 <- matrix(nrow = length(GENE_names),ncol = length(CELL))
for(i in 1:length(GENE_names)){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  for(j in 1:length(CELL)){
    index <- which(realdata[[ind]][[1]] %in% CELL_fate_names[j])
    if(length(index) == 0){
    }
    else if(length(index) == 1){
      CELL_fate_blot_velocity_53[i,j] <- realdata[[ind]][[3]][index]
    }
    else{
      d <- diff(realdata[[ind]][[3]][index])
      d <- d
      CELL_fate_blot_velocity_53[i,j] <- mean(d)
    }
  }
}

rownames(CELL_fate_blot_velocity_53) <- GENE_names
colnames(CELL_fate_blot_velocity_53) <- CELL
CELL_fate_blot_velocity_53 <- as.data.frame(CELL_fate_blot_velocity_53)
CELL_fate_blot_velocity_53_scaled <- scale(t(CELL_fate_blot_velocity_53),center = T,scale = T)
CELL_fate_blot_velocity_53_scaled <- t(CELL_fate_blot_velocity_53_scaled)
CELL_fate_blot_velocity_53_scaled <- as.data.frame(CELL_fate_blot_velocity_53_scaled)


library(DMwR2)
library(DescTools)
PlotMiss(t(CELL_fate_blot_velocity_53))
length(which(is.na(CELL_fate_blot_velocity_53)))


##双聚类
library(biclust)
set.seed(2021)
##不做标准化
##用CC双聚类，缺失值补全方法是用矩阵取值范围内的均匀分布的随机数
ma <- max(CELL_fate_blot_velocity_53,na.rm = T)
mi <- min(CELL_fate_blot_velocity_53,na.rm = T)
CELL_fate_blot_velocity_53_nona <- CELL_fate_blot_velocity_53

for(i in 1:length(GENE_names)){
  for(j in 1:length(CELL)){
    if(is.na(CELL_fate_blot_velocity_53_nona[i,j])){
      CELL_fate_blot_velocity_53_nona[i,j] <- runif(n = 1, min = mi, max = ma)
    }
  }
}
length(which(is.na(CELL_fate_blot_velocity_53_nona)))


biclust_53_velocity <- biclust(as.matrix(CELL_fate_blot_velocity_53_nona),method = BCCC())
summary(biclust_53_velocity)
bubbleplot(as.matrix(CELL_fate_blot_velocity_53_nona),biclust_53_velocity,showLabels = T)
parallelCoordinates(as.matrix(CELL_fate_blot_velocity_53_nona),bicResult = biclust_53_velocity,number = 1,plotBoth = T)
drawHeatmap(as.matrix(CELL_fate_blot_velocity_53_nona),biclust_53_velocity,5)







##用Plaid model双聚类，缺失值补全方法是行均值+列均值-全体均值

cme <- c()
rme <- c()
me <- mean(a, na.rm = T)
CELL_fate_blot_velocity_53_nona <- CELL_fate_blot_velocity_53
a <- as.matrix(CELL_fate_blot_velocity_53)
length(which(is.na(CELL_fate_blot_velocity_53)))


for(i in 1:length(GENE_names)){
  rme <- c(rme, mean(a[i,], na.rm = T))
}

for(j in 1:length(CELL)){
  cme <- c(cme, mean(a[,j], na.rm = T))
}


for(i in 1:length(GENE_names)){
  for(j in 1:length(CELL)){
    if(is.na(a[i,j])){
      x <- rme[i] + cme[j] - me
      CELL_fate_blot_velocity_53_nona[i, j] <- x
    }
  }
}


length(which(is.na(CELL_fate_blot_velocity_53_nona)))
CELL_fate_blot_velocity_53_nona_scaled <- scale(t(CELL_fate_blot_velocity_53_nona),center = T,scale = T)
CELL_fate_blot_velocity_53_nona_scaled <- t(CELL_fate_blot_velocity_53_nona_scaled)
CELL_fate_blot_velocity_53_nona_scaled <- as.data.frame(CELL_fate_blot_velocity_53_nona_scaled)


set.seed(2)
biclust_53_velocity <- biclust(as.matrix(CELL_fate_blot_velocity_53_nona_scaled),method = BCPlaid(), row.release = 0.6, col.release = 0.6, iter.layer = 20)
summary(biclust_53_velocity)
bubbleplot(as.matrix(CELL_fate_blot_velocity_53_nona_scaled),biclust_53_velocity,showLabels = T)
parallelCoordinates(as.matrix(CELL_fate_blot_velocity_53_nona_scaled),bicResult = biclust_53_velocity,number = 1,plotBoth = T)
drawHeatmap(as.matrix (CELL_fate_blot_velocity_53_nona_scaled),biclust_53_velocity,1)



setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")
biclust_53_velocity_plaid <- biclust_53_velocity
biclust_53_velocity <- biclust_53_velocity_plaid
library(pheatmap)
colorName <- c("steelblue", 'tomato')
num <- dim(biclust_53_velocity@RowxNumber)[2]

for(i in 1:num){
  gene_index <- which(biclust_53_velocity@RowxNumber[,3] == TRUE)
  cell_index <- which(biclust_53_velocity@NumberxCol[3,] == TRUE)
  output_names <- paste(i,'.pdf',sep = '')
  bicluster_expressed <- Data_3_simulation[sort(gene_index), cell_index]
  if(sum(bicluster_expressed) != length(gene_index) * length(cell_index)){
    pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName, legend = F)
  }
  else{
    pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1)
  }
}
B <- All_B[[1]]
cell_index <- B$`cell index`
cell_index <- which(biclust_53_velocity@NumberxCol[3,] == TRUE)
CatCellFate(cell_index)

a <- CatCellFate(cell_index)[[1]]
b <- CatCellFate(cell_index)[[2]]
paste0(a,'(',b,')')



