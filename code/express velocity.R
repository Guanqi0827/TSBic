
#pha-4基因的第一个文件，在E中表达的细胞里的表达速度

#挑选出数据，在一个细胞表达后，就默认一直表达

setwd("/Users/yanxianzhong/Desktop/基因聚类/hujiedata" )
which(file_newnames%in%file_gene[1,1])
which(file_names_csv%in%file_gene[1,1])
pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[1,1])]]

pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[1,1])]]
pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])

for(i in 1:dim(pha4_1_onset)[1]){
  pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
  
}


#先画有表达的E细胞的散点图
library(RColorBrewer)
display.brewer.pal(n = 8, name = 'RdBu')
colorName <- brewer.pal(n = 8, name = 'RdBu')  

pha4_1_onset_maxblot <- 0
pha4_1_onset_minblot <- 0

for(i in 1:8){
  pha4_1_onset_maxblot[i] <- max(pha4_1_onset_realdata[[i]][[3]])
  pha4_1_onset_minblot[i] <- min(pha4_1_onset_realdata[[i]][[3]])
  
}
pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
pha4_1_onset_minblot <- min(pha4_1_onset_minblot)

plot(1,1,xlim = c(130,200),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'blot')
for(i in 1:8){
  points(pha4_1_onset_realdata[[i]]$time,pha4_1_onset_realdata[[i]]$blot,pch = 4,col = colorName[i])
  
}

legend('topleft',pch = 4,col = colorName, ncol=4, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[1:8], cex = 0.8, bty = 'n')

#计算增量并画增量图

diff(pha4_1_onset_realdata[[8]]$blot)
pha4_1_onset_realdata[[1]]$time[-(1:2)]

plot(1,1,xlim = c(130,200),ylim = c(-20000,20000),type = 'n',xlab = 'time',ylab = 'blot')
for(i in 1:8){
  points(pha4_1_onset_realdata[[i]]$time[-1],abs(diff(pha4_1_onset_realdata[[i]]$blot)),pch = 4,col = colorName[i])
}

legend('topleft',pch = 4,col = colorName, ncol=4, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[1:8], cex = 0.8, bty = 'n')




#
#
#
#计算基因表达速率的矩阵
#

cell_blot_velocity_184 <- matrix(nrow = 184,ncol = cell_number)
for(i in 1:184){
  ind <- which(file_names_csv%in%file_gene[i,1])
  for(j in 1:1083){
    index <- which(realdata[[ind]][[1]] %in% cell[j])
    if(length(index) == 0){
      cell_blot_velocity_184[i,j] <- 0
    }
    else if(length(index) == 1 | length(index) == 2){
      cell_blot_velocity_184[i,j] <- mean(realdata[[ind]][[3]][index]) 
    }
    else{
      cell_blot_velocity_184[i,j] <- mean(diff(diff(realdata[[ind]][[3]][index])))
    }
  }
}


rownames(cell_blot_velocity_184) <- rownames(cell_blot_184)
colnames(cell_blot_velocity_184) <- cell

##标准化
cell_blot_velocity_184_scaled <- scale(t(cell_blot_velocity_184),center = T,scale = T)
cell_blot_velocity_184_scaled <- t(cell_blot_velocity_184_scaled)



#计算基因表达速率的矩阵174*145
#

CELL_fate_blot_velocity_174 <- matrix(nrow = 174,ncol = 145)
for(i in 1:174){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  for(j in 1:145){
    index <- which(realdata[[ind]][[1]] %in% CELL_fate_names[j])
    if(length(index) == 0){
    }
    else if(length(index) == 1 & realdata[[ind]][[3]][index] >= 0){
      CELL_fate_blot_velocity_174[i,j] <- realdata[[ind]][[3]][index]
    }
    else if(length(index) == 1 & realdata[[ind]][[3]][index] < 0){
      CELL_fate_blot_velocity_174[i,j] <- 0
    } 
    else{
      d <- diff(realdata[[ind]][[3]][index])
      d <- d[d>=0]
      CELL_fate_blot_velocity_174[i,j] <- mean(d)
    }
  }
}


rownames(CELL_fate_blot_velocity_174) <- GENE_names
colnames(CELL_fate_blot_velocity_174) <- CELL_fate_names
CELL_fate_blot_velocity_174 <- as.data.frame(CELL_fate_blot_velocity_174)
PlotMiss(t(CELL_fate_blot_velocity_174))
length(which(is.na(CELL_fate_blot_velocity_174)))



## KNN补全缺失值
library(DMwR2)
CELL_fate_blot_velocity_174 <- knnImputation(CELL_fate_blot_velocity_174)


##标准化

CELL_fate_blot_velocity_174_scaled <- scale(t(CELL_fate_blot_velocity_174),center = T,scale = T)
CELL_fate_blot_velocity_174_scaled <- t(CELL_fate_blot_velocity_174_scaled)
CELL_fate_blot_velocity_174_scaled <- as.data.frame(CELL_fate_blot_velocity_174_scaled)


##聚类
CELL_fate_blot_velocity_174_scaled_dist <- dist(t(CELL_fate_blot_velocity_174_scaled), method = 'euclidean')
hclust_CELL_fate_blot_velocity_174_scaled <- hclust(CELL_fate_blot_velocity_174_scaled_dist,method = "ward.D2")
hclust_CELL_fate_blot_velocity_174_scaled_result <- cutree(hclust_CELL_fate_blot_velocity_174_scaled,k=6)
randIndex(table(hclust_CELL_fate_blot_velocity_174_scaled_result,CELLFateOriginalFactor), original = T)

CELL_fate_blot_velocity_174_dist <- dist(t(CELL_fate_blot_velocity_174), method = 'euclidean')
hclust_CELL_fate_blot_velocity_174 <- hclust(CELL_fate_blot_velocity_174_dist,method = "ward.D2")
hclust_CELL_fate_blot_velocity_174_result <- cutree(hclust_CELL_fate_blot_velocity_174,k=8)
library(flexclust)
randIndex(table(hclust_CELL_fate_blot_velocity_174_result,CELLFateOriginalFactor), original = T)

##欧式距离热图
pheatmap(as.matrix(CELL_fate_blot_velocity_174_scaled_dist), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)
pheatmap(as.matrix(CELL_fate_blot_velocity_174_dist), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)

CELL_fate_blot_velocity_174_class<-list(6)
for(i in 1:6){
  CELL_fate_blot_velocity_174_class[[i]] <-  CELL_fate_names[hclust_CELL_fate_blot_velocity_174_scaled_result == i]
  
}


#计算基因表达速率的矩阵174*706，与上面一样的操作
#

CELL_blot_velocity_174 <- matrix(nrow = 174,ncol = 706)
for(i in 1:174){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  for(j in 1:706){
    index <- which(realdata[[ind]][[1]] %in% CELL[j])
    if(length(index) == 0){
    }
    else if(length(index) == 1 & realdata[[ind]][[3]][index] >= 0){
      CELL_blot_velocity_174[i,j] <- realdata[[ind]][[3]][index] 
    }
    else if(length(index) == 1 & realdata[[ind]][[3]][index] < 0){
      CELL_blot_velocity_174[i,j] <- 0 
    }
    else{
      d <- diff(realdata[[ind]][[3]][index])
      d <- d[d>=0]
      CELL_blot_velocity_174[i,j] <- mean(d)
    }
  }
}

rownames(CELL_blot_velocity_174) <- GENE_names
colnames(CELL_blot_velocity_174) <- CELL
CELL_blot_velocity_174 <- as.data.frame(CELL_blot_velocity_174)

library(DescTools)
par(ps = 5)
PlotMiss(t(CELL_blot_velocity_174))
length(which(is.na(CELL_blot_velocity_174)))
library(DMwR2)
CELL_blot_velocity_174 <- knnImputation(CELL_blot_velocity_174)

CELL_blot_velocity_174_scaled <- scale(t(CELL_blot_velocity_174),center = T,scale = T)
CELL_blot_velocity_174_scaled <- t(CELL_blot_velocity_174_scaled)
CELL_blot_velocity_174_scaled <- as.data.frame(CELL_blot_velocity_174_scaled)

library(pheatmap)
pheatmap(as.matrix(dist(t(CELL_blot_velocity_174), method = 'euclidean')), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)
pheatmap(as.matrix(dist(t(CELL_blot_velocity_174_scaled), method = 'euclidean')), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)


##看相关系数矩阵热图

pheatmap(cor(t(as.matrix(CELL_blot_velocity_174)),method = "pearson"),cluster_rows = F,cluster_cols = F,show_colnames = F, fontsize = 5)
pheatmap(cor(t(as.matrix(CELL_blot_velocity_174)),method = "spearman"),cluster_rows = F,cluster_cols = F,show_colnames = F, fontsize = 5)


##把前代少于10个基因表达的细胞都删掉，实际上采用的是同时删除的，基因和细胞同时删除，没有先后顺序

CELL <- cell[1:706]
gene_missing_index <- c(13,17,18,58,62,63,79,97,146,153)


##细胞只取前706个（这是之前的做法），会漏掉一些初代的细胞
Data_3 <- data_3_all[-gene_missing_index,1:706]
rownames(Data_3) <- GENE_names

cell_notexpressed_index <- as.numeric(which(colSums(Data_3) <= 17))

Data_3 <- Data_3[,-cell_notexpressed_index]

CELL_Blot_weight_velocity_174_scaled <- CELL_blot_weight_velocity_174_scaled[,-cell_notexpressed_index]

CELL_old <- CELL
CELL <- CELL[-cell_notexpressed_index]

##最后的做法，基因一样删了10个，细胞只取前724个（用看的，这之前的细胞表达多，而且都没漏掉，可以做拟合）
##但是中间删了'EMS'这个细胞，即cell[511]，所以CELL是cell[1:510] + cell[512:725]

all(CELL == c(cell[1:510], cell[512:725]))

##表达矩阵也要改一下

Data_3 <- data_3_all[-gene_missing_index,c(1:510,512:725)]
rownames(Data_3) <- GENE_names



##少删除了一些缺失比例高的

cell_missing_index <- c(707:710, 721, 722)
n[cell_missing_index]
a <- cell[cell_missing_index]

#
#
##fate的表达矩阵好像错了，重新算(差在有的基因没有某些细胞的观测，比如pha-4_1,
##在Epla细胞表达了，但没有Eplaa这个细胞的观测，那么Eplaa这个细胞要算有表达还是没表达)
##现在把他当成没表达
all(cell == colnames(data_3_all))
fate_index <- NULL
for(i in 1:length(CELL)){
  fate_index <- c(fate_index, which(cell %in% CELL[i]))
}
all(CELL == cell[fate_index])
Data_3_simulation <- data_3_all[-gene_missing_index, fate_index]
Data_3_simulation <- Data_3_simulation[1:53,]
rownames(Data_3_simulation) <- GENE_names
colnames(Data_3_simulation) <- CELL
#
#
#


