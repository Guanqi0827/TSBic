
#改用原始数据，用均值

realdata_bind <- vector('list',106)
flag <- 0
for(i in 1:106){
  flag <- flag + group[i]
  index <- flag - group[i] + 1
  for(j in index:flag){
    ind <- which(file_names_csv%in%file_gene[j,1])
    realdata_bind[[i]] <- rbind(realdata_bind[[i]],realdata[[ind]])
  }
}

cell_blot <- matrix(nrow = 106,ncol = cell_number)
for(i in 1:106){
  for(j in 1:1083){
    index <- which(realdata_bind[[i]][[1]] %in% cell[j])
    if(length(index) == 0){
      cell_blot[i,j] <- 0 
    }
    else{
      cell_blot[i,j] <- mean(realdata_bind[[i]][[3]][index]) 
    }
  }
  
}
colnames(cell_blot) <- cell
rownames(cell_blot) <- gene_names_hu
cell_blot <- data.frame(cell_blot)

cell_fate_blot <- cell_blot[,cellFateOriginalIndex]
rownames(cell_fate_blot_184) <- rownames(cell_blot_184)


#画各基因文件相关系数热图

cell_blot_184_missing <- matrix(nrow = 184,ncol = cell_number)

for(i in 1:184){
  ind <- which(file_names_csv%in%file_gene[i,1])
  for(j in 1:1083){
    index <- which(realdata[[ind]][[1]] %in% cell[j])
    if(length(index) != 0){
      cell_blot_184_missing[i,j] <- mean(realdata[[ind]][[3]][index]) 
    }
  }
}

rownames(cell_blot_184_missing) <- gene_names_184
colnames(cell_blot_184_missing) <- cell
cell_blot_184_missing <- data.frame(cell_blot_184_missing)



cell_blot_184 <- data.frame(cell_blot_184)
colnames(cell_blot_184) <- cell
rownames(cell_blot_184) <- gene_names_184  


#缺失值可视化
library(DescTools)
length(which(is.na(cell_blot_184_missing)))
m <- NULL
for(i in 711:724){
  m <- append(m,184-length(which(is.na(cell_blot_184_missing[,i]))))
}

m <- NULL
for(i in 1:1083){
  m <- c(m, length(which(is.na(cell_blot_184_missing[-gene_missing_index,i])))/174)
}
min(m[726:1083])

##画基因细胞缺失比例的柱状图
m <- NULL
for(i in 1:184){
  m <- c(m, length(which(is.na(cell_blot_184_missing[i,])))/1083)
}

h <- hist(m, col = NULL, xlab = 'proportion', main = 'Histogram of gene missing proportion', breaks = seq(0, 1, 0.05))
h$density <- h$counts/sum(h$counts)
plot(h, freq=FALSE, xlab = 'proportion', main = 'Histogram of cell missing proportion')

n <- NULL
for(i in 1:1083){
  n <- c(n, length(which(is.na(cell_blot_184_missing[,i])))/184)
}
h <- hist(n, col = NULL, xlab = 'proportion', main = 'Histogram of cell missing proportion', breaks = seq(0, 1, 0.05))
h$density <- h$counts/sum(h$counts)
plot(h, freq=FALSE, xlab = 'proportion', main = 'Histogram of cell missing proportion')

m <- m[which(m >= 0.4)]
l <- which(n > 0.6)

length(which(is.na(cell_fate_blot_184)))
cell_fate_blot_184_missing <- cell_blot_184_missing[,cellFateOriginalIndex]

missing_index <- NULL

par(ps=3)
cell_blot_184_missing[,complete.cases(cell_blot_184_missing)]
PlotMiss(t(cell_blot_184_missing))
PlotMiss(t(cell_fate_blot_184_missing))

##好看点的图（多删那几个细胞）


PlotMiss(t(cell_blot_184_missing[-gene_missing_index, c(1:510,512:706,710:720,723:725)]))


##真实情况的图（我用的数据，少删了几个细胞）
PlotMiss(t(cell_blot_184_missing[-gene_missing_index, c(1:510,512:725)]))
PlotMiss(t(cell_blot_184_missing[-gene_missing_index, fate_index]))

axis(side = 1, at = c(0, 50, 100, 145))

plot(1:184,1:184,xlim = c(1,325), ylim = c(1,184), xaxt = 'n')


plot(0,0,type = 'n', xlim = c(1,325), ylim = c(1,184))
for(j in 1:325){
  points(j,184-sum(is.na(cell_fate_blot_184_missing[,j])),pch = 4, cex = 0.5) 
  if(184-sum(is.na(cell_fate_blot_184_missing[,j])) < 10){
    missing_index <- append(missing_index,j)
  }
  
}


#325个细胞blot热图
pheatmap(cell_fate_blot,cluster_rows = F,cluster_cols = F,show_colnames = F, show_rownames = F)


#pha-4基因的第一个文件，在E中表达的细胞里的表达速度

#挑选出数据，在一个细胞表达后，就默认一直表达

setwd("/Users/yanxianzhong/Desktop/基因聚类/hujiedata" )
which(file_newnames%in%file_gene[9,1])
which(file_names_csv%in%file_gene[9,1])
pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[9,1])]]

pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[9,1])]]
pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])

for(i in 1:dim(pha4_1_onset)[1]){
  pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
  
}


#画有表达的细胞的blot值散点图
library(RColorBrewer)
display.brewer.pal(n = 4, name = 'RdBu')
colorName <- brewer.pal(n = 4, name = 'RdBu')  

#k，第几个基因文件
#m，起始细胞
#n，终止细胞
#x，x轴起点
#y，x轴终点
DrawExpressedBlot <- function(k=1,m,n,x=20,y=140){
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(pha4_1_onset_realdata[[i]][[3]])
    pha4_1_onset_minblot[i-n+1] <- min(pha4_1_onset_realdata[[i]][[3]])
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'blot')
  for(i in m:n){
    d <- diff(pha4_1_onset_realdata[[i]]$blot)
    dtime <- which(d>=0)
    points(pha4_1_onset_realdata[[i]]$time[dtime+1],pha4_1_onset_realdata[[i]]$blot[dtime+1],pch = 4,col = colorName[i-m+1])
    
  }
  
  legend('topleft',pch = 4,col = colorName, ncol=4, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}
DrawExpressedBlot(51,1,4,45,150)


#画有表达的细胞的blot值的增量的绝对值散点图

DrawExpressedDiffBlot <- function(k=1,m,n,x=20,y=140){
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(diff(pha4_1_onset_realdata[[i]][[3]]))
    pha4_1_onset_minblot[i-n+1] <- min(diff(pha4_1_onset_realdata[[i]][[3]]))
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'diff(blot)')
  for(i in m:n){
    d <- diff(pha4_1_onset_realdata[[i]]$blot)
    dtime <- which(d>=0)
    d <- d[d>=0]
    points(pha4_1_onset_realdata[[i]]$time[-1][dtime+1],d,pch = 4,col = colorName[i-m+1])
  }
  
  legend('topleft',pch = 4,col = colorName, ncol=1, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}

DrawExpressedDiffBlot(51,1,4,45,150)





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


##
##删除后的矩阵为174*706
CELL <- cell[1:706]
gene_missing_index <- c(13,17,18,58,62,63,79,97,146,153)
GENE_names <- gene_names_184[-gene_missing_index]
CELL_blot_174_missing <- cell_blot_184_missing[-gene_missing_index,1:706]

CELLFateOriginalIndex <- cellFateOriginalIndex[cellFateOriginalIndex <= 706]
CELL_fate_names <- colnames(cell_fate_blot_184)[which(colnames(cell_fate_blot_184) %in% CELL)]
CELL[CELLFateOriginalIndex] == CELL_fate_names
CELL_fate_blot_174_missing <- cell_fate_blot_184_missing[-gene_missing_index,CELLFateOriginalIndex]
library(DescTools)

par(ps=5)
PlotMiss(t(CELL_blot_174_missing))
PlotMiss(t(CELL_fate_blot_174_missing))

length(which(is.na(CELL_blot_174_missing)))/(174*706)
length(which(is.na(CELL_fate_blot_174_missing)))/(174*145)


##KNN补全缺失值

library(DMwR2)
CELL_blot_174 <- knnImputation(CELL_blot_174_missing)
length(which(is.na(CELL_blot_174)))
CELL_fate_blot_174 <- knnImputation(CELL_fate_blot_174_missing)
length(which(is.na(CELL_fate_blot_174)))


##画相关系数热图,应该用174*724（用到了缺失值补全的版本）
library(DMwR2)
CELL_blot_174_missing <- cell_blot_184_missing[-gene_missing_index,c(1:510,512:725)]
CELL_blot_174_missing <- knnImputation(CELL_blot_174_missing)
pheatmap(cor(t(as.matrix(CELL_blot_174_missing)),method = "pearson"),cluster_rows = F,cluster_cols = F,show_colnames = F, fontsize = 5)


##画相关系数热图,应该用174*724（最新版本）（不用缺失值补全，就算能对齐的）
pheatmap(cor(t(as.matrix(CELL_blot_174_missing)),method = "pearson", use = 'pairwise.complete.obs'),cluster_rows = F,cluster_cols = F,show_colnames = F, fontsize = 5)




CELL_blot_174_scaled <- scale(t(CELL_blot_174),center = T,scale = T)
CELL_blot_174_scaled <- t(CELL_blot_174_scaled)

library(pheatmap)
pheatmap(cor(t(as.matrix(CELL_blot_174_scaled)),method = "pearson"),cluster_rows = F,cluster_cols = F,show_colnames = F, fontsize = 4)

a <- cor(t(as.matrix(CELL_blot_174_scaled)),method = "spearman")

##聚类算兰德指数
library(flexclust)
CELL_fate_blot_174_dist <- dist(t(CELL_fate_blot_174), method = 'euclidean')
hclust_CELL_fate_blot_174 <- hclust(CELL_fate_blot_174_dist,method = "ward.D2")
hclust_CELL_fate_blot_174_result <- cutree(hclust_CELL_fate_blot_174,k=6)
CELLFateOriginalFactor <- as.factor(cell_fate_original_hu[colnames(cell_fate_blot_184) %in% colnames(CELL_fate_blot_174),2])
randIndex(table(hclust_CELL_fate_blot_174_result,CELLFateOriginalFactor), original = T)

CELL_fate_blot_174_scaled <- scale(t(CELL_fate_blot_174),center = T,scale = T)
CELL_fate_blot_174_scaled <- t(CELL_fate_blot_174_scaled)
CELL_fate_blot_174_scaled_dist <- dist(t(CELL_fate_blot_174_scaled), method = 'euclidean')


#####画欧式距离图（用到了缺失值补全的版本）
library(pheatmap)
CELL_fate_blot_174_scaled <- scale(t(CELL_fate_blot_174),center = T,scale = T)
CELL_fate_blot_174_scaled <- t(CELL_fate_blot_174_scaled)
CELL_fate_blot_174_scaled_dist <- dist(t(CELL_fate_blot_174_scaled), method = 'euclidean')
pheatmap(as.matrix(CELL_fate_blot_174_scaled_dist), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)

#####画欧式距离图（最新版本）（不用缺失值补全，只算能对齐的）（实际上dist函数就是这样算的，支持NA值）
library(pheatmap)
CELL_fate_blot_174_missing <- cell_blot_184_missing[-gene_missing_index,CELLFateOriginalIndex]
CELL_fate_blot_174 <- CELL_fate_blot_174_missing
CELL_fate_blot_174_scaled <- scale(t(CELL_fate_blot_174),center = T,scale = T)
CELL_fate_blot_174_scaled <- t(CELL_fate_blot_174_scaled)
CELL_fate_blot_174_scaled_dist <- dist(t(CELL_fate_blot_174_scaled), method = 'euclidean')
pheatmap(as.matrix(CELL_fate_blot_174_scaled_dist), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)



#####画欧式距离图（最新版本）（细胞排序版本）

a <- sort(CELL)
a <- c(a[59:145], a[1:58])
index <- c()
for(j in 1:145){
  index <- c(index, which(CELL %in% a[j]))
}
df <- CELL_fate_blot_174_scaled[,index]
df_dist <- dist(t(df), method = 'euclidean')
pheatmap(as.matrix(df_dist), cluster_rows = F, cluster_cols = F , show_colnames = F, fontsize = 5)













