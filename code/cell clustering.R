#读取原始数据，胡老师数据，基因名称数据
getwd()
setwd("/Users/yanxianzhong/Desktop/GeneClustering/hujiedata" )
file_names <- list.files(pattern="*.xls") 
library(readxl)
hujiedata <- lapply(file_names, function(x) read_xls(x,sheet = 3)) 

setwd("/Users/yanxianzhong/Desktop/GeneClustering/realdata" )
file_names_csv <- list.files(pattern="*.csv")
realdata <- lapply(file_names_csv, function(x) read.csv(x, header=T,stringsAsFactors=F)) 
realdata <- lapply(realdata, function(x) x=x[,c(2,3,7)])

setwd("/Users/yanxianzhong/Desktop/GeneClustering/gene_names" )
promotergenes <- read.csv('promotergenes.csv', header=T)
gene_names <- promotergenes[,1]

promoter <- read.csv("promoter.csv",header=TRUE)


#数据集去引号
ReplaceStr<-function(x){
  y<-x
  y$cell<-gsub("\"", '', y$cell)
  y$time<-gsub("\"", '', y$time)
  y$number.of.cells.when.onset<-gsub("\"", '', y$number.of.cells.when.onset)
  return(as.data.frame(y))
}

hujiedata <- lapply(hujiedata, function(x) ReplaceStr(x)) 

#文件名和基因名称对应

group <- promotergenes$NumberOfDuplicates
gene_names_rep <- NULL
for(i in 1:length(group)){
  gene_names_rep <- c(gene_names_rep,rep(gene_names[i],group[i]))
}

file_gene_ <- cbind(promoter[,1],gene_names_rep,promoter[,3])
colnames(file_gene_)[1] <- colnames(promoter)[1]
colnames(file_gene_)[3] <- colnames(promoter)[3]

#选出胡老师数据
library(dplyr)
file_newnames <- NULL
for(i in 1:length(file_names)){
  a <- sub('xls','csv',file_names[i])
  file_newnames <- c(file_newnames,a)
}
file_gene <- filter(as.data.frame(file_gene_), file_gene_[,1] %in% file_newnames)
group_hu <- table(unlist(file_gene$gene_names_rep))
group_hu <- as.data.frame(group_hu)
a <- rep(0,dim(file_gene)[1])
for(i in 1:length(group_hu_sorted_freq)){
  index <- which(file_gene$gene_names_rep == group_hu[i,1])
  a[index] <- c(group_hu[i,2],rep(NA,length(index)-1))
}
file_gene$Group <- a
sum(na.omit(file_gene$Group))
group_ <- group
group <- na.omit(file_gene$Group)

#选出胡老师数据里的所有细胞

realdata_hu <- list(0)
for(i in 1:dim(file_gene)[1]){
  realdata_hu[i] <- realdata[which(file_names_csv == file_gene[i,1])]
}

realdata_hu_cell <- list(0)
for(i in 1:length(realdata_hu)){
  realdata_hu_cell[[i]] <- realdata_hu[[i]][,1]
  
}

cell <- Reduce(union,realdata_hu_cell)
cell_number <- length(cell)

#求每个基因表达的细胞名称
file_onset <- list(0)
for(i in 1:dim(file_gene)[1]){
  index <- which(file_newnames == file_gene[i,1])
  file_onset[[i]] <- hujiedata[[index]]$cell
  
}


cell_express <- vector('list',length(file_onset))
for(i in 1:length(file_onset)){
  for(j in 1:length(file_onset[[i]])){
    cell_express[[i]] <- c(cell_express[[i]],list(cell[grep(file_onset[[i]][[j]],cell)]))
  }
}

cell_express_ <- vector('list',length(cell_express))
for(i in 1:length(cell_express)){
  for(j in 1:length(cell_express[[i]])){
    cell_express_[[i]] <- c(cell_express_[[i]],cell_express[[i]][[j]])
  }
}

cell_express_all <- vector('list',length(group))
flag <- 0 
for(i in 1:length(group)){
  flag <- flag + group[i] 
  index <- flag - group[i] + 1
  for(j in index:flag){
    cell_express_all[[i]] <- unique(c(cell_express_all[[i]],cell_express_[[j]]))
  }
}

#0-1化
cell_express_0_1 <- vector('list',length(group))
for(i in 1:length(group)){
  cell_express_0_1[[i]] <- rep(0,cell_number)
  cell_express_0_1[[i]][cell %in% cell_express_all[[i]]] <-1
  
}

gene_names_hu <- file_gene$gene_names_rep
gene_names_hu <- unique(gene_names_hu)

data_3 <- matrix(nrow = length(cell_express_0_1),ncol = cell_number)
for(i in 1:length(cell_express_0_1)){
  data_3[i,] <- cell_express_0_1[[i]]
  
}
data_3 <- data.frame(data_3)
colnames(data_3) <- cell
rownames(data_3) <- gene_names_hu

# binary距离做聚类
num_fate <- length(tissue_name)
num_class <- num_fate
cell_dist <- dist(t(data_3), method = 'binary')
hclust_ward <- hclust(cell_dist,method = "ward.D2")
hclust_ward_result <- cutree(hclust_ward,k=10)
cell_class<-list(0)
for(i in 1:10){
  cell_class[[i]] <- cell[hclust_ward_result == i]
  
}
table(hclust_ward_result)

#相关系数,带阈值的高斯核函数,k-Means方法尝试细胞命运比例聚类（效果都不好）

# r <- NULL
# wss <- NULL
# for (i in 2:19){
#   a <- kmeans(cell_fate_proportion, centers=i)
#   wss[i] <- sum(a$withinss)
#   r[i] <- randIndex(table(a$cluster,hclust_ward_result))
# }
# wss <- wss[-1]
# r <- r[-1]
# plot(2:19,r,type='b')
# plot(2:19,wss,type='b')
# text(3,wss[1],p)
# p <- paste('(2,500)')
library(corrplot)
# cell_pearson <- cor(t(cell_fate_proportion), method = 'pearson')
# c <- 1-cell_pearson
# dd <- as.dist(diss)
# a <- hclust(dd,method = "ward.D2")
# b <- cutree(a,k=10)
# randIndex(table(b,hclust_ward_result))
# 
# d <- as.matrix(cell_fate_proportion_dist)
# v <- as.vector(cell_pearson)
# a <- var(v)
# diss <- matrix(ncol = 1083,nrow = 1083)
# 
# for(i in 1:1083){
#   for(j in i:1083){
#     if(cell_pearson[i,j] >= 0.1){
#       diss[i,j] <- -a^2*log(cell_pearson[i,j])}
#     else{diss[i,j] <- 1.5}
#   }
# 
# }
# diss[lower.tri(diss)] <- t(diss)[lower.tri(diss)]

#K-Modes对细胞聚类
# library(klaR)
# a <- kmodes(t(data_3),modes = c)
# b <- NULL
# for(i in 1:12){
#   b <- append(b,cell_class[[i]][[1]])
# }
# index <- which(cell %in% b)
# c <- t(data_3)[index,]
# r <- a$cluster
# 
# randIndex(table(hclust_fate_proportion_result,r))

# cosine距离做聚类
# library(philentropy)
# 
# a <- distance(t(data_3),method = 'cosine')
# a[is.nan(a)] <- 0
# b <- hclust(as.dist(a),method = "ward.D2")
# c <- cutree(b,k=12)
# table(c)
# randIndex(table(hclust_fate_proportion_result,c))

#画图

plot(hclust_ward,hang = -1,labels = F)
rect.hclust(hclust_ward,k=9)

library(sparcl)
ColorDendrogram(hclust_ward,hclust_ward_result,branchlength = 20)

hc_ward = as.dendrogram(hclust_ward) 
plot(cut(hc_ward,h = 3)$upper,main = "Upper tree of cut at h=3")
plot(cut(hc_ward,h = 2)$lower[[9]], main = "Ninth branch of lower tree with cut at h=2")


# 读取细胞命运

setwd("/Users/yanxianzhong/Desktop/基因聚类/cell_fate" )
cell_fate <- read.csv('allcellfate.csv',header = F)
colnames(cell_fate)[1] <- 'cell'
colnames(cell_fate)[2] <- 'tissue'

cell_fate_hu <- data.frame(matrix(nrow = length(cell),ncol = 2))
colnames(cell_fate_hu)[1] <- 'cell'
colnames(cell_fate_hu)[2] <- 'tissue'
for(i in 1:length(cell)){
  if(cell[i] %in% cell_fate[,1]){
    cell_fate_hu[i,] <- cell_fate[which(cell_fate[,1] %in% cell[i]),]
  }
  
}


# 计算各细胞的细胞命运

setwd("/Users/yanxianzhong/Desktop/基因聚类/cell_fate" )
cell_fate_original <- read.csv('cellfate.csv',header = F)
colnames(cell_fate_original)[1] <- 'cell'
colnames(cell_fate_original)[2] <- 'tissue1'
colnames(cell_fate_original)[3] <- 'tissue2'
tissue_name <- Reduce(union,cell_fate_original[,2],cell_fate_original[,3])

cell_fate_hu_0_1 <- data.frame(matrix(0,nrow = length(cell),ncol = length(tissue_name)))
colnames(cell_fate_hu_0_1)[1:length(tissue_name)] <-  tissue_name
rownames(cell_fate_hu_0_1)[1:length(cell)] <- cell

for(i in 1:cell_number){
  for(j in 1:length(tissue_name)){
    index <- grep(tissue_name[j],cell_fate_hu[i,2])
    if(length(index)!=0){
      cell_fate_hu_0_1[i,j] <- 1
    }
  }
}


#只聚类原始细胞命运数据中的细胞(只有325个细胞在胡老师的数据里)
cellFateOriginalIndex <- NULL

for(i in 1:length(cell_fate_original[,1])){
  index <- which(cell == cell_fate_original[i,1])
  cellFateOriginalIndex <- append(cellFateOriginalIndex,index)
  
}

cellFateOriginalName <- cell[cellFateOriginalIndex]

cellFateHuIndex <- NULL
for(i in 1:length(cellFateOriginalIndex)){
  index <- which(cell_fate_original[,1]==cellFateOriginalName[i])
  cellFateHuIndex <- append(cellFateHuIndex,index)
}


cellFateOriginalData <- data_3[,cellFateOriginalIndex]
cellFateOriginalDataDist <- dist(t(cellFateOriginalData), method = 'binary')
hclustOriginalDataWard <- hclust(cellFateOriginalDataDist,method = "ward.D2")
hclustOriginalDataWardResult <- cutree(hclustOriginalDataWard,k=8)

#聚类树状图
plot(hclustOriginalDataWard,hang = -1,labels = F)
rect.hclust(hclustOriginalDataWard,k=8)
mtext(c('MS','ABa','ABp','E','D','C','ABp','ABa'),side = 1,at=c(25,80,150,194,207,225,262,305))

#带颜色带树状图
hclust_OriginalDataWard = as.dendrogram(hclustOriginalDataWard) 
library(factoextra)

library(RColorBrewer)
display.brewer.all(type = "qual")
display.brewer.pal(name = "Set3", n = 12)

color1 <- rgb(43,85,125, maxColorValue = 255)
color2 <- rgb(69,189,155, maxColorValue = 255)
barplot(rep(1,2),1,col = c(color1, color2))


colorName <- brewer.pal(9,"Set1")
display.brewer.pal(9,"Set1")
colorName <- brewer.pal(name = "Set3", n = 12)[c(4,6,7,8,3,1,5,10)]

fig <- fviz_dend(hclust_OriginalDataWard, k = 8, cex = 1, show_labels = F,
                 rect = T, lwd = 1, xlab = 'Cell', main = '')

fig <- fviz_dend(hclust_OriginalDataWard, k = 8, palette = hcl.colors(8, "Set 3"), cex = 1, show_labels = F,
                 rect = T, lwd = 1, xlab = 'Cell', main = '')
fviz_dend(hclust_OriginalDataWard, k = 8, palette = brewer.pal(9,"Set1")[-6], cex = 0.3, show_labels = F,
          rect = T, lwd = 0.5, xlab = 'Cell', main = '')

fviz_dend(hclust_OriginalDataWard[[2]][[2]][[1]][[2]], palette = hcl.colors(8, "Set 3")[7], cex = 0.7,
          rect = T, lwd = 0.5, xlab = 'Cell', main = '')

##数据框不能有重名的，打两次标签
mat <- data.frame('x' = c(25,194,207,225,262,305),
                  'y' = rep(-0.5, 6))
labelName <- c('MS','E','D','C','ABp','ABa')
rownames(mat) <- labelName
fig <- fviz_add(fig, mat, addlabel = T, geom = 'text', labelsize = 6)

mat <- data.frame('x' = c(80,150),
                  'y' = rep(-0.5, 2))
labelName <- c('ABa','ABp')
rownames(mat) <- labelName
fviz_add(fig, mat, addlabel = T, geom = 'text', labelsize = 6)

hclust_OriginalDataWard = as.dendrogram(hclustOriginalDataWard) 
plot(cut(hclust_OriginalDataWard,h = 1.5)$upper,main = "Upper tree of cut at h=1.5")
plot(cut(hclust_OriginalDataWard,h = 1.5)$lower[[6]], main = "Sixth branch of lower tree with cut at h=1.5")

#兰德指数

cellFateOriginalFactor <- as.factor(cell_fate_original[cellFateHuIndex,2])
randIndex(table(cellFateOriginalFactor,hclustOriginalDataWardResult),original = T) 
randIndex(table(cellFateOriginalFactor,hclust_fate_proportion_result[cellFateOriginalIndex]),original = T)

for(i in 1:length(hclustMethod)){
  for(j in 1:19){
    hclustDifferentMethod <- hclust(cellFateOriginalDataDist,method = hclustMethod[i])
    hclustDifferentMethodResult <- cutree(hclustDifferentMethod,k=j)
    tableDifferentMethodResult <- table(cellFateOriginalFactor,hclustDifferentMethodResult)
    hclust_randindex[j,i] <- randIndex(tableDifferentMethodResult)
  }
  
}

plot(1:19,hclust_randindex[,5],type='b',main='Adjusted Rand Index of Clusters',ylab='Adjusted Rand Index',xlab='Clusters',xaxt = 'n')
axis(side = 1,at=1:19,tick = 1:19)
abline(v=8,lty=2)

cell_class<-list(0)
for(i in 1:8){
  cell_class[[i]] <- cellFateOriginalName[hclustOriginalDataWardResult == i]
  
}

randIndex(table(cellFateOriginalFactor[cell_class_index],hclustOriginalDataWardResult[cell_class_index]),original = T) 

#画各类细胞组织的饼图

library(RColorBrewer)
display.brewer.pal(n = 8, name = 'RdBu')

for(i in 1:8){
  cellClassFate <- NULL
  cellClassFateIndex <- NULL
  for(j in 1:length(cell_class[[i]])){
    index <- which(cell_fate_original[,1] == cell_class[[i]][[j]])
    cellClassFateIndex <- append(cellClassFateIndex,index)
    cellClassFate <- append(cellClassFate,cell_fate_original[index,2])
    
  }
  percent <- round(100*table(cellClassFate)/length(cell_class[[i]]),2)
  percent <- paste(percent, "%", sep = "")
  plotName <- paste('Class',i)
  labelName <- NULL
  for(k in 1:length(table(cellClassFate))){
    labelName <- append(labelName,paste(names(table(cellClassFate))[k],':',
                                        as.numeric(table(cellClassFate)[k]),'/',
                                        as.numeric(cellFateNum[names(cellFateNum)==names(table(cellClassFate))[k]]),
                                        '(',round(cell_enrichment_result[[i]][[k]],3),')'))
  }
  if(length(table(cellClassFate))>=3){
    colorName <- brewer.pal(n = length(table(cellClassFate)), name = 'RdBu')
  }
  else{
    colorName <- brewer.pal(n = 4, name = 'RdBu')  
  }
  
  pdf(file = plotName, width = 17)
  pie(table(cellClassFate), labels = labelName, main = paste(plotName,'(number of cells:',length(cell_class[[i]]),')'),cex = 0.8,radius = 1,col = colorName)
  legend("bottomleft",percent, cex = 0.9, fill = colorName)
  dev.off()
}


#用于查看各类中各细胞名称及命运

cellClassFate <- NULL
cellClassFateIndex <- NULL
i <- 6
for(j in 1:length(cell_class[[i]])){
  index <- which(cell_fate_original[,1] == cell_class[[i]][[j]])
  cellClassFateIndex <- append(cellClassFateIndex,index)
  cellClassFate <- append(cellClassFate,cell_fate_original[index,2])
  
}


tissueName <- 'Nerve Tissue'
cellName <- cell_class[[i]]
cellName <- cell_class[[i]][cellClassFate == tissueName]
cellClassTissueData <- t(cellFateOriginalData)[rownames(t(cellFateOriginalData)) %in% cellName,]

a <- colSums(data_3)
which(a == 0)

# 细胞富集分析

options(digits = 7)
cell_fate_original_hu <- cell_fate_original[cellFateHuIndex,1:2]
cellFateNum <- table(cell_fate_original_hu[,2])


cellEnrichment <- function(allcell, cellclass){
  cellnum <- length(allcell)
  classcellnum <- length(cellclass)
  cellClassFate <- NULL
  cellClassFateIndex <- NULL
  for(j in 1:length(cellclass)){
    index <- which(cell_fate_original_hu[,1] == cellclass[[j]])
    cellClassFateIndex <- append(cellClassFateIndex,index)
    cellClassFate <- append(cellClassFate,cell_fate_original_hu[index,2])
    
  }
  pvalue <- NULL
  for(i in 1:length(table(cellClassFate))){
    classfatenum <- as.numeric(table(cellClassFate)[i])
    fatenum <- as.numeric(cellFateNum[names(table(cellClassFate))[i]])
    pvalue <- append(pvalue,1-phyper(classfatenum-1, fatenum, cellnum-fatenum,classcellnum))
  }
  pvalueadj <- p.adjust(pvalue,method = 'BH',n = length(pvalue))
  names(pvalueadj) <- names(table(cellClassFate))
  return(pvalueadj)
}





cell_enrichment_result <- vector('list',8)
for(i in 1:length(cell_class)){
  cell_enrichment_result[[i]] <- cellEnrichment(allcell =  cellFateOriginalName, cellclass = cell_class[[i]])
  
}


#找哪几类细胞分得不好

randIndex(table(hclust_fate_proportion_result[cellFateLastIndex],hclust_ward_result[cellFateLastIndex]),original = T)
randIndex(table(hclust_fate_proportion_result[cellFateOriginalIndex],hclust_ward_result[cellFateOriginalIndex]),original = T)

cellClassChosen <- c(1,2)
cell_class_index <- NULL
for(i in 1:length(cellClassChosen)){
  classIndex <- cellClassChosen[i]
  for(j in 1:length(cell_class[[classIndex]])){
    index <- which(cell == cell_class[[classIndex]][[j]])
    cell_class_index <- append(cell_class_index,index)
  }
}



# 确定聚类类别数

library(NbClust)

nbclust_fate_proportion <- NbClust(cell_fate_proportion, distance = "euclidean", min.nc = 2, max.nc = 19, method = "ward.D2")
barplot(table(nbclust_fate_proportion$Best.nc[1,]),xlab='Number of Clusters',ylab='Number of Criteria',main='Number of Clusters Chosen by 26 Criteria')
dev.off()
library(factoextra)
library(ggplot2)
fviz_nbclust(cell_fate_proportion,hcut, method='wss',dist = 'euclidean',k.max = 19)+geom_vline(xintercept = 12 ,linetype = 2)
dev.off()
fviz_nbclust(t(data_3),hcut, method='wss', dist = 'binary',k.max = 19)+geom_vline(xintercept = 12 ,linetype = 2)
fviz_nbclust(cell_fate_proportion,hcut, method='silhouette',dist = 'euclidean',k.max = 19)+geom_vline(xintercept = 12 ,linetype = 2)
fviz_nbclust(t(data_3),hcut, method='silhouette', dist = 'binary',k.max = 19)+geom_hline(yintercept = 0.41,linetype = 2)

library(mclust)

mclust_fate_proportion <- Mclust(cell_fate_proportion, G = 2:19)
summary(mclust_fate_proportion)
plot.Mclust(mclust_fate_proportion, what = 'BIC')


# 计算不同聚类方法不同聚类类别数的兰德指数

library(flexclust)

hclustMethod <- c('single','complete','centroid','average','ward.D2')
hclust_randindex <- data.frame(matrix(nrow = 19,ncol = length(hclustMethod)))
colnames(hclust_randindex) <- hclustMethod


for(i in 1:length(hclustMethod)){
  for(j in 1:19){
    hclustDifferentMethod <- hclust(cell_dist,method = hclustMethod[i])
    hclustDifferentMethodResult <- cutree(hclustDifferentMethod,k=j)
    hclustDifferentMethodFateResult <- cutree(hclust_fate_ward,k=j)
    tableDifferentMethodResult <- table(hclustDifferentMethodFateResult,hclustDifferentMethodResult)
    hclust_randindex[j,i] <- randIndex(tableDifferentMethodResult)
  }
  
}


#不把基因当成一组，当成全部184组测量结果，来做聚类
##新版做法，验证了一下，和之前算的是一样的
##这里要注意，要看基因中是否有该细胞，没有的话直接按不表达处理，有该细胞且表达，才算表达
realdata_hu_cell <- list(0)
for(i in 1:length(realdata_hu)){
  realdata_hu_cell[[i]] <- realdata_hu[[i]][,1]
  realdata_hu_cell[[i]] <- unique(realdata_hu_cell[[i]])
  
}

a <- matrix(0, nrow = length(cell_express_0_1_all),ncol = cell_number)
for(i in 1:length(cell_express_0_1_all)){
  for(j in 1:1083){
    if(cell[j] %in% realdata_hu_cell[[i]]){
      if(cell_express_0_1_all[[i]][j] == 1){
        a[i, j] = 1
      }
    }
  }
}
a <- data.frame(a)
all(a == data_3_all)

#####不把基因当成一组，当成全部184组测量结果，来做聚类
gene_names_hu_rep <- file_gene[,2]
cell_express_0_1_all <- vector('list',length(gene_names_hu_rep))
for(i in 1:length(gene_names_hu_rep)){
  cell_express_0_1_all[[i]] <- rep(0,cell_number)
  cell_express_0_1_all[[i]][cell %in% cell_express_[[i]]] <-1
  
}##注意，这样处理不行，可能有的基因没有某些表达的细胞


data_3_all <- matrix(nrow = length(cell_express_0_1_all),ncol = cell_number)
for(i in 1:length(cell_express_0_1_all)){
  data_3_all[i,] <- cell_express_0_1_all[[i]]
  
}

data_3_all <- data.frame(data_3_all)
colnames(data_3_all) <- cell
###

cell_dist_all <- dist(t(data_3_all), method = 'binary')
hclust_ward_all <- hclust(cell_dist_all,method = "ward.D2")
hclust_ward_result_all <- cutree(hclust_ward_all,k=num_class)

library(flexclust)
randIndex(table(hclust_ward_result_all,hclust_ward_result))

table(hclust_fate_ward_result,hclust_ward_result_all)
randIndex(table(hclust_fate_ward_result,hclust_ward_result_all))



library(flexclust)
#全部325个细胞基因聚类，选择聚类数
gene_dist <- dist(cellFateOriginalData, method = 'binary')
hclust_gene <- hclust(gene_dist,method = "ward.D2")
hclust_gene_result <- cutree(hclust_gene,k=6)
table(hclust_gene_result)
pdf(file = 'hclust_gene', width = 17)
plot(hclust_gene,hang = -1,labels = F)
rect.hclust(hclust_gene,k=6)
dev.off()


##画好看点的树状图
library(factoextra)
library(RColorBrewer)

hclust_gene_dend <- as.dendrogram(hclust_gene) 

fviz_dend(hclust_gene_dend , k = 6, cex = 0.8, label_cols = 'black', 
          rect = T, lwd = 0.5, xlab = 'Gene', main = '')


hclust_class_gene_dend <- as.dendrogram(hclust_class_gene) 

fviz_dend(hclust_class_gene_dend, k = 6, cex = 0.8, label_cols = 'black', 
          rect = T, lwd = 0.5, xlab = 'Gene', main = '')


library(ggplot2)

fviz_nbclust(data_3,hcut, method='wss', dist = 'binary',k.max = 10)
fviz_nbclust(data_3,hcut, method='silhouette', dist = 'binary',k.max = 10)

library(mclust)

mclust_gene <- Mclust(data_3, G = 2:10)
summary(mclust_gene)
plot.Mclust(mclust_gene, what = 'BIC')


#对每类细胞分别基因聚类,画聚类图

hclust_gene_randindex <- NULL
for(i in 1:8){
  cellClassFateIndexHu <- NULL
  plotName <- paste('hclust_gene_class',i,sep = "_")
  for(j in 1:length(cell_class[[i]])){
    index <- which(cell_fate_original_hu[,1] == cell_class[[i]][[j]])
    cellClassFateIndexHu <- append(cellClassFateIndexHu,index)
    
  }
  cell_class_data <- cellFateOriginalData[,cellClassFateIndexHu]
  gene_class_dist <- dist(cell_class_data, method = 'binary')
  hclust_class_gene <- hclust(gene_class_dist,method = "ward.D2")
  hclust_class_gene_result <- cutree(hclust_class_gene,k=6)
  hclust_gene_randindex <- append(hclust_gene_randindex,randIndex(table(hclust_gene_result,hclust_class_gene_result)))
  pdf(file = plotName, width = 17)
  plot(hclust_class_gene,hang = -1)
  rect.hclust(hclust_class_gene,k=6)
  dev.off()
  
}

a <- hclust_class_gene_result































