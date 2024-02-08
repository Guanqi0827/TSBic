##存全部细胞blot的系数向量和次数向量
CELL_all_blot_polySpl <- list()

for(i in 1:length(GENE_names)){
  li <- CELL_blot_polySpl[[i]]
  li_all <- list()
  for(j in 1:length(CELL)){
    l <- list()
    cell_name <- CELL[j]
    s <- substr(cell_name, 1, 1)
    id <- grep(s, names(li))
    if(length(id) == 1){
      lii <- li[[id]]
      ind <- which(names(lii[[1]]) %in% cell_name)
      b_index <- lii[[1]][[ind]]
      if(!is.na(b_index[1])){
        beta_h <- lii$beta_h[b_index]
        degree <- 0:(length(b_index) - 1)
        l <- list('coefficients' = beta_h, 'degree' = degree)
      }
      li_all[[cell_name]] <- l
    }
    else{
      li_all[[cell_name]] <- l
    }

  }
  CELL_all_blot_polySpl[[GENE_names[i]]] <- li_all
}

save(CELL_all_blot_polySpl, file = 'CELL_all_blot_polySpl.Rdata')


save(realdata_hu, file = 'realdata_hu.Rdata')



#######对全部模型求导
CELL_all_diff_polySpl_deri <- list()

for(i in 1:length(GENE_names)){
  li <- CELL_all_blot_polySpl[[i]]
  li_all <- list()
  for(j in 1:length(CELL)){
    cell_name <- CELL[j]
    l <- derivative(li[[j]])
    li_all[[cell_name]] <- l
  }
  CELL_all_diff_polySpl_deri[[GENE_names[i]]] <- li_all
}

save(CELL_all_diff_polySpl_deri, file = 'CELL_all_diff_polySpl_deri.Rdata')


#######对全部导数模型进行平移缩放

##

CELL_all_diff_polySpl_deri_scale <- list()

for(i in 1:length(GENE_names)){
  li_all <- list()
  for(j in 1:length(CELL)){
    mode <- CELL_all_diff_polySpl_deri[[i]][[j]]
    li <- deriv_scale(mode, i, CELL[j])
    li_all[[CELL[j]]] <- li
  }
  CELL_all_diff_polySpl_deri_scale[[GENE_names[i]]] <- li_all
}

save(CELL_all_diff_polySpl_deri_scale, file = 'CELL_all_diff_polySpl_deri_scale.Rdata')

##这里已经提前算好corr，直接读取即可
##写两组基因的corr，之前按基因下标定义，现在按基因名称定义
#gene_names是定义的112个基因名称，注意不要重复，这里统一用gene_name表示两个基因名称
#细胞照样用下标即可，method有三种：single（最短），complete（最长），average（平均）
two_gene_corr <- function(gene_name_1, gene_name_2, cell_index, method = 'average'){
  gene_index_1 <- GENE_group[[gene_name_1]]
  gene_index_2 <- GENE_group[[gene_name_2]]
  gtotal <- c()
  for(j in 1:length(gene_index_1)){
    g_ind_1 <- gene_index_1[j]
    for(k in 1:length(gene_index_2)){
      g_ind_2 <- gene_index_2[k]
      gcorr <- c()
      for(i in 1:length(cell_index)){
        ind <- cell_index[i]
        li <- all_corr[[ind]]
        if(g_ind_2 < g_ind_1){
          gcorr <- c(gcorr, li[g_ind_2, g_ind_1])
        }
        else{
          gcorr <- c(gcorr, li[g_ind_1, g_ind_2])
        }
      }
      gcorr <- as.numeric(na.omit(gcorr))
      if(length(gcorr) > 0){
        gcorr <- mean(gcorr)
        gtotal <- c(gtotal, gcorr)
      }
    }
  }
  gtotal <- as.numeric(na.omit(gtotal))
  if(length(gtotal) > 0){
    if(method == 'single'){
      gcorr_score <- min(gtotal)
    }
    if(method == 'average'){
      gcorr_score <- mean(gtotal)
    }
    if(method == 'complete'){
      gcorr_score <- max(gtotal)
    }
  }
  else{
    gcorr_score <- 0
  }
  return(gcorr_score)
}
  
grep('Eala', CELL)
two_gene_corr('pha-4', 'B0310.2', grep('E', CELL), method = 'average')
two_gene_corr('pha-4', 'B0310.2', grep('E', CELL), method = 'single')
two_gene_corr('pha-4', 'B0310.2', grep('E', CELL), method = 'complete')


##从基因名称集合转化为基因下标集合的函数
##输入为gene_name，输出为gene_index

Getgene_index <- function(gene_name){
  gene_index <- c()
  for(i in 1:length(gene_name)){
    gene_index <- c(gene_index, GENE_group[[gene_name[i]]])
  }
  return(gene_index)
}

Getgene_index(c('pha-4','tbx-11'))


##corr这项的函数

gene_corr <- function(gene_name, cell_index, method = 'average'){
  gtotal <- NULL
  for(j in 1:(length(gene_name) - 1)){
    g_name_1 <- gene_name[j]
    for(k in (j + 1):length(gene_name)){
      g_name_2 <- gene_name[k]
      gcorr <- two_gene_corr(g_name_1, g_name_2, cell_index, method)
      gtotal <- c(gtotal, gcorr)
    }
  }
  gtotal <- na.omit(gtotal)
  if(length(gtotal) > 0){
    gcorr_score <- mean(gtotal)
  }
  else{
    gcorr_score <- 0
  }
  gene_index <- Getgene_index(gene_name)
  back <- list('Gene Names'= gene_name,'Cell Names'=CELL[cell_index],
               'Gcorr' = gtotal,'Gcorr_Score'= gcorr_score)
  return(back)
  
}

gene_corr(c('pha-4', 'tbx-11', 'B0310.2'), grep('E', CELL))






###用p值表示细胞间距离，创建p值距离分函数(p值越大距离越近)

cell_pdist <- function(gene_name, cell_index){
  gene_index <- Getgene_index(gene_name)
  ptotal <- NULL
  for(j in 1:(length(cell_index)-1)){
    c_ind_1 <- cell_index[j]
    for(k in (j+1):length(cell_index)){
      c_ind_2 <- cell_index[k]
      pvalue <- NULL
      for(i in 1:length(gene_index)){
        ind <- gene_index[i]
        li <- all_pvalue[[ind]]
        if(c_ind_2 < c_ind_1){
          pvalue <- c(pvalue, li[c_ind_2, c_ind_1])
        }
        else{
          pvalue <- c(pvalue, li[c_ind_1, c_ind_2])
        }
      }
      pvalue <- as.numeric(na.omit(pvalue))
      if(length(pvalue) > 0){
        pvalue <- min(pvalue)
        ptotal <- c(ptotal, pvalue)
      }
    }
  }
  if(length(ptotal) > 0){
    pdist_score <- mean(ptotal)
  }
  else{
    pdist_score <- 0
  }
  back <- list('Gene Names'= gene_name,'Cell Names'=CELL[cell_index],
               'Pdist' = ptotal,'Pdist_Score'=pdist_score)
  return(back)
}

cell_pdist(c('pha-4','tbx-11'), sample(grep('E',CELL), 31))
cell_pdist(GENE_names[gene_index], c(cell_index))


##计算有效面积


area <- function(gene_name, cell_index){
  gene_index <- Getgene_index(gene_name)
  bicluster_expressed <- Data_3[gene_index,cell_index]
  area_ <- sum(bicluster_expressed)
  return(area_)
  
}

area(c('pha-4','tbx-11'), grep('E',CELL))

##score函数

score <- function(gene_name, cell_index, alpha=1, lambda = 1, beta = 1, method = 'average'){
  corr_gene <- gene_corr(gene_name, cell_index, method)
  corr_gene_all <- corr_gene$'Gcorr'
  corr_gene_s <- (corr_gene$'Gcorr_Score')^alpha
  dist_cell <- cell_pdist(gene_name, cell_index)
  dist_cell_all <- dist_cell$Pdist
  dist_cell_s <- (dist_cell$Pdist_Score)^lambda
  area_score <- area(gene_name, cell_index)^beta
  score <- corr_gene_s*dist_cell_s*area_score
  gene_index <- Getgene_index(gene_name)
  back <- list('Gene Names'= gene_name,'Cell Names'=CELL[cell_index],
               'Score'=score,'Gene Corr Score'= corr_gene_s,'Cell Dist Score'=dist_cell_s,
               'Gene Corr'= corr_gene_all,'Cell Dist'=dist_cell_all,'Area Score'=area_score)
  return(back)
}

score(c('pha-4', 'tbx-11', 'B0310.2'), grep('E', CELL),
      alpha=1, lambda = 1, beta = 1, method = 'single')






##画表达图
drawBiclusterExpressed <- function(bicluster = B){
  library(pheatmap)
  colorName <- c(rgb(146,197,222, maxColorValue = 255), rgb(217,83,79, maxColorValue = 255))
  gene_index <- B$`gene index`
  Cell_names <- B$`cell names`
  bicluster_expressed <- Data_3[sort(gene_index), which(colnames(Data_3) %in% Cell_names)]
  if(sum(bicluster_expressed) < length(gene_index) * length(Cell_names)){
    pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName, legend = F)
  }
  else{
    pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1)
  }
}

##画新表达图，合并同基因的不同文件，表达按check_express中的，即百分之60以上的基因文件有表达才算有表达


drawBiclusterExpressed <- function(bicluster = B){
  library(pheatmap)
  colorName <- c(rgb(146,197,222, maxColorValue = 255), rgb(217,83,79, maxColorValue = 255))
  gene_name <- B$`gene names`
  gene_name <- sort(gene_name)
  Cell_names <- B$`cell names`
  Cell_names <- sort(Cell_names)
  bicluster_expressed <- c()
  cindex <- c()
  for(j in 1:length(Cell_names)){
    cindex <- c(cindex, which(CELL %in% Cell_names[j]))
  }
  for(i in 1:length(gene_name)){
    gene_index <- Getgene_index(gene_name[i])
    bicluster_row <- Data_3[sort(gene_index), cindex]
    r <- colMeans(bicluster_row)
    for(j in 1:length(Cell_names)){
      if(r[j] >= 3/5){
        r[j] <- 1
      }
      else{
        r[j] <- 0
      }
    }
    bicluster_expressed <- rbind(bicluster_expressed, r)
  }
  rownames(bicluster_expressed) <- gene_name
  colnames(bicluster_expressed) <- Cell_names
  if(sum(bicluster_expressed) < length(gene_name) * length(Cell_names)){
    pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName, legend = F)
  }
  else{
    pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1)
  }
}

B <- All_B[[2]]
drawBiclusterExpressed(B)

###


palette.colors(n = 7)

B <- GABiclustering_initialize(10, 20, 30)[[1]]
B <- B[[1]]
B <- All_B[[2]]
drawBiclusterExpressed(B)
B <- B[[1]]
B <- GAbiclustering_all_result[[1]]

pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName,
         filename = '1.pdf')

pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1)

pheatmap(Data_3, show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F, color = colorName)

bicluster_expressed <- Data_3[, grep('E', CELL)]
bicluster_expressed <- bicluster_expressed[which(rowMeans(bicluster_expressed) != 0), ]
pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName)
setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")


pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName,
         filename = '1.pdf')

pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName,
         filename = '1.png')

