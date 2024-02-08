
##All代表所有找到的双聚类的基因，细胞指标
B <- GAbiclustering_all_result[[1]]
gene_index <- B[[1]]
gene_name <- B[[2]]
cell_index <- B[[3]]
All_gene_index <- gene_index
All_cell_index <- cell_index


##单独写惩罚项函数


overlapping <- function(gene_name, cell_index, delta = 0.1, 
                        All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  gene_index <- Getgene_index(gene_name)
  g_o <- intersect(gene_index,All_gene_index)
  c_o <- intersect(cell_index,All_cell_index)
  o <- exp(-delta * (length(g_o) * length(c_o)))
  back <- list('penalty' = o, 'overlapping_gene_index' = g_o, 'overlapping_cell_index' = c_o)
  return(back)
  
}

overlapping(gene_name, 1:10, delta = 0.1, All_gene_index = All_gene_index, All_cell_index = All_cell_index)

## 目标函数加惩罚项(找全部双聚类的目标函数)
##如果相关系数为负数，则算相关系数的次方值是NaN，无所谓，反正这个双聚类也不好
score <- function(gene_name, cell_index, alpha=1, lambda = 1, beta = 1, delta = 0.1,
                  method = 'average', All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  corr_gene <- gene_corr(gene_name, cell_index, method)
  corr_gene_all <- corr_gene$'Gcorr'
  corr_gene_s <- (corr_gene$'Gcorr_Score')^alpha
  dist_cell <- cell_pdist(gene_name, cell_index)
  dist_cell_all <- dist_cell$Pdist
  dist_cell_s <- (dist_cell$Pdist_Score)^lambda
  area_score <- area(gene_name, cell_index)^beta
  o <- overlapping(gene_name, cell_index, delta, All_gene_index, All_cell_index)$penalty
  score <- corr_gene_s * dist_cell_s * area_score * o
  back <- list('Gene Names'= gene_name,'Cell Names'= CELL[cell_index],
               'Score'=score,'Gene Corr Score'= corr_gene_s,'Cell Dist Score'=dist_cell_s,
               'Gene Corr'= corr_gene_all,'Cell Dist'=dist_cell_all,'Area Score'=area_score, 'penalty' = o)
  return(back)
}

score(gene_name,cell_index,alpha=1,lambda = 1.1,beta = 1.35, delta = 0.1, method = 'average',
      All_gene_index = All_gene_index, All_cell_index = All_cell_index)


##fitness函数也要改
GABiclustering_fitness <- function(B, alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                                   All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  fitness <- vector('list',length(B))
  s <- NULL
  c <- NULL
  d <- NULL
  a <- NULL
  o <- NULL
  for(i in 1:length(B)){
    li <- B[[i]]
    fitness[[i]] <- score(gene_name = li$`gene names`,cell_index = li$`cell index`, 
                          alpha = alpha, lambda = lambda, beta = beta, delta = delta, 
                          method = 'average', All_gene_index, All_cell_index)
    s <- append(s,fitness[[i]]$'Score')
    c <- append(c,fitness[[i]]$'Gene Corr Score')
    d <- append(d,fitness[[i]]$'Cell Dist Score')
    a <- append(a,fitness[[i]]$'Area Score')
    o <- append(o,fitness[[i]]$'penalty')
  }
  s_ordered <- sort(s,decreasing = T)
  s_ordered_index <- order(s,decreasing = T)
  c_ordered <- c[s_ordered_index]
  d_ordered <- d[s_ordered_index]
  a_ordered <- a[s_ordered_index]
  o_ordered <- o[s_ordered_index]
  back <- list('s_ordered' = round(s_ordered,5), 's_ordered_index' = s_ordered_index,
               'c_ordered' = round(c_ordered,5), 'd_ordered' = round(d_ordered,5),
               'a_ordered' = round(a_ordered,5), 'o_ordered' = round(o_ordered,5))
  return(back)
}

B <- GABiclustering_initialize()
GABiclustering_fitness(B, alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                       All_gene_index = All_gene_index, All_cell_index = All_cell_index)

###FindAllBiclusterings的GABiclustering函数
#max_iteration，最大迭代次数
#N，初始种群个数
#m，初始种群中的最大基因个数
#n，初始种群中的最大细胞个数
#alpha，细胞面积权重
#beta，基因相关系数权重
#lambda，细胞距离权重
#delta，惩罚项权重
#good_num，每次保留适应度最大值的个数，称为好的组，剩下中等组和坏的组大致各一半
#add_p_coef调整面积分数转为0-1的权重
#gene_result_address，基因结果保存地址
#cell_result_address，细胞结果保存地址
#k，初始定义的交叉互换基因或细胞的最大个数
#l，初始定义的一次添加基因或细胞的最大个数
#h，初始定义的一次删除基因或细胞的最大个数
#seed，随机种子
#continue，看是否从上次的结果开始
##和fate数据集不一样的，这里nice组取前5组，而不是前10组


GABiclustering <- function(max_iteration = 1000, N = 100, m = 4, n = 9, alpha = 1, 
                           lambda = 1, beta = 1, delta = 0.1, good_num = 20,
                           All_gene_index = All_gene_index, All_cell_index = All_cell_index,
                           Biclustering_num = Biclustering_num, 
                           continue = 0, seed = 2021){
  set.seed(seed)
  same_times <- 0
  mid_num <- ceiling((N-good_num)/2)
  save_names <- paste('B_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.Rdata',sep = '')
  if(continue == 0){
    B <- GABiclustering_initialize(N, m, n)
  }
  # else{
  #   load(save_names)
  # }
  output_names <- paste('output_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.txt',sep = '')
  fitness_old <- GABiclustering_fitness(B, alpha, lambda, beta, delta, All_gene_index, All_cell_index)
  sink(output_names)
  cat('Initial Score:',fitness_old$s_ordered[1:5],'\n')
  for(i in 1:max_iteration){
    fitness_old <- GABiclustering_fitness(B, alpha, lambda, beta, delta, All_gene_index, All_cell_index)
    best_index_old <- fitness_old$s_ordered_index[1]
    good_index_old <- fitness_old$s_ordered_index[1:good_num]
    nice_index_old <- fitness_old$s_ordered_index[1:5]
    mid_index <- fitness_old$s_ordered_index[(good_num+1):(good_num+mid_num)]
    bad_index <- fitness_old$s_ordered_index[(good_num+mid_num+1):N]
    good_b_old <- B[good_index_old]
    mid_b_old <- B[mid_index]
    bad_b_old <- B[bad_index]
    nice_b_old <- B[nice_index_old]
    good_m_new <- lapply(good_b_old,GABiclustering_mutation)
    mid_m_new <- lapply(mid_b_old,GABiclustering_mutation)
    bad_m_new <- lapply(bad_b_old,GABiclustering_mutation)
    good_c_new <- lapply(good_b_old,GABiclustering_crossover,goodbicluster = good_b_old)
    mid_c_new <- lapply(mid_b_old,GABiclustering_crossover,goodbicluster = good_b_old)
    bad_c_new <- lapply(bad_b_old,GABiclustering_crossover,goodbicluster = good_b_old)
    B_new <- c(B,good_m_new,mid_m_new,bad_m_new,good_c_new,mid_c_new,bad_c_new)
    B_new <- unique(B_new)
    fitness_new <- GABiclustering_fitness(B_new, alpha, lambda, beta, delta, All_gene_index, All_cell_index)
    new_index <- fitness_new$s_ordered_index[1:N]
    B <- B_new[new_index]
    fitness_new <- GABiclustering_fitness(B, alpha, lambda, beta, delta, All_gene_index, All_cell_index)
    nice_index_new <- fitness_new$s_ordered_index[1:5]
    nice_b_new <- B[nice_index_new]
    if(identical(nice_b_old, nice_b_new)){
      same_times <- same_times+1
    }
    else{ same_times <- 0 }
    if(same_times == 1){
      print('Local maximum reached.')
      sink()
      return(B)
    }
    best_index_new <- fitness_new$s_ordered_index[1]
    if(i %% 1 == 0){
      best_index <- fitness_new$s_ordered_index[1]
      cat('iteration:',i,'\n')
      cat('Score:',fitness_new$s_ordered[1:5],'\n')  
      cat('Gene Corr Score:',fitness_new$c_ordered[1:5],'\n')
      cat('Cell Dist Score:',fitness_new$d_ordered[1:5],'\n')
      cat('Area Score:',fitness_new$a_ordered[1:5],'\n')
      cat('Penalty:',fitness_new$o_ordered[1:5],'\n')
      cat('best gene index:',B[[best_index]]$`gene index`,'\n')
      cat('best gene names:',B[[best_index]]$`gene names`,'\n')
      cat('best cell index:',B[[best_index]]$`cell index`,'\n')
      cat('best cell names:',B[[best_index]]$`cell names`,'\n')
    }
    save(B, file = save_names)
  }
  print('Maximum number of iterations reached.')
  return(B)
  sink()
}

##Biclustering_num 从第几个双聚类开始找
##start = 1 算法开始条件
##注意设置alpha，beta，lambda，delta，seed，B，All_gene_index, All_cell_index
##注意设置continue，Biclustering_num，seed参数，看是否有创建新的All_B
B <- GAbiclustering_all_result[[1]]
B <- All_B[[1]]
B <- B[[1]]
All_gene_index <- c()
All_cell_index <- c()
All_gene_index <- B$`gene index`
All_cell_index <- B$`cell index`
alpha <- 0.35
lambda <- 0.2
beta <- 0.4
delta <- 0.05
seed <- 737

FindAllBiclusterings <- function(Biclustering_num = 2, start = 1, delta = delta,
                                 All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  # All_B <- vector('list', 1)
  # All_B[[1]] <- B
  # 
  while(start == 1){
    seed <- sample(1:1000,1)
    GAbiclustering_all_result <- GABiclustering(max_iteration = 1000, N = 100, m = 10, n = 15, alpha = alpha, 
                                                lambda = lambda, beta = beta, delta = delta, good_num = 20,
                                                All_gene_index = All_gene_index, All_cell_index = All_cell_index,
                                                Biclustering_num = Biclustering_num, 
                                                continue = 0, seed = seed)

    B <- GAbiclustering_all_result[[1]]
    gene_index <- B[[1]]
    cell_index <- B[[3]]
    Cell_names <- B$`cell names`
    if(Biclustering_num <= 30){
      All_gene_index <- c(All_gene_index, gene_index)
      All_gene_index <- unique(All_gene_index)
      All_cell_index <- c(All_cell_index, cell_index)
      All_cell_index <- unique(All_cell_index)
      All_B[[Biclustering_num]] <- B
      All_gene_index <- sort(All_gene_index)
      All_cell_index <- sort(All_cell_index)
      save(All_gene_index,file=paste('All_gene_index_',alpha,'_',lambda,'_',beta,'_',seed,
                                     '_',delta,'_',Biclustering_num,'.Rdata',sep = ''))
      save(All_cell_index,file=paste('All_cell_index_',alpha,'_',lambda,'_',beta,'_',seed,
                                     '_',delta,'_',Biclustering_num,'.Rdata',sep = ''))
      save(All_B,file=paste('All_B_',alpha,'_',lambda,'_',beta,'_',seed,
                            '_',delta,'.Rdata',sep = ''))
      g_length <- length(gene_index)
      c_length <- length(Cell_names)
      bicluster_expressed <- Data_3[sort(gene_index), which(colnames(Data_3) %in% Cell_names)]
      filename <- paste('express_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.png',sep = '')
      if(sum(bicluster_expressed) !=  g_length * c_length){
        pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName, legend = F,
                 filename = filename)
      }
      else{
        pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1,
                 filename = filename)
      }
      
      Biclustering_num <- Biclustering_num + 1
    }
    else{
      start = 0
    }
  }
}

##continue为1，单独写

FindAllBiclusterings <- function(Biclustering_num = 2, start = 1, delta = delta,
                                 All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  All_B <- vector('list', 1)
  All_B[[1]] <- B
  
  while(start == 1){
    if(Biclustering_num == 1){
      GAbiclustering_all_result <- GABiclustering(max_iteration = 1000, N = 100, m = 10, n = 15, alpha = alpha, 
                                                  lambda = lambda, beta = beta, delta = delta, good_num = 20,
                                                  All_gene_index = All_gene_index, All_cell_index = All_cell_index,
                                                  Biclustering_num = Biclustering_num, 
                                                  continue = 1, seed = seed)
      
    }
    else{
      seed <- sample(1:1000,1)
      GAbiclustering_all_result <- GABiclustering(max_iteration = 1000, N = 100, m = 10, n = 15, alpha = alpha, 
                                                  lambda = lambda, beta = beta, delta = delta, good_num = 20,
                                                  All_gene_index = All_gene_index, All_cell_index = All_cell_index,
                                                  Biclustering_num = Biclustering_num, 
                                                  continue = 0, seed = seed)
    }
    B <- GAbiclustering_all_result[[1]]
    gene_index <- B[[1]]
    cell_index <- B[[3]]
    Cell_names <- B$`cell names`
    if(Biclustering_num <= 30){
      All_gene_index <- c(All_gene_index, gene_index)
      All_gene_index <- unique(All_gene_index)
      All_cell_index <- c(All_cell_index, cell_index)
      All_cell_index <- unique(All_cell_index)
      All_B[[Biclustering_num]] <- B
      All_gene_index <- sort(All_gene_index)
      All_cell_index <- sort(All_cell_index)
      save(All_gene_index,file=paste('All_gene_index_',alpha,'_',lambda,'_',beta,'_',seed,
                                     '_',delta,'_',Biclustering_num,'.Rdata',sep = ''))
      save(All_cell_index,file=paste('All_cell_index_',alpha,'_',lambda,'_',beta,'_',seed,
                                     '_',delta,'_',Biclustering_num,'.Rdata',sep = ''))
      save(All_B,file=paste('All_B_',alpha,'_',lambda,'_',beta,'_',seed,
                            '_',delta,'.Rdata',sep = ''))
      g_length <- length(gene_index)
      c_length <- length(Cell_names)
      bicluster_expressed <- Data_3[sort(gene_index), which(colnames(Data_3) %in% Cell_names)]
      filename <- paste('express_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.pdf',sep = '')
      if(sum(bicluster_expressed) !=  g_length * c_length){
        pheatmap(bicluster_expressed, cluster_rows = F, cluster_cols = F, color = colorName,
                 filename = filename)
      }
      else{
        pheatmap(bicluster_expressed, legend = F, cluster_rows = F, cluster_cols = F, color = colorName[2], breaks = 1,
                 filename = filename)
      }
      
      Biclustering_num <- Biclustering_num + 1
    }
    else{
      start = 0
    }
  }
}



ALL_B <- FindAllBiclusterings(Biclustering_num = 2, start = 1, delta = delta, All_gene_index = All_gene_index, All_cell_index = All_cell_index)


###输出各个双聚类的基因和细胞名称，用txt文档
setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")
output_names <- 'All_B.txt'
sink(output_names)
for(i in 1:30){
  cat(i,'\n')
  gene_index <- All_B[[i]]$`gene index`
  gene_name <- All_B[[i]]$`gene names`
  cell_name <- All_B[[i]]$`cell names`
  output <- paste(sort(gene_name), collapse = ',')
  cat('best gene names:',output,'\n')
  output <- paste(sort(cell_name), collapse = ',')
  cat('best cell names:',output,'\n')
  g_length <- length(gene_index)
  c_length <- length(cell_name)
  bicluster_expressed <- Data_3[sort(gene_index), which(colnames(Data_3) %in% cell_name)]
  num <- sum(bicluster_expressed) / (g_length * c_length) * 100
  cat('gene expressed area:',num,'%', '\n')

}
sink()


s <- c()
for(i in 1:length(All_gene_index)){
  ch <- GENE_names[i]
  x <- strsplit(ch, '_')[[1]][1]
  s <- c(s, x)
}
s <- unique(s)
s


which(All_cell_index %in% cell_missing_index)

All_cell_index[447]
