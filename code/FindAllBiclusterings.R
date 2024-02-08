
##All代表所有找到的双聚类的基因，细胞指标

All_gene_index <- gene_index
All_cell_index <- cell_index

##单独写惩罚项函数


overlapping <- function(gene_index = 1:10, cell_index = 1:10, delta = 0.1, 
                        All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  
  g_o <- intersect(gene_index,All_gene_index)
  c_o <- intersect(cell_index,All_cell_index)
  o <- exp(-delta * (length(g_o) * length(c_o)))
  back <- list('penalty' = o, 'overlapping_gene_index' = g_o, 'overlapping_cell_index' = c_o)
  return(back)

}

overlapping(1:10, 2:20, delta = 0.1)


## 目标函数加惩罚项(找全部双聚类的目标函数)

score <- function(gene_index,cell_index,alpha=1,lambda = 1,beta = 1, delta = 0.1,
                  All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  corr_gene <- (gene_corr_poly(gene_index, cell_index)$Gcorr_Score)^alpha
  dist_cell <- (cell_pdist(gene_index,cell_index)$Pdist_Score)^lambda
  area_score <- area(gene_index, cell_index)^beta
  o <- overlapping(gene_index,cell_index, delta, All_gene_index, All_cell_index)$penalty
  score <- corr_gene * dist_cell * area_score * o
  back <- list('Gene Names'=GENE_names[gene_index],'Cell Names'=CELL[cell_index],
               'Score'=score,'Gene Corr Score'=corr_gene,'Cell Dist Score'=dist_cell,
               'Area Score'=area_score, 'penalty' = o)
  
  return(back)
}

score(gene_index,cell_index,alpha=1.4,lambda = 1.1,beta = 1, delta = 0,
      All_gene_index = All_gene_index, All_cell_index = All_cell_index)
score(c(34,5,35,41,7,1,2),c(78,52,68,70:77))

score(sample(51,5),sample(145,15),alpha=1,lambda = 4.2,beta = 1.04, delta = 0.1,
      All_gene_index = All_gene_index, All_cell_index = All_cell_index)


#计算适应度，并排序，只按得分排序
GABiclustering_fitness <- function(B, alpha = 0.004, lambda = 1,beta = 1, delta = 0.1,
                                   All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  fitness <- vector('list',length(B))
  s <- NULL
  c <- NULL
  d <- NULL
  a <- NULL
  o <- NULL
  for(i in 1:length(B)){
    fitness[[i]] <- score(gene_index = B[[i]]$`gene index`,cell_index = B[[i]]$`cell index`, 
                          alpha = alpha,lambda = lambda,beta = beta,delta = delta, All_gene_index, All_cell_index)
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





##找全部双聚类的函数


delta = 0.1


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
#result_num，迭代多少次保存一次结果

GABiclustering <- function(max_iteration = 1000, N = 100, m = 4, n = 9, alpha = 0.8, 
                           lambda = 1, beta = 1, delta = 0.1, good_num = 20,
                           gene_result_address = 'GABiclustering_gene_result_',
                           cell_result_address = 'GABiclustering_cell_result_',
                           All_gene_index = All_gene_index, All_cell_index = All_cell_index,
                           Biclustering_num = Biclustering_num, 
                           result_num = 20, seed = 2021){
  set.seed(seed)
  same_times <- 0
  mid_num <- ceiling((N-good_num)/2)
  B <- GABiclustering_initialize(N,m,n)
  biclustering_gene_result_names <- paste(gene_result_address,alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.csv',sep = '')
  biclustering_cell_result_names <- paste(cell_result_address,alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.csv',sep = '')
  output_names <- paste('output_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.txt',sep = '')
  fitness_old <- GABiclustering_fitness(B, alpha, lambda, beta, delta, All_gene_index, All_cell_index)
  sink(output_names)
  cat('Initial Score:',fitness_old$s_ordered[1:10],'\n')
  for(i in 1:max_iteration){
    fitness_old <- GABiclustering_fitness(B, alpha, lambda, beta, delta, All_gene_index, All_cell_index)
    best_index_old <- fitness_old$s_ordered_index[1]
    good_index_old <- fitness_old$s_ordered_index[1:good_num]
    nice_index_old <- fitness_old$s_ordered_index[1:10]
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
    nice_index_new <- fitness_new$s_ordered_index[1:10]
    nice_b_new <- B[nice_index_new]
    if(identical(nice_b_old, nice_b_new)){
      same_times <- same_times+1
    }
    else{ same_times <- 0 }
    if(same_times == 5){
      print('Local maximum reached.')
      sink()
      return(B)
    }
    best_index_new <- fitness_new$s_ordered_index[1]
    if(i %% result_num == 0){
      best_index <- fitness_new$s_ordered_index[1]
      cat('iteration:',i,'\n')
      cat('Score:',fitness_new$s_ordered[1:10],'\n')  
      cat('Gene Corr Score:',fitness_new$c_ordered[1:10],'\n')
      cat('Cell Dist Score:',fitness_new$d_ordered[1:10],'\n')
      cat('Area Score:',fitness_new$a_ordered[1:10],'\n')
      cat('Penalty:',fitness_new$o_ordered[1:10],'\n')
      cat('best gene index:',B[[best_index]]$`gene index`,'\n')
      cat('best gene names:',B[[best_index]]$`gene names`,'\n')
      cat('best cell index:',B[[best_index]]$`cell index`,'\n')
      cat('best cell names:',B[[best_index]]$`cell names`,'\n')
      biclustering_gene_result <- data.frame(B[[best_index]]$`gene index`,B[[best_index]]$`gene names`)
      biclustering_cell_result <- data.frame(B[[best_index]]$`cell index`,B[[best_index]]$`cell names`)
      colnames(biclustering_gene_result) <- c('gene index','gene name')
      colnames(biclustering_cell_result) <- c('cell index','cell name')
      write.csv(biclustering_gene_result,biclustering_gene_result_names)
      write.csv(biclustering_cell_result,biclustering_cell_result_names)
    }
  }
  print('Maximum number of iterations reached.')
  return(B)
  sink()
}

##bic_num = 100
##Biclustering_num 从第几个双聚类开始找
##start = 1 算法开始条件
##注意设置alpha，beta，lambda，delta，seed，B，All_gene_index, All_cell_index
##

All_gene_index <- gene_index
All_cell_index <- cell_index


FindAllBiclusterings <- function(bic_num = 100, Biclustering_num = 2, start = 1, delta = delta,
                                 All_gene_index = All_gene_index, All_cell_index = All_cell_index){
  All_B <- vector('list', 1)
  All_B[[1]] <- B
  
  while(start == 1){
    seed <- sample(500,1)
    GAbiclustering_fate_result <- GABiclustering(max_iteration = 1000, N = 100, m = 10, n = 15, alpha = alpha, 
                                                 lambda = lambda, beta = beta, delta = delta, good_num = 20,
                                                 gene_result_address = 'GABiclustering_fate_gene_result_',
                                                 cell_result_address = 'GABiclustering_fate_cell_result_',
                                                 All_gene_index = All_gene_index, All_cell_index = All_cell_index,
                                                 Biclustering_num = Biclustering_num, 
                                                 result_num = 1, seed = seed)
    
    save(GAbiclustering_fate_result,file=paste('GAbiclustering_fate_result_',alpha,'_',lambda,'_',beta,'_',seed,
                                               '_',delta,'_',Biclustering_num,'.Rdata',sep = ''))
    B <- GAbiclustering_fate_result[[1]]
    gene_index <- B[[1]]
    cell_index <- B[[3]]
    if(Biclustering_num <= 10){
      pdf(paste('express_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.pdf',sep = ''),width = 12,height = 9)
      drawBiclusterExpressed(B)
      dev.off()
      pdf(paste('cell_',alpha,'_',lambda,'_',beta,'_',seed,'_',delta,'_',Biclustering_num,'.pdf',sep = ''),width = 12,height = 7)
      CatCellFate(cell_index)
      dev.off()
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
      Biclustering_num <- Biclustering_num + 1
    }
    else{
      start = 0
    }

  }
  All_gene_index <- sort(All_gene_index)
  All_cell_index <- sort(All_cell_index)
  back <- list('All_Biclusterings' = All_B, 'All_gene_index' = All_gene_index, 'All_gene_names' = GENE_names[All_gene_index], 
               'All_cell_index' = All_cell_index, 'All_cell_names' = CELL[All_cell_index])
  return(back)
  
}




