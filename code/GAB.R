##计算有效面积

gene_index <- B[[1]]
B[[2]]
B[[4]]
cell_index <- B[[3]]

area <- function(gene_index, cell_index){

  bicluster_expressed <- Data_3_simulation[gene_index,cell_index]
  area_ <- sum(bicluster_expressed)
  return(area_)

}

area(gene_index, cell_index)

###用p值表示细胞间距离，创建p值距离分函数(p值越大距离越近)
cell_pdist <- function(gene_index, cell_index){
  ptotal <- NULL
  for(j in 1:(length(cell_index)-1)){
    c_ind_1 <- cell_index[j]
    for(k in (j+1):length(cell_index)){
      c_ind_2 <- cell_index[k]
      pvalue <- NULL
      for(i in 1:length(gene_index)){
        ind <- gene_index[i]
        flag_1 <- which(CELL_fate_diff[[ind]][,1] == CELL[c_ind_1])
        flag_2 <- which(CELL_fate_diff[[ind]][,1] == CELL[c_ind_2])
        if(length(flag_1) >=1 & length(flag_2) >=1){
          d_1 <- CELL_fate_diff[[ind]][flag_1,3]
          d_2 <- CELL_fate_diff[[ind]][flag_2,3]
          p_ <- ks.test(d_1, d_2, alternative = 'two.sided', exact = F)$p.value
          pvalue <- c(pvalue, p_)
        }
      }
      if(!is.null(pvalue)){
        pvalue <- min(pvalue)
        ptotal <- c(ptotal, pvalue)
      }
    }
  }
  pdist_score <- mean(ptotal)
  back <- list('Gene Names'=GENE_names[gene_index],'Cell Names'=CELL[cell_index],
               'Pdist' = ptotal,'Pdist_Score'=pdist_score)
  return(back)
}

###最新版目标函数（用函数相关系数，理论值，而不是取点）(实际上这个函数没有变，只是更新了gene_corr_poly)


score <- function(gene_index,cell_index,alpha=1,lambda = 1,beta = 1){
  corr_gene <- gene_corr_poly(gene_index,cell_index)
  corr_gene_all <- corr_gene$'Gcorr'
  corr_gene_s <- (corr_gene$'Gcorr_Score')^alpha
  dist_cell <- cell_pdist(gene_index,cell_index)
  dist_cell_all <- dist_cell$Pdist
  dist_cell_s <- (dist_cell$Pdist_Score)^lambda
  area_score <- area(gene_index, cell_index)^beta
  score <- corr_gene_s*dist_cell_s*area_score
  back <- list('Gene Names'=GENE_names[gene_index],'Cell Names'=CELL[cell_index],
               'Score'=score,'Gene Corr Score'= corr_gene_s,'Cell Dist Score'=dist_cell_s,
               'Gene Corr'= corr_gene_all,'Cell Dist'=dist_cell_all,'Area Score'=area_score)
  return(back)
}

###
###
###
###新增一个检查函数，检查当前双聚类是否满足一定高的表达比例
###输入和遗传算法中一致，是一个四维列表，分别为基因指标，基因名称，细胞指标，细胞名称
###参数pro表示要求各细胞在基因的表达比例最低阈值，各基因在细胞的表达比例最低阈值
##基因和细胞个数都超过3个才执行，如果删除导致基因或细胞个数小于2，直接重新随机初始化同维度双聚类
##删除顺序：不分细胞和基因，每次删除不表达比例最高的，迭代删除


check_express <- function(bicluster, pro = 3/5){
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index`
  l_g <- length(gene_index)
  l_c <- length(cell_index)
  if(l_g < 3 && l_c < 3){
    return(bicluster)
  }
  dat <- Data_3[gene_index, cell_index, drop = F]
  c_expr <- colMeans(dat)
  g_expr <- rowMeans(dat)
  expr <- c(c_expr, g_expr)
  min_expr <- min(expr)
  if(min_expr >= pro){
    return(bicluster)
  }
  else{
    ind <- which(expr %in% min_expr)[1]
    if(ind <= l_c){
      cell_index <- cell_index[-ind]
    }
    else{
      ind <- ind - l_c
      gene_index <- gene_index[-ind]
    }
  }
  if(length(gene_index) < 2 || length(cell_index) < 2){
    gene_index <- sample(1:length(GENE_names), l_g)
    dat <- Data_3[gene_index,]
    c_expressed <- as.numeric(which(colMeans(dat) != 0))
    if(length(c_expressed) >= l_c){
      cell_index <- sample(c_expressed, l_c)
    }
    else{
      cell_index <- sample(1:length(CELL), l_c)
    }
    bicluster <- list('gene index' = gene_index, 'gene names' = GENE_names[gene_index],
                      'cell index' = cell_index, 'cell names' = CELL[cell_index])
    return(bicluster)
  }
  else{
    bicluster <- list('gene index' = gene_index, 'gene names' = GENE_names[gene_index],
                      'cell index' = cell_index, 'cell names' = CELL[cell_index])
    bicluster <- check_express(bicluster, pro)
    return(bicluster)
  }
}


B <- GABiclustering_initialize(1, 10, 40)[[1]]
drawBiclusterExpressed(B)
B <- check_express(B)
drawBiclusterExpressed(B)

#初始化种群，N为初始种群个数，
#m为初始种群中的最大基因个数
#n为初始种群中的最大细胞个数
#双聚类中不含无表达的细胞
GABiclustering_initialize <- function(N = 10, m = 8, n = 12){
  
  B <- vector('list',N)
  m <- c(2:m)  
  n <- c(2:n)   
  for(i in 1:N){
    g <- sample(1:length(GENE_names),sample(m,1))    
    #sample(X,size,replace=T/F)：从数据X中随机抽取size个样本， 当replace=T时为放回取样，为F时为不放回抽样(默认无放回F)。
    dat <- Data_3[g,]
    c_expressed <- as.numeric(which(colMeans(dat) != 0))   ##不等于0
    c_unexpressed <- as.numeric(which(colMeans(dat) == 0))   ##加
    num <- sample(n,1)
    if(length(c_expressed) >= num){
      c <- sample(c_expressed, num)
    }
    else{
      c <- c(c_expressed, sample(c_unexpressed, num-length(c_expressed)))
    }
    B[[i]][['gene index']] <- g
    B[[i]][['gene names']] <- GENE_names[g]
    B[[i]][['cell index']] <- c
    B[[i]][['cell names']] <- CELL[c]
  }
  
  return(B)
}

li <- GABiclustering_initialize()
###
###

#计算适应度，并排序，只按得分排序
GABiclustering_fitness <- function(B, alpha = 0.004, lambda = 1,beta = 1){
  fitness <- vector('list',length(B))
  s <- NULL
  c <- NULL
  d <- NULL
  a <- NULL
  for(i in 1:length(B)){
    fitness[[i]] <- score(gene_index = B[[i]]$`gene index`,cell_index = B[[i]]$`cell index`, alpha = alpha,lambda = lambda,beta = beta)
    s <- append(s,fitness[[i]]$'Score')
    c <- append(c,fitness[[i]]$'Gene Corr Score')
    d <- append(d,fitness[[i]]$'Cell Dist Score')
    a <- append(a,fitness[[i]]$'Area Score')
  }
  s_ordered <- sort(s,decreasing = T)
  s_ordered_index <- order(s,decreasing = T)
  c_ordered <- c[s_ordered_index]
  d_ordered <- d[s_ordered_index]
  a_ordered <- a[s_ordered_index]
  back <- list('s_ordered' = round(s_ordered,5), 's_ordered_index' = s_ordered_index,
               'c_ordered' = round(c_ordered,5), 'd_ordered' = round(d_ordered,5),
               'a_ordered' = round(a_ordered,5))
  return(back)
}
###
#随机添加基因，优先添加表达细胞多的
GABiclustering_gene_add <- function(bicluster = B[[1]], l = 1){
  
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index`
  gene_add_index <- c(1:length(GENE_names))[-which(c(1:length(GENE_names)) %in% gene_index)]
  dat <- Data_3[,cell_index]
  g_expressed <- as.numeric(which(rowMeans(dat) >= 1/2))
  gene_add_index <- intersect(gene_add_index, g_expressed)
  if(length(gene_add_index) == 0){
    return(bicluster)
  }
  else if(length(gene_add_index) >= 1 & length(gene_add_index) <= l){
    gene_index <- c(bicluster[[1]],gene_add_index)
    gene_names <- c(bicluster[[2]],GENE_names[gene_add_index])  
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_names
    return(bicluster)
  }
  else{
    gene_add_index <- sample(gene_add_index,l)
    gene_index <- c(bicluster[[1]],gene_add_index)
    gene_names <- c(bicluster[[2]],GENE_names[gene_add_index])  
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_names
    return(bicluster)
  }
  
}

###
###


#随机添加细胞，优先添加有超过一半表达的
GABiclustering_cell_add <- function(bicluster = B[[1]], h = 1){
  
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index` 
  cell_add_index <- c(1:length(CELL))[-which(c(1:length(CELL)) %in% cell_index)]
  dat <- Data_3[gene_index,]
  c_expressed <- as.numeric(which(colMeans(dat) >= 1/2))
  cell_add_index <- intersect(cell_add_index,c_expressed)
  if(length(cell_add_index) == 0){
    return(bicluster)
  }
  else if(length(cell_add_index) >= 1 & length(cell_add_index) <= h){
    cell_index <- c(bicluster[[3]],cell_add_index)
    cell_names <- c(bicluster[[4]],CELL[cell_add_index])  
    bicluster$`cell index` <- cell_index
    bicluster$`cell names` <- cell_names
    return(bicluster)
  }
  else{
    cell_add_index <- sample(cell_add_index,h)
    cell_index <- c(bicluster[[3]],cell_add_index)
    cell_names <- c(bicluster[[4]],CELL[cell_add_index])  
    bicluster$`cell index` <- cell_index
    bicluster$`cell names` <- cell_names
    return(bicluster)
    
  }
}

#随机删除基因。优先删除表达少的

GABiclustering_gene_remove <- function(bicluster = B[[1]], l = 1){
  
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index`  
  dat <- Data_3[gene_index,cell_index]  
  g_notexpressed <- as.numeric(which(rowMeans(dat) < 1/2))
  if(length(gene_index) <= l+2){
    return(bicluster)
  }
  else if(length(g_notexpressed) == 0){
    gene_remove_index <- sample(gene_index,l)
    gene_index <- gene_index[-which(gene_index %in% gene_remove_index)]
    gene_names <- GENE_names[gene_index]  
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_names
    return(bicluster)
  }
  else if(length(g_notexpressed) <= l & length(g_notexpressed) >= 1){
    gene_remove_index <- which(GENE_names %in% rownames(dat)[g_notexpressed])
    gene_index <- gene_index[-which(gene_index %in% gene_remove_index)]
    gene_names <- GENE_names[gene_index]  
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_names
    return(bicluster)
  }
  else{
    gene_remove_index <- which(GENE_names %in% rownames(dat)[g_notexpressed])
    gene_remove_index <- sample(gene_remove_index,l)
    gene_index <- gene_index[-which(gene_index %in% gene_remove_index)]
    gene_names <- GENE_names[gene_index] 
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_names
    return(bicluster)
  }
  
}

###
####随机删除细胞,优先删除表达少的

GABiclustering_cell_remove <- function(bicluster = B[[1]], h = 1){
  
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index`
  dat <- Data_3[gene_index,cell_index]
  c_notexpressed <- as.numeric(which(colMeans(dat) < 1/2))
  if(length(cell_index) <= h+2){
    return(bicluster)
  }
  else if(length(c_notexpressed) == 0){
    cell_remove_index <- sample(cell_index,h)
    cell_index <- cell_index[-which(cell_index %in% cell_remove_index)]
    cell_names <- CELL[cell_index]  
    bicluster$`cell index` <- cell_index
    bicluster$`cell names` <- cell_names
    return(bicluster)
  }
  else if(length(c_notexpressed) <= h & length(c_notexpressed) >= 1){
    cell_remove_index <- which(CELL %in% colnames(dat)[c_notexpressed])
    cell_index <- cell_index[-which(cell_index %in% cell_remove_index)]
    cell_names <- CELL[cell_index]  
    bicluster$`cell index` <- cell_index
    bicluster$`cell names` <- cell_names
    return(bicluster)
  }
  else{
    cell_remove_index <- which(CELL %in% colnames(dat)[c_notexpressed])
    cell_remove_index <- sample(cell_remove_index,h)
    cell_index <- cell_index[-which(cell_index %in% cell_remove_index)]
    cell_names <- CELL[cell_index]  
    bicluster$`cell index` <- cell_index
    bicluster$`cell names` <- cell_names
    return(bicluster)
  }
  
}

#交叉互换，将部分好的基因随机给其他组
GABiclustering_gene_exchange <- function(badbicluster = B[[2]], goodbicluster = B){
  
  samp <- sample(length(goodbicluster),1)
  goodbicluster <- goodbicluster[[samp]]
  good_gene_index <- goodbicluster$`gene index`
  bad_gene_index <- badbicluster$`gene index`
  
  k <- 0
  while (k == 0) {
    k <- min(rpois(1, lambda = 1), min(nrow(good_gene_index), nrow(bad_gene_index)))
  }
  
  if(length(bad_gene_index) <= 3*k | length(good_gene_index) <= 3*k){
    badbicluster <- GABiclustering_mutation(badbicluster)
    return(badbicluster)
  }
  else{
    good_gene_exchange_index <- sample(1:length(good_gene_index),k)
    add_value <- good_gene_index[good_gene_exchange_index]
    if(length(bad_gene_index) >= max(good_gene_exchange_index)){
      bad_gene_index <- c(bad_gene_index[-good_gene_exchange_index],add_value)
      bad_gene_index <- unique(bad_gene_index)
      gene_names <- GENE_names[bad_gene_index]  
      badbicluster$`gene index` <- bad_gene_index
      badbicluster$`gene names` <- gene_names
      return(badbicluster)
    }
    else{
      remove_index <- sample(1:length(bad_gene_index),k)
      bad_gene_index <- c(bad_gene_index[-remove_index],add_value)
      bad_gene_index <- unique(bad_gene_index)
      gene_names <- GENE_names[bad_gene_index]  
      badbicluster$`gene index` <- bad_gene_index
      badbicluster$`gene names` <- gene_names
      return(badbicluster)
    }
  }
}

b <- GABiclustering_gene_exchange()


#交叉互换，将部分好的细胞随机给其他组

GABiclustering_cell_exchange <- function(badbicluster = B[[2]], goodbicluster = B){
  
  samp <- sample(length(goodbicluster),1)
  goodbicluster <- goodbicluster[[samp]]
  good_cell_index <- goodbicluster$`cell index`
  bad_cell_index <- badbicluster$`cell index`
  k <- 0
  while (k == 0) {
    k <- min(rpois(1, lambda = 1), min(nrow(good_cell_index), nrow(bad_cell_index)))
  }
  if(length(bad_cell_index) <= 3*k | length(good_cell_index) <= 3*k){
    badbicluster <- GABiclustering_mutation(badbicluster)
    return(badbicluster)
  }
  else{
    good_cell_exchange_index <- sample(1:length(good_cell_index),k)
    add_value <- good_cell_index[good_cell_exchange_index]
    if(length(bad_cell_index) >= max(good_cell_exchange_index)){
      bad_cell_index <- c(bad_cell_index[-good_cell_exchange_index],add_value)
      bad_cell_index <- unique(bad_cell_index)
      cell_names <- CELL[bad_cell_index]  
      badbicluster$`cell index` <- bad_cell_index
      badbicluster$`cell names` <- cell_names
      return(badbicluster)
    }
    else{
      remove_index <- sample(1:length(bad_cell_index),k)
      bad_cell_index <- c(bad_cell_index[-remove_index],add_value)
      bad_cell_index <- unique(bad_cell_index)
      cell_names <- CELL[bad_cell_index]  
      badbicluster$`cell index` <- bad_cell_index
      badbicluster$`cell names` <- cell_names
      return(badbicluster)
    }
  }
}

#交叉互换，将部分好的基因或细胞随机给其他组
GABiclustering_crossover <- function(badbicluster = B[[2]], goodbicluster = B){
  
  if(runif(1)>=0.5){
    badbicluster <- GABiclustering_cell_exchange(badbicluster, goodbicluster)
  }
  else{
    badbicluster <- GABiclustering_gene_exchange(badbicluster, goodbicluster)
  }
  badbicluster <- check_express(badbicluster)
  return(badbicluster)
  
}

# 定义截断泊松分布的概率质量函数
truncated_poisson_pmf <- function(x) {
  prob_zero <- 0.1
  prob_other <- 0.9
  
  if (x == 0) {
    return(prob_zero)
  } else {
    peak_lambda <- 1
    prob_x <- (exp(-peak_lambda) * peak_lambda^x) / factorial(x)
    return(prob_other * prob_x)
  }
}

# 从截断泊松分布中取一个随机样本
sample_from_truncated_poisson <- function() {
  prob_other <- 0.9  # 定义 prob_other 在函数内部
  
  while (TRUE) {
    x <- rpois(1, lambda = 1)
    u <- runif(1)
    prob_x <- truncated_poisson_pmf(x)
    if (x == 0) {
      if (u <= prob_x) {
        return(x)
      }
    } else {
      if (u <= prob_x / prob_other) {
        return(x)
      }
    }
  }
}

#检查行列数是否同时为0
check_zero_mutations <- function(l, h) {
  if (l == 0 && h == 0) {
    return(TRUE)  # 添加和删除的行列数同时为0，需要重新进行变异
  } else {
    return(FALSE) # 添加和删除的行列数不同时为0，可以进行下一步操作
  }
}

#变异，增加或删除基因或细胞，增加或删除的概率和行数有关
GABiclustering_mutation <- function(bicluster = B[[1]]){
  # 添加行的概率
  add_row_prob <- exp(- 0.1* length(bicluster[[1]]))
  add_col_prob <- exp(-0.06* length(bicluster[[3]]))
  #delete_col_prob <- 1 - exp(- 0.06* length(bicluster[[3]]))
  # 删除行的概率
  #delete_row_prob <- 1 - exp(- 0.1* length(bicluster[[1]]))
  # 从截断泊松分布中取一个随机样本
  l <- sample_from_truncated_poisson()
  h <- sample_from_truncated_poisson()
  if (check_zero_mutations(l, h)) {
    # 如果同时为0，则重新进行变异
    return(GABiclustering_mutation(bicluster))
  } else {
    if(runif(1)<add_row_prob){
      if (l != 0) {
        bicluster <- GABiclustering_gene_add(bicluster, l)
      }
    }
    else{
      if (l != 0) 
        bicluster <- GABiclustering_gene_remove(bicluster, l)
    }
    
    
    if(runif(1)<add_col_prob){
      if (h != 0){ 
        bicluster <- GABiclustering_cell_add(bicluster, h)
      }
    }
    else{
      if (h != 0){ 
        bicluster <- GABiclustering_cell_remove(bicluster, h)
      }
    }
    
    bicluster <- check_express(bicluster)
    return(bicluster)
  }
}

a <- GABiclustering_mutation(bicluster = B[[1]])



##在findall之前最终版的GABiclustering
#max_iteration，最大迭代次数
#N，初始种群个数
#m，初始种群中的最大基因个数
#n，初始种群中的最大细胞个数
#alpha，细胞面积权重
#beta，基因相关系数权重
#lambda，细胞距离权重
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
                           lambda = 1, beta = 1, good_num = 20,
                           gene_result_address = 'GABiclustering_gene_result_',
                           cell_result_address = 'GABiclustering_cell_result_',
                           result_num = 20, seed = 2021){
  set.seed(seed)
  same_times <- 0
  mid_num <- ceiling((N-good_num)/2)
  B <- GABiclustering_initialize(N,m,n)
  biclustering_gene_result_names <- paste(gene_result_address,alpha,'_',lambda,'_',beta,'_',seed,'.csv',sep = '')
  biclustering_cell_result_names <- paste(cell_result_address,alpha,'_',lambda,'_',beta,'_',seed,'.csv',sep = '')
  output_names <- paste('output_',alpha,'_',lambda,'_',beta,'_',seed,'.txt',sep = '')
  fitness_old <- GABiclustering_fitness(B, alpha, lambda, beta)
  sink(output_names)
  cat('Initial Score:',fitness_old$s_ordered[1:10],'\n')
  for(i in 1:max_iteration){
    fitness_old <- GABiclustering_fitness(B, alpha, lambda,beta)
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
    fitness_new <- GABiclustering_fitness(B_new, alpha, lambda, beta)
    new_index <- fitness_new$s_ordered_index[1:N]
    B <- B_new[new_index]
    fitness_new <- GABiclustering_fitness(B, alpha, lambda, beta)
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

getwd()
setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")

GAbiclustering_result <- GABiclustering(max_iteration = 1000, N = 100, m = 4, n = 9, alpha = 0.8, 
                                        lambda = 1, beta = 1, good_num = 20,
                                        gene_result_address = 'GABiclustering_gene_result_',
                                        cell_result_address = 'GABiclustering_cell_result_',
                                        result_num = 20, seed = 2021)

save(GAbiclustering_result,file=paste('GAbiclustering_result_',0.002,'_',1,'_',1,'_',2020,'.Rdata',sep = ''))



B <- GAbiclustering_result[[1]]