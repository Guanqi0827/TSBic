
##把每个基因文件按基因名称分组，同一个基因的不同文件当成一个组，按整体来聚类

GENE_group <- list()

##共104个不同的基因
names_GENE <- c()
for(i in 1:length(GENE_names)){
  ch <- GENE_names[i]
  x <- strsplit(ch, '_')[[1]][1]
  if(!is.element(x, names_GENE)){
    names_GENE <- c(names_GENE, x)
  }
}
##存每个基因的文件下标

for(i in 1:length(names_GENE)){
  ch <- names_GENE[i]
  x <- grep(ch, GENE_names)
  GENE_group[[ch]] <- x
}

###初始化随机矩阵

GABiclustering_initialize <- function(N = 10, m = 8, n = 12){
  
  B <- vector('list',N)
  m <- c(2:m)
  n <- c(2:n)
  for(i in 1:N){
    gene_name <- sample(names_GENE, sample(m,1))
    g <- Getgene_index(gene_name)
    dat <- Data_3[g,]
    c_expressed <- as.numeric(which(colMeans(dat) != 0))
    c_unexpressed <- as.numeric(which(colMeans(dat) == 0))
    num <- sample(n,1)
    if(length(c_expressed) >= num){
      c <- sample(c_expressed, num)
    }
    else{
      c <- c(c_expressed, sample(c_unexpressed, num-length(c_expressed)))
    }
    B[[i]][['gene index']] <- g
    B[[i]][['gene names']] <- gene_name
    B[[i]][['cell index']] <- c
    B[[i]][['cell names']] <- CELL[c]
  }
  
  return(B)
}

B <- GABiclustering_initialize(m = 20, n = 40)


#计算适应度，并排序，只按得分排序
GABiclustering_fitness <- function(B, alpha = 1, lambda = 1,beta = 1){
  fitness <- vector('list',length(B))
  s <- NULL
  c <- NULL
  d <- NULL
  a <- NULL
  for(i in 1:length(B)){
    li <- B[[i]]
    fitness[[i]] <- score(gene_name = li$`gene names`,cell_index = li$`cell index`, 
                          alpha = alpha, lambda = lambda, beta = beta)
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

GABiclustering_fitness(B)

##添加基因，按名称添加
GABiclustering_gene_add <- function(bicluster = B[[1]], l = 1){
  
  gene_name <- bicluster$`gene names`
  gene_add_index <- c(1:length(names_GENE))[-which(names_GENE %in% gene_name)]
  if(length(gene_add_index) == 0){
    return(bicluster)
  }
  else if(length(gene_add_index) >= 1 & length(gene_add_index) <= l){
    gene_newname <- names_GENE[gene_add_index]
    gene_newindex <- Getgene_index(gene_newname)
    gene_index <- c(bicluster[[1]], gene_newindex)
    gene_name <- c(bicluster[[2]], gene_newname)  
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_name
    return(bicluster)
  }
  else{
    gene_add_index <- sample(gene_add_index,l)
    gene_newname <- names_GENE[gene_add_index]
    gene_newindex <- Getgene_index(gene_newname)
    gene_index <- c(bicluster[[1]], gene_newindex)
    gene_name <- c(bicluster[[2]], gene_newname)  
    bicluster$`gene index` <- gene_index
    bicluster$`gene names` <- gene_name
    return(bicluster)
  }
  
}

li <- GABiclustering_gene_add(B[[1]])


##添加细胞

GABiclustering_cell_add <- function(bicluster = B[[1]], h = 1){
  
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index` 
  cell_add_index <- c(1:length(CELL))[-which(c(1:length(CELL)) %in% cell_index)]
  dat <- Data_3[gene_index,]
  c_expressed <- as.numeric(which(colMeans(dat) > 0))
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


##删除基因，按名称删除

GABiclustering_gene_remove <- function(bicluster = B[[1]], l = 1){
  
  gene_name <- bicluster$`gene names`
  if(length(gene_name) <= l+2){
    return(bicluster)
  }
  else{
    gene_remove_index <- sample(gene_name,l)
    gene_newname <- gene_name[-which(gene_name %in% gene_remove_index)]
    gene_newindex <- Getgene_index(gene_newname)
    bicluster$`gene index` <- gene_newindex
    bicluster$`gene names` <- gene_newname
    return(bicluster)
  }
  
}

GABiclustering_gene_remove()

##删除细胞
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


GABiclustering_cell_remove()

##检查要求表达占百分之60以上
##如果删除后数量太少（<3），会按生成双聚类的方式随机打乱，重新得到一个初始双聚类
check_express <- function(bicluster, pro = 3/5){
  gene_index <- bicluster$`gene index`
  cell_index <- bicluster$`cell index`
  gene_name <- bicluster$`gene names`
  l_g <- length(gene_name)
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
      gene_drop_index <- gene_index[ind]
      gene_drop_name <- GENE_names[gene_drop_index]
      gene_drop_name <- strsplit(gene_drop_name, '_')[[1]][1]
      gene_drop_index <- which(gene_name %in% gene_drop_name)
      gene_name <- gene_name[-gene_drop_index]
    }
  }
  if(length(gene_name) < 2 || length(cell_index) < 2){
    gene_name <- sample(names_GENE, l_g)
    gene_index <- Getgene_index(gene_name)
    dat <- Data_3[gene_index,]
    c_expressed <- as.numeric(which(colMeans(dat) != 0))
    if(length(c_expressed) >= l_c){
      cell_index <- sample(c_expressed, l_c)
    }
    else{
      cell_index <- sample(1:length(CELL), l_c)
    }
    bicluster <- list('gene index' = gene_index, 'gene names' = gene_name,
                      'cell index' = cell_index, 'cell names' = CELL[cell_index])
    return(bicluster)
  }
  else{
    gene_index <- Getgene_index(gene_name)
    bicluster <- list('gene index' = gene_index, 'gene names' = gene_name,
                      'cell index' = cell_index, 'cell names' = CELL[cell_index])
    bicluster <- check_express(bicluster, pro)
    return(bicluster)
  }
}

li <- GABiclustering_initialize(100, 20, 30)
for(i in 1:100){
  B <-li[[1]]
  B <- check_express(B)
}
B <- GABiclustering_initialize(10, 20, 30)[[1]]
drawBiclusterExpressed(B)
B <- check_express(B)
drawBiclusterExpressed(B)

B <- All_B[[19]]
check_express(B)

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

#变异，增加或删除基因或细胞，增加或删除的概率和面积有关
#area_score为双聚类面积分数
#add_p_coef调整面积分数转为0-1的权重
GABiclustering_mutation <- function(bicluster = B[[1]]){
  
  #area <- length(bicluster[[1]])*length(bicluster[[3]])
  #add_p <- exp(-add_p_coef*area)
  add_row_prob <- exp(- 0.1* length(bicluster[[1]]))
  add_col_prob <- exp(-0.06* length(bicluster[[3]]))
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


#交叉互换，将部分好的基因随机给其他组


GABiclustering_gene_exchange <- function(badbicluster = B[[2]], goodbicluster = B){
  
  samp <- sample(length(goodbicluster),1)
  goodbicluster <- goodbicluster[[samp]]
  good_gene_name <- goodbicluster$`gene names`
  bad_gene_name <- badbicluster$`gene names`
  
  k <- 0
  while (k == 0) {
    k <- min(rpois(1, lambda = 1), min(nrow(good_gene_name), nrow(bad_gene_name)))
  }
  
  if(length(bad_gene_name) <= 3*k | length(good_gene_name) <= 3*k){
    badbicluster <- GABiclustering_mutation(badbicluster)
    return(badbicluster)
  }
  else{
    good_gene_exchange_index <- sample(1:length(good_gene_name), k)
    add_value <- good_gene_name[good_gene_exchange_index]
    if(length(bad_gene_name) >= max(good_gene_exchange_index)){
      bad_gene_name <- c(bad_gene_name[-good_gene_exchange_index],add_value)
      bad_gene_name <- unique(bad_gene_name)
      gene_index <- Getgene_index(bad_gene_name)  
      badbicluster$`gene index` <- gene_index
      badbicluster$`gene names` <- bad_gene_name
      return(badbicluster)
    }
    else{
      remove_index <- sample(1:length(bad_gene_name),k)
      bad_gene_name <- c(bad_gene_name[-remove_index],add_value)
      bad_gene_name <- unique(bad_gene_name)
      gene_index <- Getgene_index(bad_gene_name)  
      badbicluster$`gene index` <- gene_index
      badbicluster$`gene names` <- bad_gene_name
      return(badbicluster)
    }
  }
}

li <- GABiclustering_initialize(10, 20, 30)
B <- li[[1]]
GABiclustering_gene_exchange(B, li[2])


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

li <- GABiclustering_initialize(10, 20, 30)
B <- li[[1]]
GABiclustering_cell_exchange(B, li[2])

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

##在findall之前最终版的GABiclustering，因为跑程序的时间很长，为了
##意外中止后，能从上次停止的地方继续跑，多一个参数continue，
##continue默认为0，从新开始，即随机初始化B，如为1，则直接取当前B

GABiclustering <- function(max_iteration = 1000, N = 100, m = 4, n = 9, alpha = 0.8, 
                           lambda = 1, beta = 1, good_num = 20,
                           continue = 0, seed = 2021){
  set.seed(seed)
  same_times <- 0
  mid_num <- ceiling((N-good_num)/2)
  save_names <- paste('B_',alpha,'_',lambda,'_',beta,'_',seed,'.Rdata',sep = '')
  if(continue == 0){
    B <- GABiclustering_initialize(N, m, n)
  }
  else{
    load(save_names)
  }
  output_names <- paste('output_',alpha,'_',lambda,'_',beta,'_',seed,'.txt',sep = '')
  fitness_old <- GABiclustering_fitness(B, alpha, lambda, beta)
  sink(output_names)
  cat('Initial Score:',fitness_old$s_ordered[1:5],'\n')
  for(i in 1:max_iteration){
    fitness_old <- GABiclustering_fitness(B, alpha, lambda,beta)
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
    fitness_new <- GABiclustering_fitness(B_new, alpha, lambda, beta)
    new_index <- fitness_new$s_ordered_index[1:N]
    B <- B_new[new_index]
    fitness_new <- GABiclustering_fitness(B, alpha, lambda, beta)
    nice_index_new <- fitness_new$s_ordered_index[1:5]
    nice_b_new <- B[nice_index_new]
    if(identical(nice_b_old, nice_b_new)){
      same_times <- same_times+1
    }
    else{ same_times <- 0 }
    if(same_times == 2){
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



getwd()
setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")

GAbiclustering_result <- GABiclustering(max_iteration = 1000, N = 100, m = 10, n = 15, alpha = 0.8, 
                                        lambda = 1, beta = 1, good_num = 20, seed = 2021)
save(GAbiclustering_result,file=paste('GAbiclustering_result_',0.002,'_',1,'_',1,'_',2020,'.Rdata',sep = ''))


B <- GAbiclustering_result[[1]]

