

B <- B[[1]]
B <- All_B[[2]]
m <- length(B$`gene names`)
n <- length(B$`cell index`)
###产生一定维度的随机矩阵，和初始化双聚类步骤类似

random_initialize <- function(N = 10, m = 8, n = 12){
  
  B <- vector('list',N)
  for(i in 1:N){
    gene_name <- sample(names_GENE, m)
    g <- Getgene_index(gene_name)
    dat <- Data_3[g,]
    c_expressed <- as.numeric(which(colMeans(dat) != 0))
    if(length(c_expressed) >= n){
      c <- sample(c_expressed, n)
    }
    else{
      c <- sample(1:length(CELL), n)
    }
    B[[i]][['gene index']] <- g
    B[[i]][['gene names']] <- gene_name
    B[[i]][['cell index']] <- c
    B[[i]][['cell names']] <- CELL[c]
  }
  return(B)
}

###产生随机维度的随机矩阵，和初始化双聚类步骤类似

random_initialize <- function(N = 10){
  
  B <- vector('list',N)
  for(i in 1:N){
    m <- sample(2:20, 1)
    n <- sample(2:200, 1)
    gene_name <- sample(names_GENE, m)
    g <- Getgene_index(gene_name)
    dat <- Data_3[g,]
    c_expressed <- as.numeric(which(colMeans(dat) != 0))
    if(length(c_expressed) >= n){
      c <- sample(c_expressed, n)
    }
    else{
      c <- sample(1:length(CELL), n)
    }
    B[[i]][['gene index']] <- g
    B[[i]][['gene names']] <- gene_name
    B[[i]][['cell index']] <- c
    B[[i]][['cell names']] <- CELL[c]
  }
  return(B)
}




m <- 3
n <- 14
a <- random_initialize(N = 100)
B <- a
s <- NULL
fitness <- vector('list',length(B))
for(i in 1:length(B)){
  li <- B[[i]]
  fitness[[i]] <- score(gene_name = li$`gene names`,cell_index = li$`cell index`, 
                        alpha = 0.35, lambda = 0.15, beta = 0.6, delta = 0.04,
                        method = 'average', All_gene_index, All_cell_index)
  s <- append(s,fitness[[i]]$'Score')
}

li <- GABiclustering_fitness(a, alpha = 0.35, lambda = 0.15, beta = 0.6, delta = 0.04,
                             All_gene_index = All_gene_index, All_cell_index = All_cell_index)
li$s_ordered
random_score <- li

setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")

mean(li$s_ordered)
save(random_score, file = 'random_score2.Rdata')




