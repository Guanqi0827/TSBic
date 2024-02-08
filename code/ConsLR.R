##构造的时候考虑每隔几个时间点，多加一次（用BIC准则定最高次数）
##输入gene_index，基因文件的下标，tree，树的根节点名称，interval，每隔几个时间点多一次（多一个要求，最少
## 3次，10个点以下统一3次，10个点以上多几个点多一次）
##输出每个细胞对应的beta的向量的位置，用起始下标，和结束下标表示，同时记录样本个数num, 和beta长度beta_length


beta_index <- function(gene_index = 1, tree = 'C', interval = 5){
  li <- list()
  cell_name <- grep(tree, CELL, value = T)
  ind <- 1
  num <- 0
  beta_length <- 0
  not_na_num <- 0
  for(i in 1:length(cell_name)){
    flag <- which(CELL_blot[[gene_index]][,1] == cell_name[i])
    t <- NA
    q <- length(flag)
    if(q > 4){
      not_na_num <- not_na_num + 1
      if(q <= 10){
        k <- 3
      }
      else{
        k <- (q - 10) %/% interval + 3
      }
      t <- ind:(ind + k)
      ind <- ind + k + 1
      num <- num + q
      beta_length <- beta_length + k
    }
    li[[cell_name[i]]] <- t
  }
  back <- list(li, num, beta_length + not_na_num)
  names(back) <- c(tree, 'sample_size', 'beta_length')
  return(back)
}


Beta_index <- beta_index(gene_index = 1, tree = 'C', interval = 12)


##构造数据矩阵X
X <- list()

##输入gene_index，基因文件的下标， li，即beta_index函数返回的列表
##beta_length，beta向量的长度，是beta_index函数返回的列表中的'beta_length'，
##n，设计矩阵的行数（样本个数），是beta_index函数返回的列表中的'sample_size'，
##输出设计矩阵X，和响应变量y, 还有各细胞样本在X行中的起始和终止位置

designmat <- function(gene_index = 1, li = Beta_index){
  beta_length <- li[['beta_length']]
  n <- li[['sample_size']]
  tree <- names(li)[[1]]
  mat <- matrix(0, nrow = n, ncol = beta_length)
  y <- c()
  cell_name <- grep(tree, CELL, value = T)
  ind <- 0
  lis <- list()
  ind_x <- 0
  for(i in 1:length(cell_name)){
    flag <- which(CELL_blot[[gene_index]][,1] == cell_name[i])
    if(length(flag) > 4){
      b_index <- li[[1]][[cell_name[i]]]
      ctime <- CELL_blot[[gene_index]][flag,2]
      y <- append(y, CELL_blot[[gene_index]][flag,3])
      l <- length(ctime)
      for(j in 1:l){
        x <- ctime[j]
        for(k in 1:length(b_index)){
          mat[ind + j, b_index[k]] <- x ^ (k - 1)
        }
      }
      ind <- ind + l
      lis[[cell_name[i]]] <- c(ind_x + 1, ind_x + l)
      ind_x <- ind_x + l
    }
  }
  back <- list('X' = mat, 'y' = y, 'X_ind' = lis)
  return(back)
}

lii <- designmat(gene_index = 1, li = Beta_index)
X <- lii[['X']]
y <- lii[['y']]




##输入gene_index，基因文件的下标， li，即beta_index函数返回的列表
##输出约束条件矩阵H
CELL_root

constrainedmat <- function(gene_index = 1, li = Beta_index){
  beta_length <- li[['beta_length']]
  tree <- names(li)[[1]]
  cell_name <- grep(tree, CELL, value = T)
  r <- length(cell_name) - 1
  mat <- matrix(0, nrow = r, ncol = beta_length)
  ind <- 1
  for(i in 1:length(cell_name)){
    if(is.element(cell_name[i], CELL_root)){
      next
    }
    else{
      flag <- which(CELL_blot[[gene_index]][,1] == cell_name[i])
      if(length(flag) > 4){
        s <- substr(cell_name[i], 1, nchar(cell_name[i]) - 1)
        flg <- which(CELL_blot[[gene_index]][,1] == s)
        if(length(flg) > 1){
          x <- CELL_blot[[gene_index]][flg,2][length(flg)]
          b_index <- li[[1]][[s]]
          for(k in 1:length(b_index)){
            mat[ind, b_index[k]] <- x ^ (k - 1)
          }
          b_index <- li[[1]][[cell_name[i]]]
          for(k in 1:length(b_index)){
            mat[ind, b_index[k]] <- - x ^ (k - 1)
          }
          ind <- ind + 1
        }
      }
    }
    
  }
  ct <- 0
  for(i in 1:r){
    if(all(mat[i,] == 0)){
      ct <- 1
      break
    }
  }
  if(ct == 1){
    H <- mat[1:(i-1),,drop = F]
    return(H)
  }
  else{
    return(mat)
  }
}


H <- constrainedmat(gene_index = 1, li = Beta_index)


##求Beta的参数估计
##输入gene_index，基因文件的下标，tree，树的根节点名称，interval，每隔几个时间点多一次
##输出beta, 模型的拟合值，R方，BIC值

library(Matrix)
ConsLR <- function(gene_index = 1, tree = 'C', interval = 12){
  Beta_index <- beta_index(gene_index, tree, interval)
  li <- designmat(gene_index, Beta_index)
  l <- li$'X_ind'
  X <- li[['X']]
  X <- Matrix(X, sparse = T)
  y <- li[['y']]
  n <- Beta_index$sample_size
  beta_length <- Beta_index$beta_length
  H <- constrainedmat(gene_index, Beta_index)
  H <- Matrix(H, sparse = T)
  x_ <- solve(t(X) %*% X, tol = 1e-1000)
  beta_l <- x_ %*% t(X) %*% y
  if(all(H == 0)){
    beta_h <- beta_l
    d <- beta_length
  }
  else{
    h_ <- solve(H %*% x_ %*% t(H), tol = 1e-1000)
    beta_h <- beta_l - x_ %*% t(H) %*% h_ %*% H %*% beta_l
  }
  d <- beta_length - nrow(H)
  y_hat <- X %*% beta_h
  TSS <- var(y) * (n - 1)
  RSS <- sum((y - y_hat) ^ 2)
  r.squared <- 1 - RSS / TSS
  adj.r.squared <- 1 - (RSS/(n - d)) / (TSS/(n - 1))
  sigma_h <- RSS / (n - d)
  bic = 1/n * (RSS + log(n) * d * sigma_h)
  aic = 1/n * (RSS + 2 * d * sigma_h)
  back <- list(Beta_index[[1]], 'sample_size' = Beta_index[[2]],
               'beta_length' = Beta_index[[3]], 'X' = X, 'y' = y, 
               'H' = H, 'beta_h' = beta_h, 'x_ind' = l,  'y_hat' = y_hat, 
               'sigma_h' = sigma_h, 'r.squared' = r.squared, 
               'adj.r.squared' = adj.r.squared, 'BIC' = bic, 'AIC' = aic)
  names(back)[1] = tree
  return(back)
  
}

li <- ConsLR(gene_index = 1, tree = 'C', interval = 12)

X <- li[['X']]


##拟合每个基因的全部树，用bic准则来定次数

CELL_blot_polySpl <- list()

for(i in 1:length(GENE_names)){
  li_all <- list()
  for(j in 1:length(CELL_root)){
    bic = c()
    for(k in 4:20){
      li <- ConsLR(gene_index = i, tree = CELL_root[j], interval = k)
      bic <- append(bic, li$'BIC')
    }
    bic_min <- min(bic)
    interval <- which(bic %in% bic_min)[1] + 3
    pspl <- ConsLR(gene_index = i, tree = CELL_root[j], interval = interval)
    li_all[[CELL_root[j]]] <- pspl
  }
  CELL_blot_polySpl[[GENE_names[i]]] <- li_all
}

save(CELL_blot_polySpl, file = 'CELL_blot_polySpl.Rdata')
print(object.size(CELL_blot_polySpl), unit = 'MB')

##看一个基因，一个细胞分支用bic，aic准则的差别
j = 'E'
li_all <- list()
bic = c()
for(k in 4:10){
  li <- ConsLR(gene_index = i, tree = j, interval = k)
  bic <- append(bic, li$'BIC')
}
bic_min <- min(bic)
interval <- which(bic %in% bic_min)[1] + 3
pspl <- ConsLR(gene_index = i, tree = j, interval = interval)
li_all[['BIC']] <- pspl

aic = c()
for(k in 4:10){
  li <- ConsLR(gene_index = i, tree = j, interval = k)
  aic <- append(aic, li$'AIC')
}
aic_min <- min(aic)
interval <- which(aic %in% aic_min)[1] + 3
pspl <- ConsLR(gene_index = i, tree = j, interval = interval)
li_all[['AIC']] <- pspl





##输入细胞名称， 输出该细胞所在树的根细胞名称

find_root <- function(cell_name){
  for(i in 0:nchar(cell_name)){
    s <- substr(cell_name, 1, nchar(cell_name) - i)
    if(is.element(s, CELL_root)){
      break
    }
  }
  return(s)
}

find_root('E')

CELL_blot_polySpl <- list()
CELL_blot_polySpl[[1]] <- list()
CELL_blot_polySpl[[1]][['E']]<- li

##画图（单单画一个细胞）

draw_polySpl <- function(gene_index, cell_name){
  cell_root <- find_root(cell_name)
  mode <- CELL_blot_polySpl[[gene_index]][[cell_root]]
  flag <- which(CELL_blot[[gene_index]][,1] == cell_name)
  if(length(flag) > 4){
    ctime <- CELL_blot[[gene_index]][flag,2]
    ind <- which(names(mode$x_ind) %in% cell_name)
    st <- mode$x_ind[[ind]][1]
    ed <- mode$x_ind[[ind]][2]
    lab <- paste(GENE_names[gene_index], 'in', cell_name)
    plot(ctime, mode$y[st:ed], ylab = 'blot', xlab = 'time', main = lab)
    lines(ctime, mode$y_hat[st:ed], lwd = 2, col = 2)
  }
  else{
    return('NA')
  }
}

draw_polySpl(27, 'Eplpa')


##画R^2图(一棵数一个R方，总共174*5棵数)


r <- c()
for(i in 1:length(GENE_names)){
  for(j in 1:length(CELL_root)){
    mode <- CELL_blot_polySpl[[i]][[j]]
    if(mode$r.squared > 0){
      r <- c(r, mode$r.squared)
    }
    
  }
}

summary(r)
hist(r, col = NULL, main = 'Histogram of R^2', xlab = 'R^2', breaks = seq(0,1,0.05), ylim = c(0, 250))
sd(r)
mean(r)
length(which(r > 0.7))/869
length(which(r > 0.8))/869

