
##单独算每个细胞的R方

CELL_fate_blot_r <- list()


for(i in 1:length(GENE_names)){
  li <- CELL_blot_polySpl[[i]]
  li_all <- list()
  r.squared <- NA
  for(j in 1:length(CELL)){
    cell_name <- CELL[j]
    s <- substr(cell_name, 1, 1)
    id <- grep(s, names(li))
    lii <- li[[id]]
    ind <- which(names(lii$'x_ind') %in% cell_name)
    if(length(ind) == 1){
      st <- lii$'x_ind'[[ind]][1]
      ed <- lii$'x_ind'[[ind]][2]
      y <- lii$y[st:ed]
      y_hat <- lii$y_hat[st:ed]
      n <- ed - st + 1
      TSS <- var(y) * (n - 1)
      RSS <- sum((y - y_hat) ^ 2)
      r.squared <- 1 - RSS / TSS
    }
    li_all[[cell_name]] <- r.squared
  }
  CELL_fate_blot_r[[GENE_names[i]]] <- li_all
}


##只算有表达的细胞的R方，且只存大于0的R^2，画R^2图，看看效果，效果不好，不画了
##直接用全部数的R方的均值和中位数就行

r <- c()
for(i in 1:length(GENE_names)){
  for(j in 1:length(CELL)){
    if(Data_3[i, j] == 1){
      a <- CELL_fate_blot_r[[i]][[j]]
      if(!is.na(a) && a > 0){
        r <- c(r, a)
      }
    }
  }
}

summary(r)
hist(r)

##存fate中细胞blot的系数向量和次数向量


CELL_fate_blot_polySpl <- list()

for(i in 1:length(GENE_names)){
  li <- CELL_blot_polySpl[[i]]
  li_all <- list()
  for(j in 1:length(CELL)){
      l <- list()
      cell_name <- CELL[j]
      s <- substr(cell_name, 1, 1)
      id <- grep(s, names(li))
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
  CELL_fate_blot_polySpl[[GENE_names[i]]] <- li_all
}

save(CELL_fate_blot_polySpl, file = 'CELL_fate_blot_polySpl.Rdata')

save(CELL_blot, file = 'CELL_blot.Rdata')

####blot模型中的多项式函数求导
##输入多项式的系数向量和次数向量（列表的形式）
##输出多项式求导后的系数向量和次数向量


derivative <- function(mode){
  if(length(mode) > 0){
    coef <- mode$coefficients
    if(length(coef) > 1){
      len <- length(coef)
      coef_new <- rep(0, len - 1)
      for(i in 2:len){
        coef_new[i - 1] <- coef[i] * (i - 1)
      }
      degree_new <- 0:(len - 2)
    }
    else{
      coef_new <- 0
      degree_new <- 0
    }
    back <- list('coefficients' = coef_new, 'degree' = degree_new)
    return(back)
  }
  else{
    back <- list()
  }
}
#######
#######对全部模型求导
CELL_fate_diff_polySpl_deri <- list()

for(i in 1:length(GENE_names)){
  li <- CELL_fate_blot_polySpl[[i]]
  li_all <- list()
  for(j in 1:length(CELL)){
    cell_name <- CELL[j]
    l <- derivative(li[[j]])
    li_all[[cell_name]] <- l
  }
  CELL_fate_diff_polySpl_deri[[GENE_names[i]]] <- li_all
}

save(CELL_fate_diff_polySpl_deri, file = 'CELL_fate_diff_polySpl_deri.Rdata')

print(object.size(CELL_fate_blot_polySpl), unit = 'MB')
print(object.size(CELL_fate_diff_polySpl_deri), unit = 'MB')


####导数模型中的多项式函数平移缩放定义域到[0, 1]之间
##输入导数多项式的系数向量和次数向量（列表的形式）,基因下标，细胞名称
##输出导数多项式平移缩放后的系数向量和次数向量

deriv_scale <- function(mode, gene_index, cell_name){
  
  if(length(mode) > 0){
    flag <- which(CELL_blot[[gene_index]][,1] == cell_name)
    ctime <- CELL_blot[[gene_index]][flag,2]
    t <- length(ctime) - 1
    st <- ctime[1] 
    coef <- mode$coefficients
    len <- length(coef)
    coef_new <- rep(0, len)
    degree_new <- 0:(len - 1)
    for(i in 1:len){
      for(j in i:len){
        coef_new[i] <- coef_new[i] + coef[j] * choose(j - 1, i - 1) * st ^ (j - i) * t ^ (i - 1)
      }
    }
    back <- list('coefficients' = coef_new, 'degree' = degree_new)
  }
  else{
    back <- list()
  }
  return(back)
}

li <- deriv_scale(CELL_fate_diff_polySpl_deri[[1]][[1]], 1, 'Cpppaa')

#######对全部导数模型进行平移缩放

##

CELL_fate_diff_polySpl_deri_scale <- list()

for(i in 1:length(GENE_names)){
  li_all <- list()
  for(j in 1:length(CELL)){
    mode <- CELL_fate_diff_polySpl_deri[[i]][[j]]
    li <- deriv_scale(mode, i, CELL[j])
    li_all[[CELL[j]]] <- li
  }
  CELL_fate_diff_polySpl_deri_scale[[GENE_names[i]]] <- li_all
}

save(CELL_fate_diff_polySpl_deri_scale, file = 'CELL_fate_diff_polySpl_deri_scale.Rdata')








score(1:8, 34:46)
score(9:12, 134:140)
gene_index <- All_B[[2]]$`gene index`
cell_index <- All_B[[2]]$`cell index`
alpha <- 1.4
lambda <- 1.1
beta <- 1
biclusterplot(sort(gene_index),sort(cell_index),alpha,lambda,beta)





