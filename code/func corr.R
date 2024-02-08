
##模型中的多项式函数求导，并且缩放定义域到[0, 1]区间上，
##mode：线性回归模型, T：原区间长度，输出为多项式函数的系数和次数

deriv_scale <- function(mode, T){
  coef <- as.numeric(mode$coefficients)
  mode_term <- mode$terms
  if(attributes(mode_term)$intercept == 1){
    coef <- coef[-1]
  }
  len <- length(coef)
  coef_new <- rep(0, len)
  for(i in 1:len){
    if(!is.na(coef[i])){
      coef_new[i] <- coef[i] * i
    }
  }
  degree_new <- 0:(len - 1)
  for(i in 1:len){
    coef_new[i] <- coef_new[i] * T ^ (i - 1)
  }
  back <- list('coefficients' = coef_new, 'degree' = degree_new)
  return(back)
}

##存所有模型求导并缩放后的系数向量和次数向量

CELL_fate_diff_polySpl_deri <- list()
for(i in 1:length(GENE_names)){
  a <- CELL_fate_blot_polySpl[[i]]
  li <- list()
  for(j in 1:length(CELL)){
    b <- a[[j]]
    lii <- list()
    if(length(b) > 0){
      l <- b[[3]]
      der <- deriv_scale(b[[1]], l)
      lii[['coefficients']] <- der$'coefficients'
      lii[['degree']] <- der$'degree'
      li[[CELL[j]]] <- lii
    }
    else{
      li[[CELL[j]]] <- lii
    }
  }
  CELL_fate_diff_polySpl_deri[[GENE_names[i]]] <- li
}

save(CELL_fate_diff_polySpl_deri, file = 'CELL_fate_diff_polySpl_deri.Rdata')

##两个多项式函数的乘积，输入和输出均是系数向量和次数向量

Poly_multiple <- function(coef1, degree1, coef2, degree2){
  l1 <- length(degree1)
  l2 <- length(degree2)
  coef <- rep(0, l1 + l2 - 1)
  degree <- 0:(l1 + l2 - 2)
  for(i in 1:l1){
    for(j in 1:l2){
      d <- degree1[i] + degree2[j]
      coef[d+1] <- coef[d+1] + coef1[i] * coef2[j]
    }
  }
  back <- list('coefficients' = coef, 'degree' = degree)
  return(back)
}
##验证正确性
a <- CELL_fate_diff_polySpl_deri[[1]][[1]]
b <- CELL_fate_diff_polySpl_deri[[1]][[2]]
Poly_multiple(a[[1]], a[[2]], b[[1]], b[[2]])

##在[0,1]上的多项式函数的定积分，输入是系数向量和次数向量，输出是定积分值（实数）

Poly_integral <- function(coef, degree){
  integral <- 0 
  l <- length(degree)
  for(i in 1:l){
    integral <- integral + coef[i] / (degree[i] + 1)
  }
  return(integral)

}
##验证正确性
coef <- c(1, 1, 4)
degree <- c(0, 1, 3)
b <- function(x) 1 + 1 * x + 4 * x ^ 3
a <- integrate(b, 0, 1)
Poly_integral(coef, degree)
  

##计算两个多项式函数的相关系数，输入两个多项式函数的系数向量和次数向量，输出相关系数（实数）

corr_poly <- function(coef1, degree1, coef2, degree2){
  pol <- Poly_multiple(coef1, degree1, coef2, degree2)
  coef <- pol$'coefficients'
  degree <- pol$'degree'
  ex <- Poly_integral(coef1, degree1)
  ey <- Poly_integral(coef2, degree2)
  cov_poly <- Poly_integral(coef, degree) - ex * ey
  pol <- Poly_multiple(coef1, degree1, coef1, degree1)
  coef <- pol$'coefficients'
  degree <- pol$'degree'
  varx <- Poly_integral(coef, degree) - ex ^ 2
  pol <- Poly_multiple(coef2, degree2, coef2, degree2)
  coef <- pol$'coefficients'
  degree <- pol$'degree'
  vary <- Poly_integral(coef, degree) - ey ^ 2
  cor_poly <- cov_poly / (sqrt(varx) * sqrt(vary))
  if(is.infinite(cor_poly) || is.nan(cor_poly)){
    return(NA)
  }
  else{
    return(cor_poly)
  }

}

##验证正确性
coef1 <- CELL_fate_diff_polySpl_deri[["pha-4_2"]][["Ealpp"]][[1]]
degree1 <- CELL_fate_diff_polySpl_deri[["pha-4_2"]][["Ealpp"]][[2]]
coef2 <- CELL_fate_diff_polySpl_deri[["B0310.2_4"]][["Ealpp"]][[1]]
degree2 <- CELL_fate_diff_polySpl_deri[["B0310.2_4"]][["Ealpp"]][[2]]

a <- CELL_fate_diff_polySpl[[1]][[45]][[2]]
b <- CELL_fate_diff_polySpl[[1]][[40]][[2]]
cor(a, b)
corr_poly(coef1, degree1, coef2, degree2)  
  
##计算双聚类的相关系数分值
##注意：方差为0时，相关系数为Inf，应该删除，不参与计算


gene_corr_poly <- function(gene_index, cell_index){
  gtotal <- NULL
  for(j in 1:(length(gene_index) - 1)){
    g_ind_1 <- gene_index[j]
    for(k in (j + 1):length(gene_index)){
      g_ind_2 <- gene_index[k]
      gcorr <- NULL
      for(i in 1:length(cell_index)){
        ind <- cell_index[i]
        g_1 <- CELL_fate_diff_polySpl_deri_scale[[g_ind_1]][[ind]]
        g_2 <- CELL_fate_diff_polySpl_deri_scale[[g_ind_2]][[ind]]
        if(length(g_1) > 0 && length(g_2) > 0){
          coef1 <- g_1[[1]]
          degree1 <- g_1[[2]]
          coef2 <- g_2[[1]]
          degree2 <- g_2[[2]]
          gcorr_ <- corr_poly(coef1, degree1, coef2, degree2) 
          gcorr <- append(gcorr, gcorr_)
        }
      }
      gcorr <- na.omit(gcorr)
      if(length(gcorr) > 0){
        gcorr <- mean(gcorr)
        gtotal <- append(gtotal, gcorr)
      }
    }
  }
  if(!is.null(gtotal)){
    gcorr_score <- mean(na.omit(gtotal))
  }
  else{
    gcorr_score <- 0
  }
  back <- list('Gene Names'=GENE_names[gene_index],'Cell Names'=CELL[cell_index],
               'Gcorr' = gtotal,'Gcorr_Score'= gcorr_score)
  return(back)
}

##验证正确性
gene_corr_poly(1:8, 34:46)
gene_corr_poly(9:12, 134:140)
gene_corr_poly(gene_index, cell_index)



