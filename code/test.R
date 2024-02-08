###取出每个细胞(已知fate的细胞)各时间点的原始数据值，进行随机性检验


CELL_fate_blot <- vector('list',174)

for(i in 1:174){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  cellname_ <- NULL
  time_ <- NULL
  blot_ <- NULL
  for(j in 1:145){
    index <- which(realdata[[ind]][[1]] %in% CELL[j])
    if(length(index) == 0){
    }
    else{
      cellname_  <- append(cellname_ ,rep(CELL[j],length(index)))
      time_ <- append(time_,realdata[[ind]][[2]][index])
      blot_  <- append(blot_ ,realdata[[ind]][[3]][index])
    }
  }
  CELL_fate_blot[[i]] <- data.frame('cell' = cellname_, 'time' = time_, 'blot' = blot_)
}


#取出每个细胞(已知fate的细胞)各时间点的差分值(全部增量)，进行随机性检验


CELL_fate_diff <- vector('list',174)

for(i in 1:174){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  cellname_ <- NULL
  time_ <- NULL
  diff_ <- NULL
  for(j in 1:145){
    index <- which(realdata[[ind]][[1]] %in% CELL[j])
    if(length(index) == 0){
    }
    else if(length(index) == 1){
      cellname_ <- append(cellname_, CELL[j])
      time_ <- append(time_,realdata[[ind]][[2]][index])
      diff_ <- append(diff_,realdata[[ind]][[3]][index])
    }
    else{
      d <- diff(realdata[[ind]][[3]][index])
      cellname_  <- append(cellname_ ,rep(CELL[j],length(d)))
      time_ <- append(time_,realdata[[ind]][[2]][index][-1])
      diff_  <- append(diff_ ,d)
    }
  }
  CELL_fate_diff[[i]] <- data.frame('cell' = cellname_, 'time' = time_, 'diff blot' = diff_)
}


##只取各个基因中有表达的细胞，算各时间点的blot值，准备做假设检验

CELL_blot_expressed <- vector('list',174)

for(i in 1:length(GENE_names)){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  cellname_ <- NULL
  time_ <- NULL
  blot_ <- NULL
  for(j in 1:length(CELL)){
    if(Data_3[i,j] == 1){
      index <- which(realdata[[ind]][[1]] %in% CELL[j])
      cellname_  <- append(cellname_ ,rep(CELL[j],length(index)))
      time_ <- append(time_,realdata[[ind]][[2]][index])
      blot_  <- append(blot_ ,realdata[[ind]][[3]][index])
    }
  }
  CELL_blot_expressed[[i]] <- data.frame('cell' = cellname_, 'time' = time_, 'blot' = blot_)
}


##只取各个基因中有表达的细胞，算各时间点的diff值，准备做假设检验


CELL_diff_expressed <- vector('list',174)


for(i in 1:length(GENE_names)){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  cellname_ <- NULL
  time_ <- NULL
  diff_ <- NULL
  for(j in 1:length(CELL)){
    if(Data_3[i,j] == 1){
      index <- which(realdata[[ind]][[1]] %in% CELL[j])
      if(length(index) == 1){
        diff_  <- append(diff_ ,realdata[[ind]][[3]][index])
        cellname_  <- append(cellname_ ,CELL[j])
        time_ <- append(time_,realdata[[ind]][[2]][index])
      }
      else{
        d <- diff(realdata[[ind]][[3]][index])
        cellname_  <- append(cellname_ ,rep(CELL[j],length(d)))
        time_ <- append(time_,realdata[[ind]][[2]][index][-1])
        diff_  <- append(diff_ ,d)
      }
      
    }
  }
  CELL_diff_expressed[[i]] <- data.frame('cell' = cellname_, 'time' = time_, 'diff blot' = diff_)
}


##取各个基因中的细胞，算各时间点的diff值，准备拟合

CELL_diff <- vector('list',174)

for(i in 1:length(GENE_names)){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  cellname_ <- NULL
  time_ <- NULL
  diff_ <- NULL
  for(j in 1:length(CELL)){
    index <- which(realdata[[ind]][[1]] %in% CELL[j])
    if(length(index) == 1){
      diff_  <- append(diff_ ,realdata[[ind]][[3]][index])
      cellname_  <- append(cellname_ ,CELL[j])
      time_ <- append(time_,realdata[[ind]][[2]][index])
    }
    else{
      d <- diff(realdata[[ind]][[3]][index])
      cellname_  <- append(cellname_ ,rep(CELL[j],length(d)))
      time_ <- append(time_,realdata[[ind]][[2]][index][-1])
      diff_  <- append(diff_ ,d)
    }
  }
  CELL_diff[[i]] <- data.frame('cell' = cellname_, 'time' = time_, 'diff blot' = diff_)
}


##取各个基因中的细胞，算各时间点的blot值，准备拟合

CELL_blot<- vector('list',174)

for(i in 1:length(GENE_names)){
  k <- seq(1,184)[-gene_missing_index][i]
  ind <- which(file_names_csv%in%file_gene[k,1])
  cellname_ <- NULL
  time_ <- NULL
  blot_ <- NULL
  for(j in 1:length(CELL)){
    index <- which(realdata[[ind]][[1]] %in% CELL[j])
    cellname_  <- append(cellname_ ,rep(CELL[j],length(index)))
    time_ <- append(time_,realdata[[ind]][[2]][index])
    blot_  <- append(blot_ ,realdata[[ind]][[3]][index])
  }
  CELL_blot[[i]] <- data.frame('cell' = cellname_, 'time' = time_, 'blot' = blot_)
}

##
# 游程检验

runtest_res <- vector('list',174)
for(i in 1:174){
  cellname_ <- NULL
  p_ <- NULL
  for(j in 1:145){
    flag <- which(CELL_fate_diff[[i]][,1] == CELL[j])
    if(length(flag) >= 3){
      cellname_ <- append(cellname_, CELL[j])
      p_ <- append(p_, round(runs.test(CELL_fate_diff[[i]][flag,3])$p.value,3))
    }
    else{
      cellname_ <- append(cellname_, CELL[j])
      p_ <- append(p_, NA)
    }
  }
  runtest_res[[i]] <- data.frame('cell' = cellname_, 'runtest_pvalue' = p_)
}

plot(1,1,xlim = c(1,145),ylim = c(0,1),type = 'n',xlab = 'cell',ylab = 'pvalue',xaxt = 'n')

#runtest_percent 为p值小于0.05的占比（174个基因分开算）
runtest_percent <- NULL
counts <- 0
total <- 0
for(i in 1:174){
  for(j in 1:145){
    if(!is.na(runtest_res[[i]][j,2])){
      total <- total + 1
      points(j,runtest_res[[i]][j,2],pch = 4, cex= 0.5)
      if(runtest_res[[i]][j,2] <= 0.05){
        counts <- counts + 1
      }
    }
  }
  runtest_percent <- append(runtest_percent,round(1-counts/total,3))
}
abline(h = 0.05, col = 'red')

runs.test(rnorm(100),plot.it = TRUE)
runs.test(CELL_fate_diff[[i]][flag,3],plot.it = TRUE)

# Ljung-Box 检验(自相关检验)

library(tseries)

boxtest_res <- vector('list',174)
for(i in 1:174){
  cellname_ <- NULL
  p_ <- NULL
  for(j in 1:145){
    flag <- which(CELL_fate_diff[[i]][,1] == CELL[j])
    if(length(flag) >= 3){
      cellname_ <- append(cellname_, CELL[j])
      p_ <- append(p_, round(Box.test(CELL_fate_diff[[i]][flag,3],lag = 1, type = 'Ljung-Box')$p.value,3))
    }
    else{
      cellname_ <- append(cellname_, CELL[j])
      p_ <- append(p_, NA)
    }
  }
  boxtest_res[[i]] <- data.frame('cell' = cellname_, 'boxtest_pvalue' = p_)
}

#boxtest_percent 为p值小于0.05的占比（174个基因分开算）
boxtest_percent <- NULL
counts <- 0
total <- 0
for(i in 1:174){
  for(j in 1:145){
    if(!is.na(boxtest_res[[i]][j,2])){
      total <- total + 1
      points(j,boxtest_res[[i]][j,2],pch = 4, cex= 0.5)
      if(boxtest_res[[i]][j,2] <= 0.05){
        counts <- counts + 1
      }
    }
  }
  boxtest_percent <- append(boxtest_percent,round(1-counts/total,3))
}
abline(h = 0.05, col = 'red')

Box.test(rnorm(100),lag = 1, type = 'Ljung-Box')

acf(CELL_fate_diff[[i]][flag,3])

Box.test(CELL_fate_diff[[i]][flag,3],lag = 1, type = 'Ljung-Box')$p.value


# Dickey-Fuller检验(单位根检验，看是否平稳)
library(fUnitRoots)
library(tseries)

adftest_res <- vector('list',174)
for(i in 1:174){
  cellname_ <- NULL
  p_ <- NULL
  for(j in 1:145){
    flag <- which(CELL_fate_diff[[i]][,1] == CELL[j])
    if(length(flag) >= 5){
      cellname_ <- append(cellname_, CELL[j])
      p_ <- append(p_, round(adfTest(CELL_fate_diff[[i]][flag,3])@test$p.value,3))
    }
    else{
      cellname_ <- append(cellname_, CELL[j])
      p_ <- append(p_, NA)
    }
  }
  adftest_res[[i]] <- data.frame('cell' = cellname_, 'adftest_pvalue' = p_)
}

adfTest(CELL_fate_diff[[i]][flag,3])
adf.test(CELL_fate_diff[[i]][flag,3])

#adftest_percent 为p值小于0.05的占比（174个基因分开算）
adftest_percent <- NULL
counts <- 0
total <- 0
for(i in 1:174){
  for(j in 1:145){
    if(!is.na(adftest_res[[i]][j,2])){
      total <- total + 1
      points(j,adftest_res[[i]][j,2],pch = 4, cex= 0.5)
      if(adftest_res[[i]][j,2] <= 0.05){
        counts <- counts + 1
      }
    }
  }
  adftest_percent <- append(adftest_percent,round(1-counts/total,3))
}
abline(h = 0.05, col = 'red')


###将以上操作包装成函数


##一起检验


# test_pvalue <- function(d = CELL_diff_expressed){
#   res <- vector('list',174)
#   for(i in 1:length(GENE_names)){
#     cellname_ <- NULL
#     p_run <- NULL
#     p_lb <- NULL
#     p_adf <- NULL
#     for(j in 1:length(CELL)){
#       if(Data_3[i,j] == 1){
#         flag <- which(d[[i]][,1] == CELL[j])
#         if(length(flag) >= 5){
#           cellname_ <- append(cellname_, CELL[j])
#           p_run <- append(p_run, round(runs.test(d[[i]][flag,3])$p.value,3))
#           p_lb <- append(p_lb, round(Box.test(d[[i]][flag,3],lag = 1, type = 'Ljung-Box')$p.value,3))
#           p_adf <- append(p_adf, round(adfTest(d[[i]][flag,3])@test$p.value,3))
#         }
#         else{
#           cellname_ <- append(cellname_, CELL[j])
#           p_run <- append(p_run, NA)
#           p_lb <- append(p_lb, NA)
#           p_adf <- append(p_adf, NA)
#         }
#       }
#       res[[i]] <- data.frame('cell' = cellname_, 'runtest_pvalue' = p_run,
#                              'lbtest_pvalue' = p_lb, 'adftest_pvalue' = p_adf)
#     }
#   }
#   return(res)
# }

library(lawstat)
test_blot_expressed <- vector('list',174)
test_blot_expressed <- test_pvalue(d = CELL_blot_expressed)
test_diff_expressed <- vector('list',174)
test_diff_expressed <- test_pvalue(d = CELL_diff_expressed)


# test_pvalue_percent <- function(d = test_blot_expressed){
#   runtest_percent <- NULL
#   lbtest_percent <- NULL
#   adftest_percent <- NULL
# 
#   for(i in 1:length(GENE_names)){
#     counts_run <- 0
#     counts_lb <- 0
#     counts_adf <- 0
#     total <- 0
#     for(j in 1:length(CELL)){
#       if(!is.na(d[[i]][j,2]) & !is.na(d[[i]][j,3]) & !is.na(d[[i]][j,4])){
#         total <- total + 1
#         if(d[[i]][j,2] >= 0.05){
#           counts_run <- counts_run + 1
#         }
#         if(d[[i]][j,3] >= 0.05){
#           counts_lb <- counts_lb + 1
#         }
#         if(d[[i]][j,4] <= 0.05){
#           counts_adf <- counts_adf + 1
#         }
#       }
#     }
#     runtest_percent <- append(runtest_percent,round(counts_run/total,3))
#     lbtest_percent <- append(lbtest_percent,round(counts_lb/total,3))
#     adftest_percent <- append(adftest_percent,round(counts_adf/total,3))
#   }
#   res <- list('runtest_percent' = runtest_percent, 'lbtest_percent' = lbtest_percent,
#               'adftest_percent' = adftest_percent)
#   return(res)
# }

test_percent_blot_expressed <- test_pvalue_percent(test_blot_expressed)
test_percent_diff_expressed <- test_pvalue_percent(test_diff_expressed)


###用上细胞内的全部数据，用两细胞内数据的假设检验的p值作为细胞间距离
##其中包括KS检验，Wilcoxon检验，Mann-Whitney检验.

ks.test(x,y,alternative = 'two.sided')
wilcox.test(x, y, paired = F, alternative = 'two.sided', exact = T)


CELL_fate_diff_pdist <- matrix(nrow = length(CELL), ncol = length(CELL))
for(j in 1:length(CELL)){
  for(k in j:length(CELL)){
    pvalue <- NULL
    for(i in 1:length(GENE_names)){
      flag_1 <- which(CELL_fate_diff[[i]][,1] == CELL[j])
      flag_2 <- which(CELL_fate_diff[[i]][,1] == CELL[k])
      if(length(flag_1) >=1 & length(flag_2) >=1){
        d_1 <- CELL_fate_diff[[i]][flag_1,3]
        d_2 <- CELL_fate_diff[[i]][flag_2,3]
        p_ <- ks.test(d_1, d_2, alternative = 'two.sided', exact = F)$p.value
        pvalue <- append(pvalue, p_)
      }
    }
    CELL_fate_diff_pdist[j,k] <- mean(pvalue)
  }
}

###用时间序列ARIMA模型对原始数据和增量分别定阶
library(tseries) 
library(TSA)
library(forecast)

arima_order_diff_expressed <- vector('list',174)

for(i in 1:length(GENE_names)){
  cellname_ <- NULL
  p <- NULL
  d <- NULL
  q <- NULL
  for(j in 1:length(CELL)){
    if(Data_3[i,j] == 1){
      flag <- which(CELL_diff_expressed[[i]][,1] == CELL[j])
      if(length(flag) >= 5){
        cellname_ <- append(cellname_, CELL[j])
        x <- ts(CELL_diff_expressed[[i]][flag,3])
        ar_ <- arimaorder(auto.arima(x, max.q = 0))
        p <- append(p, ar_[1])
        d <- append(d, ar_[2])
        q <- append(q, ar_[3])
      }
      else{
        cellname_ <- append(cellname_, CELL[j])
        p <- append(p, NA)
        d <- append(d, NA)
        q <- append(q, NA)
      }
    }
  }
  arima_order_diff_expressed[[i]] <- data.frame('cell' = cellname_, 'ar_order' = p,
                                                'diff_order' = d, 'ma_order' = q)
}


arima_order_blot_expressed <- vector('list',174)

for(i in 1:length(GENE_names)){
  cellname_ <- NULL
  p <- NULL
  d <- NULL
  q <- NULL
  for(j in 1:length(CELL)){
    if(Data_3[i,j] == 1){
      flag <- which(CELL_blot_expressed[[i]][,1] == CELL[j])
      if(length(flag) >= 5){
        cellname_ <- append(cellname_, CELL[j])
        x <- ts(CELL_blot_expressed[[i]][flag,3])
        ar_ <- arimaorder(auto.arima(x, max.q = 0))
        p <- append(p, ar_[1])
        d <- append(d, ar_[2])
        q <- append(q, ar_[3])
      }
      else{
        cellname_ <- append(cellname_, CELL[j])
        p <- append(p, NA)
        d <- append(d, NA)
        q <- append(q, NA)
      }
    }

  }
  arima_order_blot_expressed[[i]] <- data.frame('cell' = cellname_, 'ar_order' = p,
                                                'diff_order' = d, 'ma_order' = q)
}




flag <- which(CELL_diff_expressed[[1]][,1] %in% CELL[59])
CELL_diff_expressed[[1]][flag,]
x <- ts(CELL_blot_expressed[[1]][flag,3])
pacf(x)
acf(x)
a <- auto.arima(x, max.q = 0, max.d = 1)
arimaorder(a)
a <- arima(x, order = c(1, 0, 0), include.mean = TRUE)
a <- arma(x, order = c(1, 0), include.intercept = TRUE)
s <- summary(a)
##统计定阶数据
##差分数据
d <- arima_order_diff_expressed
p_percent <- NULL
d_percent <- NULL
q_percent <- NULL

for(i in 1:length(GENE_names)){
  counts_p <- 0
  counts_d <- 0
  counts_q <- 0
  total <- 0
  for(j in 1:length(d[[i]][,1])){
    if(!is.na(d[[i]][j,2]) & !is.na(d[[i]][j,3]) & !is.na(d[[i]][j,4])){
      total <- total + 1
      if(d[[i]][j,2] == 0){
        counts_p <- counts_p + 1
      }
      if(d[[i]][j,3] == 0){
        counts_d <- counts_d + 1
      }
      if(d[[i]][j,4] == 0){
        counts_q <- counts_q + 1
      }
    }
  }
  p_percent <- append(p_percent,round(counts_p/total,3))
  d_percent <- append(d_percent,round(counts_d/total,3))
  q_percent <- append(q_percent,round(counts_q/total,3))
}
arima_order_percent_diff_expressed <- list('p=0_percent' = p_percent, 'd=0_percent' = d_percent,'q=0_percent' = q_percent)


##原始数据
d <- arima_order_blot_expressed
p_percent <- NULL
d_percent <- NULL
q_percent <- NULL

for(i in 1:length(GENE_names)){
  counts_p <- 0
  counts_d <- 0
  counts_q <- 0
  total <- 0
  for(j in 1:length(d[[i]][,1])){
    if(!is.na(d[[i]][j,2]) & !is.na(d[[i]][j,3]) & !is.na(d[[i]][j,4])){
      total <- total + 1
      if(d[[i]][j,2] >= 1){
        counts_p <- counts_p + 1
      }
      if(d[[i]][j,3] == 1){
        counts_d <- counts_d + 1
      }
      if(d[[i]][j,4] == 0){
        counts_q <- counts_q + 1
      }
    }
  }
  p_percent <- append(p_percent,round(counts_p/total,3))
  d_percent <- append(d_percent,round(counts_d/total,3))
  q_percent <- append(q_percent,round(counts_q/total,3))
}
arima_order_percent_blot_expressed <- list('p=0_percent' = p_percent, 'd=1_percent' = d_percent,'q=0_percent' = q_percent)

##画箱线图
library(ggplot2)
library(reshape2)
library(ggsci)

##假设检验的
r <- test_percent_blot_expressed$runtest_percent
l <- test_percent_blot_expressed$lbtest_percent
a <- test_percent_blot_expressed$adftest_percent
x <- test_percent_diff_expressed$runtest_percent
y <- test_percent_diff_expressed$lbtest_percent
z <- test_percent_diff_expressed$adftest_percent
group <- c(rep('blot',174),rep('diff blot',174))
dt <- data.frame('Runtest' = c(r,x), 'LBtest' = c(l,y), 'ADFtest' = c(a,z), 'group' = group)
df <- melt(dt, variable.name = 'test', value.name = 'pass_percent')
p <- ggplot(data=df, aes(x=test, y=pass_percent, fill = group)) + 
     geom_boxplot(alpha = 0.7) + scale_fill_lancet()+
     ylab('pass percent')
p
getwd()
ggsave(p,filename = "test_boxplot.pdf",width = 5,height = 3)


##arima定阶的箱线图

r <- arima_order_percent_blot_expressed$`p=0_percent`
l <- arima_order_percent_blot_expressed$`d=1_percent`

x <- arima_order_percent_diff_expressed$`p=0_percent`
y <- arima_order_percent_diff_expressed$`d=0_percent`

group <- c(rep('blot',174),rep('diff blot',174))
dt <- data.frame('p_is_lager_than_1___p_is_0' = c(r,x), 'd_is_1___d_is_0' = c(l,y), 'group' = group)
df <- melt(dt, variable.name = 'arima_order', value.name = 'order_percent')
p <- ggplot(data=df, aes(x=arima_order, y=order_percent, fill = group)) + geom_boxplot(alpha = 0.7) + scale_fill_lancet()
p

##AR(1)模型系数假设检验，含均值项，后面发现均值项不显著，改为不含均值项

ar_1_test_diff_expressed_no_intercept <- vector('list',174)

for(i in 1:length(GENE_names)){
  cellname_ <- NULL
  p <- NULL
  for(j in 1:length(CELL)){
    if(Data_3[i,j] == 1){
      flag <- which(CELL_diff_expressed[[i]][,1] == CELL[j])
      if(length(flag) >= 5){
        cellname_ <- append(cellname_, CELL[j])
        x <- ts(CELL_diff_expressed[[i]][flag,3])
        ar_ <- arma(x, order = c(1, 0), include.intercept = FALSE)
        ar_ <- summary(ar_)
        p <- append(p, round(ar_$coef[1,4], 3))
      }
      else{
        cellname_ <- append(cellname_, CELL[j])
        p <- append(p, NA)
      }
    }
  }
  ar_1_test_diff_expressed_no_intercept[[i]] <- data.frame('cell' = cellname_, 'ar_pvalue' = p)
                                                
}




##统计系数项、截距项显著性比例

d <- ar_1_test_diff_expressed_no_intercept
p_percent <- NULL
intercept_percent <- NULL


for(i in 1:length(GENE_names)){
  counts_p <- 0
  total <- 0
  for(j in 1:length(d[[i]][,1])){
    if(!is.na(d[[i]][j,2])){
      total <- total + 1
      if(d[[i]][j,2] >= 0.05){
        counts_p <- counts_p + 1
      }
    }
  }
  p_percent <- append(p_percent,round(counts_p/total,3))

}
ar_1_test_percent_diff_expressed_no_intercept <- list('ar_p>=0.05_percent' = p_percent)



##箱线图

r <- ar_1_test_percent_diff_expressed$`ar_p>=0.05_percent`
l <- ar_1_test_percent_diff_expressed$`intercept_p>=0.05_percent`

x <- ar_1_test_percent_diff_expressed_no_intercept$`ar_p>=0.05_percent`
y <- rep(0,174)

group <- c(rep('ar_include_intercept',174),rep('ar_not_include_intercept',174))
dt <- data.frame('ar' = c(r,x), 'intercept' = c(l,y), 'group' = group)
df <- melt(dt, variable.name = 'coefficient', value.name = 'Pass_the_significance_test_percent')
p <- ggplot(data=df, aes(x=coefficient, y=Pass_the_significance_test_percent, fill = group)) + geom_boxplot(alpha = 0.7) + scale_fill_lancet()
p
