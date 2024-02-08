library(snowfall)
library(parallel)
library(Matrix)
detectCores()
setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")

ptm <- proc.time()
B <- GABiclustering_initialize(N = 40)
result <- GABiclustering_fitness(B,alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                       All_gene_index = All_gene_index, All_cell_index = All_cell_index)
proc.time() - ptm

sfCpus()
sfInit(parallel = TRUE, cpus = 4)
sfCpus()
library(Matrix)

sfExportAll()
sfExport('B')
sfExport('all_corr')
sfExport('all_pvalue')


sfStop()


GABiclustering_score <- function(li = B[[1]], alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                                   All_gene_index = All_gene_index, All_cell_index = All_cell_index){
    fitness <- score(gene_name = li$`gene names`,cell_index = li$`cell index`, 
                     alpha = alpha, lambda = lambda, beta = beta, delta = delta, 
                     method = 'average', All_gene_index, All_cell_index)
    s <- fitness$'Score'
    c <- fitness$'Gene Corr Score'
    d <- fitness$'Cell Dist Score'
    a <- fitness$'Area Score'
    o <- fitness$'penalty'
    back <- c(s, c, d, a, o)
    return(back)
}

GABiclustering_score(B, alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                       All_gene_index = All_gene_index, All_cell_index = All_cell_index)

ptm <- proc.time()
result <- sfSapply(B, GABiclustering_fitness, alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                   All_gene_index = All_gene_index, All_cell_index = All_cell_index)
proc.time() - ptm

ptm <- proc.time()
result <- sapply(B, GABiclustering_fitness, alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                   All_gene_index = All_gene_index, All_cell_index = All_cell_index)
proc.time() - ptm
sfStop()


##实际上只需要修改GABiclustering_fitness的计算方法，让他用多个核跑
GABiclustering_fitness <- function(B, alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                                   All_gene_index = All_gene_index, All_cell_index = All_cell_index){

    result <- sfSapply(B, GABiclustering_score, alpha = alpha, lambda = lambda, beta = beta, delta = delta,
                         All_gene_index = All_gene_index, All_cell_index = All_cell_index)
    s <- result[1,]
    c <- result[2,]
    d <- result[3,]
    a <- result[4,]
    o <- result[5,]
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
ptm <- proc.time()
result2 <- GABiclustering_fitness(B,alpha = 1, lambda = 1,beta = 1, delta = 0.1,
                                 All_gene_index = All_gene_index, All_cell_index = All_cell_index)
proc.time() - ptm

##library(Matrix)有问题，改一下函数，把稀疏矩阵转换为稠密矩阵
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
                li <- as.matrix(li)
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


##library(Matrix)有问题，改一下函数，把稀疏矩阵转换为稠密矩阵
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
                li <- as.matrix(li)
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