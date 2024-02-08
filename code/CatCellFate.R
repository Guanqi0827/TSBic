
##整理145个细胞命运

cell_index <- B[[3]]
Cell_names <- B[[4]]
Cell_names %in% CELL_fate_names
CELL_fate_names %in% cell_fate_original_hu[,1]
CELL_fate_original_hu <- data.frame(nrow = 145,ncol = 2)
which(CELL_fate_names %in% cell_fate_original_hu[,1])
for(i in 1:145){
  index <- which(cell_fate_original_hu[,1] %in% CELL_fate_names[i])
  CELL_fate_original_hu[i,1] <- CELL_fate_names[i]
  CELL_fate_original_hu[i,2] <- cell_fate_original_hu[index,2]
  
}

##看双聚类的细胞命运
CatCellFate <- function(cell_index = c(1:24)){
  library(RColorBrewer)
  Cell_names <- CELL[cell_index]
  CELL_fate <- NULL
  for (i in 1:length(cell_index)){
    index <- which(CELL_fate_original_hu[,1] == Cell_names[i])
    CELL_fate[i] <- CELL_fate_original_hu[,2][index]
  }
  if(length(table(CELL_fate))>=3){
    colorName <- brewer.pal(n = length(table(CELL_fate)), name = 'RdBu')
  }
  else{
    colorName <- brewer.pal(n = 4, name = 'RdBu')  
  }
  labelName <- names(table(CELL_fate))
  percent <- round(100*table(CELL_fate)/length(cell_index),2)
  percent <- paste(percent, "%", sep = "")
  pie(table(CELL_fate), labels = labelName, main = paste('number of cells:',length(cell_index)),cex = 0.8,radius = 1,col = colorName)
  legend("bottomleft",percent, cex = 0.5,fill = colorName)
  back <- list('Cell name' <- Cell_names, 'Cell fate' = CELL_fate,
               'all fate' = table(CELL_fate), 'Percent' = percent)
  return(back)

}
##325个细胞的，为了画细胞聚类的那个饼状图

cellFateNum <- table(cell_fate_original_hu[,2])

cell_index <- which(cellFateOriginalName %in% cell_class[[7]])
CatCellFate(cell_index)


cellEnrichment(allcell = cellFateOriginalName, cellclass = cell_class[[7]])

##这个函数是画细胞聚类饼图的
CatCellFate <- function(Cell_names = "MSpppaa"){

  CELL_fate <- NULL
  colorName <- c('#FFF1C9', '#F7B7A3', '#ea5f90', '#9b3192', '#57167e', '#2b0b3f')
  
  for(i in 1:length(Cell_names)){
    index <- which(cellFateOriginalName == Cell_names[i])
    CELL_fate[i] <- cell_fate_original_hu[,2][index]
  }
  Tab <- table(CELL_fate)
  labelName <- NULL
  cell_enrichment_label <- cellEnrichment(allcell = cellFateOriginalName, cellclass = Cell_names)
  nam <- names(cell_enrichment_label)
  percent <- NULL
  tab_new <- c()
  for(k in 1:length(nam)){
    a <- Tab[which(names(Tab) %in% nam[k])]
    b <- cellFateNum[names(cellFateNum) == names(a)]
    c <- round(as.numeric(cell_enrichment_label[k]),3)
    labelName <- c(labelName,paste(names(cell_enrichment_label[k]),':',as.numeric(a),
                                   '/',as.numeric(b), '(',c,')'))
    perc <- round(100*as.numeric(a)/length(Cell_names),2)
    perc <- paste(perc, "%", sep = "")
    percent <- c(percent, perc)
    tab_new <- c(tab_new, as.numeric(a))
  }
  names(tab_new) <- nam  
  pie(tab_new, labels = percent, main = paste('Cell Number:',length(Cell_names)),
      cex = 0.8,radius = 1,col = colorName)
  legend('bottomleft',labelName, cex = 0.6, fill = colorName, bty = 'n')
  back <- list('Cell name' <- Cell_names, 'Cell fate' = CELL_fate,
               'all fate' = table(CELL_fate), 'Percent' = percent)
  return(back)
  
}

CatCellFate(cell_class[[7]])



##看双聚类的细胞命运,改进一下，做细胞富集分析
##重新更新一下fate的数据，现在只有145个细胞
options(digits = 7)
cell_fate_original_hu <- cell_fate_original[cellFateHuIndex,1:2]
tab <- NULL
for(i in 1:length(CELL)){
  index <- which(cell_fate_original_hu[,1] == CELL[i])
  tab <- c(tab, cell_fate_original_hu[index,2])
}

cellFateNum <- table(tab)


###重写一下，希望p值按从小到大排序

cellEnrichment <- function(allcell, cellclass){
  cellnum <- length(allcell)
  classcellnum <- length(cellclass)
  cellClassFate <- NULL
  cellClassFateIndex <- NULL
  for(j in 1:length(cellclass)){
    index <- which(cell_fate_original_hu[,1] == cellclass[[j]])
    cellClassFateIndex <- append(cellClassFateIndex,index)
    cellClassFate <- append(cellClassFate,cell_fate_original_hu[index,2])
    
  }
  pvalue <- NULL
  for(i in 1:length(table(cellClassFate))){
    classfatenum <- as.numeric(table(cellClassFate)[i])
    fatenum <- as.numeric(cellFateNum[names(table(cellClassFate))[i]])
    pvalue <- append(pvalue,1-phyper(classfatenum-1, fatenum, cellnum-fatenum,classcellnum))
  }
  pvalueadj <- p.adjust(pvalue,method = 'BH',n = length(pvalue))
  names(pvalueadj) <- names(table(cellClassFate))
  pvalueadj <- sort(pvalueadj)
  return(pvalueadj)
}

cellEnrichment(allcell = cellFateOriginalName, cellclass = cell_class[[7]])



##这里重新创建table，为了按让p值从小到大来画饼图
library(RColorBrewer)
display.brewer.pal(n = 11, name = 'PiYG')
display.brewer.pal(n = 11, name = 'RdBu')

colorName <- c('#FFF1C9', '#F7B7A3', '#ea5f90', '#9b3192', '#57167e', '#2b0b3f')

colorName <- c('#D47a72', '#eac8c2', '#535439', '#f0e5dc', '#563f40')
colorName <- c('#fb9489', '#a9ddd4', '#9ec3db', '#cbc7de', '#fdfcc9')
colorName <- c('#58539f', '#bbbbd6', '#d86967', '#eebabb')

colorName <- c('#cea9bc', '#323232', '#2085ec', '#72b4eb', '#0a417a')
colorName <- c('#33539e', '#7facd6', '#bfb8da', '#e8b7d4', '#a5678e')
colorName <- c('#a65f87', '#bc74a4', '#51264c', '#f6a479', '#fa943c', '#613260')

barplot(rep(1,6), 1, col = colorName)
CatCellFate <- function(cell_index = c(1:24)){

  Cell_names <- CELL[cell_index]
  CELL_fate <- NULL
  colorName <- c('#FFF1C9', '#F7B7A3', '#ea5f90', '#9b3192', '#57167e', '#2b0b3f')
 
  for(i in 1:length(cell_index)){
    index <- which(CELL_fate_original_hu[,1] == Cell_names[i])
    CELL_fate[i] <- CELL_fate_original_hu[,2][index]
  }
  Tab <- table(CELL_fate)
  labelName <- NULL
  cell_enrichment_label <- cellEnrichment(allcell = CELL, cellclass = CELL[cell_index])
  nam <- names(cell_enrichment_label)
  percent <- NULL
  tab_new <- c()
  for(k in 1:length(nam)){
    a <- Tab[which(names(Tab) %in% nam[k])]
    b <- cellFateNum[names(cellFateNum) == names(a)]
    c <- round(as.numeric(cell_enrichment_label[k]),3)
    labelName <- c(labelName,paste(names(cell_enrichment_label[k]),':',as.numeric(a),
                                   '/',as.numeric(b), '(',c,')'))
    perc <- round(100*as.numeric(a)/length(cell_index),2)
    perc <- paste(perc, "%", sep = "")
    percent <- c(percent, perc)
    tab_new <- c(tab_new, as.numeric(a))
  }
  names(tab_new) <- nam  
  pie(tab_new, labels = percent, main = paste('Cell Number:',length(cell_index)),
      cex = 0.8,radius = 1,col = colorName, clockwise = T)
  legend('topright',labelName, cex = 0.6, fill = colorName, bty = 'n')
  back <- list('Cell name' <- Cell_names, 'Cell fate' = CELL_fate,
               'all fate' = table(CELL_fate), 'Percent' = percent)
  return(back)
  
}




###空心圆环图

##颜色
##  "#F07971","#CC6969","#D69037","#ACA33E","#64B24B","#30B35E","#499CD4"

rgb(78,159,215,maxColorValue = 255)
colorName <- c("#C69339","#CC6969","#D69037","#ACA33E","#64B24B","#F07971","#30B35E","#9B85BB","#A7CEE2","#14B7DB","#499CD4")
colorName <- c("#64B24B","#F07971","#30B35E","#9B85BB","#14B7DB","#4E9FD7")

barplot(1:10,1, col = colorName)


library(ggplot2)
library(ggforce)

###新的画图函数, 画细胞聚类的时候直接输入细胞名称，
##然后把cellFateNum换过来，先跑上面那个细胞聚类函数，得到df，再用下面的ggplot
Cell_names <- cell_class[[7]]
cell_index <- All_B[[2]]$`cell index`

CatCellFate <- function(cell_index = c(1:24)){
  library(ggplot2)
  library(ggforce)
  Cell_names <- CELL[cell_index]
  CELL_fate <- NULL
  colorName <- c("#F07971","#30B35E","#14B7DB","#9B85BB","#C69339","#4E9FD7")
  for(i in 1:length(cell_index)){
    index <- which(CELL_fate_original_hu[,1] == Cell_names[i])
    CELL_fate[i] <- CELL_fate_original_hu[,2][index]
  }
  Tab <- table(CELL_fate)
  labelName <- NULL
  cell_enrichment_label <- cellEnrichment(allcell = CELL, cellclass = CELL[cell_index])
  nam <- names(cell_enrichment_label)
  percent <- NULL
  tab_new <- c()
  for(k in 1:length(nam)){
    a <- Tab[which(names(Tab) %in% nam[k])]
    b <- cellFateNum[names(cellFateNum) == names(a)]
    c <- round(as.numeric(cell_enrichment_label[k]),3)
    labelName <- c(labelName,paste(names(cell_enrichment_label[k]),':',as.numeric(a),
                                   '/',as.numeric(b), '(',c,')'))
    perc <- round(100*as.numeric(a)/length(cell_index),2)
    perc <- paste(perc, "%", sep = "")
    percent <- c(percent, perc)
    tab_new <- c(tab_new, as.numeric(a))
  }
  names(tab_new) <- nam  
  df <- data.frame('label' = percent, 'num' = as.numeric(tab_new), 'fate' = labelName)
  df$fate <- factor(df$fate,levels = df$fate)
  ggplot(df, aes(x = 1.5, y = num, fill = fate)) +
    geom_col(color = "white")+
    scale_fill_manual(values = colorName)+
    coord_polar(theta = "y", direction = -1) +
    xlim(c(0, 2)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(color = 'black', size = 12, hjust = 0.5, vjust = -4,face="bold")
    )+
    ggtitle(paste(length(cell_index),'cells',sep=' '))
}

CatCellFate(All_B[[2]]$`cell index`)

barplot(rep(1,6),1, col = colorName)


