#画有表达的细胞的blot值散点图
library(RColorBrewer)
display.brewer.pal(n = 9, name = 'Set1')
colorName <- brewer.pal(n = 4, name = 'RdBu')  
colorName <- brewer.pal(n = 9, name = 'Set1')[c(2,3,5,8)]
barplot(rep(1,6),1, col = colorName)
#k，第几个基因文件
#m，起始细胞
#n，终止细胞
#x，x轴起点
#y，x轴终点
DrawExpressedBlot <- function(k=1,m,n,x=20,y=140){
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(pha4_1_onset_realdata[[i]][[3]])
    pha4_1_onset_minblot[i-n+1] <- min(pha4_1_onset_realdata[[i]][[3]])
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'blot')
  for(i in m:n){
    d <- pha4_1_onset_realdata[[i]]$blot
    points(pha4_1_onset_realdata[[i]]$time,d,pch = 4,col = colorName[i-m+1])
    
  }
  
  legend('topleft',pch = 4,col = colorName, ncol=2, adj = 0, text.width = 2, y.intersp = 2,legend = pha4_1_onset$cell[m:n], cex = 1, bty = 'n')
  
}
# ##
# DrawExpressedBlot(30,9,12,65,140)
# 
# DrawExpressedBlot(38,15,18,50,100)

DrawExpressedBlot(30,9,12,65,140)

DrawExpressedBlot(38,15,18,50,100)


#画有表达的细胞的blot值的增量散点图(正增量)

DrawExpressedDiffBlot <- function(k=1,m,n,x=20,y=140){
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(diff(pha4_1_onset_realdata[[i]][[3]]))
    pha4_1_onset_minblot[i-n+1] <- min(diff(pha4_1_onset_realdata[[i]][[3]]))
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'diff(blot)')
  for(i in m:n){
    d <- diff(pha4_1_onset_realdata[[i]]$blot)
    dtime <- which(d>=0)
    d <- d[d>=0]
    points(pha4_1_onset_realdata[[i]]$time[-1][dtime+1],d,pch = 4,col = colorName[i-m+1])
  }
  
  legend('topleft',pch = 4,col = colorName, ncol=1, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}

DrawExpressedDiffBlot(1,1,8,130,200)



##画折线图
DrawExpressedBlot_line <- function(k=1,m,n,x=20,y=140,cell_set){
  cat('负增量的时间点占比','\n')
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(pha4_1_onset_realdata[[i]][[3]])
    pha4_1_onset_minblot[i-n+1] <- min(pha4_1_onset_realdata[[i]][[3]])
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'blot')
  for(i in m:n){
    
    for(val in cell_set){
      flag <- which(pha4_1_onset_realdata[[i]]$cell == val)
      if(length(flag) > 0){
        d <- diff(pha4_1_onset_realdata[[i]]$blot[flag])
        dtime <- which(d>=0)
        lines(pha4_1_onset_realdata[[i]]$time[flag][dtime+1],pha4_1_onset_realdata[[i]]$blot[flag][dtime+1],col = colorName[i-m+1])
        cat(val,':',length(dtime)/length(d),'\n')
      }
    }
  }
  legend('topleft',pch = 4,col = colorName, ncol=4, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}

DrawExpressedBlot_line(1,1,8,130,200,c(pha4_1_onset[1:8,1],CELL[33:48]))



#画增量折线图(正增量)
DrawExpressedDiffBlot_line <- function(k=1,m,n,x=20,y=140,cell_set){
  cat('负增量的时间点占比','\n')
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(pha4_1_onset_realdata[[i]][[3]])
    pha4_1_onset_minblot[i-n+1] <- min(pha4_1_onset_realdata[[i]][[3]])
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(-20000,20000),type = 'n',xlab = 'time',ylab = 'blot')
  for(i in m:n){
    
    for(val in cell_set){
      flag <- which(pha4_1_onset_realdata[[i]]$cell == val)
      if(length(flag) > 1){
        d <- diff(pha4_1_onset_realdata[[i]]$blot[flag])
        dtime <- which(d>=0)
        lines(pha4_1_onset_realdata[[i]]$time[flag][-1][dtime],d[dtime],col = colorName[i-m+1])
        cat(val,':',length(dtime)/length(d),'\n')
      }
    }
  }
  legend('topleft',pch = 4,col = colorName, ncol=4, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}

DrawExpressedDiffBlot_line(1,1,8,130,200,c(pha4_1_onset[1:8,1],CELL[33:48]))




#画增量折线图(正，负增量)
DrawExpressedDiffBlot_line <- function(k=1,m,n,x=20,y=140,cell_set){
  cat('负增量的时间点占比','\n')
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(pha4_1_onset_realdata[[i]][[3]])
    pha4_1_onset_minblot[i-n+1] <- min(pha4_1_onset_realdata[[i]][[3]])
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(-20000,20000),type = 'n',xlab = 'time',ylab = 'blot')
  for(i in m:n){
    
    for(val in cell_set){
      flag <- which(pha4_1_onset_realdata[[i]]$cell == val)
      if(length(flag) > 1){
        d <- diff(pha4_1_onset_realdata[[i]]$blot[flag])
        # dtime <- which(d>=0)
        lines(pha4_1_onset_realdata[[i]]$time[flag][-1],d,col = colorName[i-m+1])
        cat(val,':',length(dtime)/length(d),'\n')
      }
    }
  }
  legend('topleft',pch = 4,col = colorName, ncol=4, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}

DrawExpressedDiffBlot_line(1,1,8,130,200,c(pha4_1_onset[1:8,1],CELL[33:48]))



#画有表达的细胞的blot值的增量散点图(正,负增量)

DrawExpressedDiffBlot <- function(k=1,m,n,x=20,y=140){
  dif <- NULL
  pha4_1_onset <- hujiedata[[which(file_newnames%in%file_gene[k,1])]]
  pha4_1_realdata <- realdata[[which(file_names_csv%in%file_gene[k,1])]]
  pha4_1_onset_realdata <- vector('list',length = dim(pha4_1_onset)[1])
  for(i in 1:dim(pha4_1_onset)[1]){
    pha4_1_onset_realdata[[i]] <- subset(pha4_1_realdata,pha4_1_realdata$cell %in% grep(pha4_1_onset[i,1],cell,value = T) & pha4_1_realdata$time >= pha4_1_onset[i,2])
    
  }
  pha4_1_onset_maxblot <- 0
  pha4_1_onset_minblot <- 0
  for(i in m:n){
    pha4_1_onset_maxblot[i-m+1] <- max(diff(pha4_1_onset_realdata[[i]][[3]]))
    pha4_1_onset_minblot[i-n+1] <- min(diff(pha4_1_onset_realdata[[i]][[3]]))
    
  }
  pha4_1_onset_maxblot <- max(pha4_1_onset_maxblot)
  pha4_1_onset_minblot <- min(pha4_1_onset_minblot)
  plot(1,1,xlim = c(x,y),ylim = c(pha4_1_onset_minblot,pha4_1_onset_maxblot),type = 'n',xlab = 'time',ylab = 'diff(blot)')
  for(i in m:n){
    d <- diff(pha4_1_onset_realdata[[i]]$blot)
    # dtime <- which(d>=0)
    points(pha4_1_onset_realdata[[i]]$time[-1],d,pch = 4,col = colorName[i-m+1])
    dif <- append(dif,d) 
  }
  
  legend('topleft',pch = 4,col = colorName, ncol=1, adj = 1, text.width = 2, y.intersp = 0.3,legend = pha4_1_onset$cell[m:n], cex = 0.8, bty = 'n')
  
}

DrawExpressedDiffBlot(1,1,8,130,200)




