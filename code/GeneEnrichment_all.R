library(clusterProfiler)
library(org.Ce.eg.db)


columns(org.Ce.eg.db)
#背景基因的ID转换
entrez_id <- mapIds( x = org.Ce.eg.db,
                     keys = names_GENE,
                     keytype = "SYMBOL",
                     column = "ENTREZID")
entrez_id <- na.omit(entrez_id)
entrez_id <- data.frame(entrez_id)

###各基因转化为entrez_id，这个文件是基因间相关性的表格
setwd("/Users/yanxianzhong/Desktop/GeneClustering/code")
a <- read.csv(file="BIOGRID_C_elegans_tab3.txt",header = T,sep="\t")

save(entrez_id, file = 'entrez_id.Rdata')


##有些基因的名称不一致，建立一一对应表格，虽然NA但实际上有注释的基因
gene_name_old <- names(entrez_id)[which(is.na(entrez_id))]
gene_name_new <- c("egl-27.b",'nhr-69','swsn-7',"D1081.8",'attf-2','myrf-2','aptf-4',
                   'nfyb-1','saeg-2','lsy-27','rbm-26','sup-37',"D05D10.1",'eef-1A.1','F17C11.26',
                   'zip-8','hinf-1','gadr-1','dhhc-10',"tag-185",'madf-10')
gene_id <- mapIds( x = org.Ce.eg.db,
                   keys = gene_name_new,
                   keytype = "SYMBOL",
                   column = "ENTREZID")
gene_id_change <- data.frame(gene_name_old, gene_name_new, gene_id, row.names = NULL)
for(i in 1:21){
  num <- which(names_GENE %in% gene_name_old[i])
  entrez_id[num, 1] <- gene_id_change[i, 3]
}

which(names_GENE %in% 'egl-27.b')
entrez_id[66,] <- 175975
entrez_id[15,] <- 174121
entrez_id[29,] <- 172639

#全部基因的富集结果
columns(org.Ce.eg.db)

select(org.Ce.eg.db,'180357',c('ENTREZID','GENENAME','ENSEMBL','ONTOLOGY'))

geneEnrichmentGO <- enrichGO(gene = entrez_id[,1],
                             OrgDb = org.Ce.eg.db,
                             keyType = "ENTREZID",
                             ont = 'ALL',
                             pAdjustMethod = 'BH',
                             pvalueCutoff = 0.1,
                             qvalueCutoff = 0.2,
                             readable = T)

geneResultGO <- as.data.frame(geneEnrichmentGO@result)

##部分基因

gene_name <- B$`gene names`
gene_name <- c("F09G2.9","F39B2.1", "F16B12.6")
gene_name <- c("pha-4","B0310.2", "tbx-11",'end-3')
gene_name <- c("sma-9","skr-8", "sdc-2")
gene_name <- c("sma-9","skr-8", "sdc-2", "F16B12.6")
gene_name <- c('dpy-31','mel-28','pgp-2')
gene_name <- c("pha-4",'tbx-9')
gene_name <- c("skr-8", "sdz-38")
gene_name <- c("skr-8", "F58D2.1")
gene_name <- c("sdc-2", "F39B2.1")
gene_name <- c('nhr-79','mel-28','pgp-2')
gene_name <- c("pha-4","B0310.2", 'hlh-1')
gene_name <- c("tbx-11","B0310.2", 'tlp-1')
gene_name <- c("cnd-1","B0310.2", 'ceh-43')
gene_name <- c("F58D2.1", "F39B2.1")
gene_name <- c("tbx-38", "tbx-37", 'his-72', 'sdz-38')
gene_name <- c("tbx-38", "tbx-37", 'skr-8')



names_GENE %in% 'lin-25' 
gene_name <- All_B[[1]]$`gene names`
B <- GAbiclustering_all_result[[1]]
gene_name <- B$`gene names`
gene_name <- c("pha-4", "ceh-43","elt-6",'F16B12.6')

##直接从我建立的表格中读取就行, 注意填充顺序
##搞个函数，输入是基因名称，输出是基因转化为entrez后的id
gene_entrez <- function(gene_name){
  l <- length(gene_name)
  gene_id_cl1 <- c()
  for(i in 1:l){
    num <- which(names_GENE %in% gene_name[i])
    gene_id_cl1 <- c(gene_id_cl1, entrez_id[num, 1])
  }
  names(gene_id_cl1) <- gene_name
  gene_id_cl1 <- data.frame(gene_id_cl1)
  return(gene_id_cl1)
}

gene_name <- names_GENE[c(1,11,19,32,42,63,70,79,88,6,100,104,97,72,18,94,76)]
gene_name <- c('pgp-2', 'nhr-79', 'end-1', 'egl-27.b', 'F17C11.1', 'tbx-11', 'tps-2',
               'sdc-2', 'hmg-11', 'pha-4')
gene_name <- c('mel-28','dve-1','dpy-31','his-72','tps-2','F39B2.1')
gene_name <- c('nhr-79','tbx-11','pha-4','F16B12.6','eft-3','mel-28',
               'F17C11.1', 'F39B2.1', 'end-1', 'F09G2.9', 'end-3', 'elt-7',
               'nhr-68', 'sma-9')
gene_name <- c('nhr-79','tbx-11','pha-4','F16B12.6','eft-3','mel-28',
               'F17C11.1', 'F39B2.1', 'end-1', 'F09G2.9')
gene_name <- c('mel-28','F17C11.1','ZK185.1','nhr-68','T23H4.2','elt-7',
               'sdc-2','tbx-11','nhr-79','skr-8','end-3','his-72',
               'nhr-2','pha-4','dpy-31','pgp-2','sma-9','F39B2.1','F16B12.6')
##GO富集分析

gene_name <- All_B[[2]]$`gene names`
gene_id_cl1 <- gene_entrez(gene_name)

gene_name <- B[[1]]$`gene names`
gene_id_cl1 <- gene_entrez(gene_name)

##各个双聚类的阈值p：0.11,0.15,0.17,0.05,0.16,0.055,0.07,0.024,0.02,0.26
gene_name <- All_B[[10]]$`gene names`
gene_id_cl1 <- gene_entrez(gene_name)


gene_enrich_cl1 <- enrichGO(gene = gene_id_cl1[,1],
                            OrgDb = org.Ce.eg.db,
                            keyType = "ENTREZID",
                            ont = 'ALL',
                            universe = entrez_id[,1],
                            pAdjustMethod = "none",
                            pvalueCutoff = 0.26,
                            qvalueCutoff = 1,
                            readable = T)

gene_enrich_cl1_result <- as.data.frame(gene_enrich_cl1@result)
gene_enrich_cl1_result['GO:0008270',9]


##全部富集分析结果放一起
library(tidyverse)
gene_enrich_cl1_result %>% mutate(group = 10) -> gene_enrich_cl1_result
all_gene_enrich_result %>% mutate(group = 1) -> all_gene_enrich_result

all_gene_enrich_result <- rbind(all_gene_enrich_result, gene_enrich_cl1_result)

##把字符串分数转为数值
x <- all_gene_enrich_result$GeneRatio


b <- c()
for(i in 1:51){
  b <- c(b, a[[i]][1])
}
d <- c()
for(i in 1:51){
  d <- c(d, a[[i]][2])
}
a <- strsplit(x, '/')

b <- as.numeric(b)
d <- as.numeric(d)
x <- b/d

all_gene_enrich_result$GeneRatio <- x
##
##0.2 0.4 42五组连续4次
##0.15 0.6 67五组连续2次
#画图
dotplot(gene_enrich_cl1, color = 'pvalue', orderBy = 'x', showCategory = 6, 
        size = 'Count', x = 'GeneRatio/BgRatio', title = 'Top 6 GO Enrichment')

barplot(gene_enrich_cl1,title = 'geneEnrichmentGO')
plotGOgraph(gene_enrich_cl1)
cnetplot(gene_enrich_cl1, circular = T, colorEdge = T, size = 'Count', showCategory = 6)
heatplot(gene_enrich_cl1)

##自己画dotplot,单个双聚类的
library(tidyverse)
gene_enrich_cl1_result %>% as_tibble() -> gene_enrich_cl1_result

GeneRatio <- c(4/7, 4/7, 4/7, 3/7, 3/9)
BgRatio <- c(21/88, 22/88, 22/88, 18/88, 12/91)
FoldEnrichment <- GeneRatio/BgRatio

gene_enrich_cl1_result %>% 
  mutate(FoldEnrichment = FoldEnrichment) %>% 
  arrange(pvalue, desc = T)



library(ggplot2)

ggplot(data = gene_enrich_cl1_result, mapping = aes(x=FoldEnrichment,y=Description)) + 
  geom_point(aes(size= Count,color=-log10(pvalue),shape = factor(ONTOLOGY))) +
  scale_color_gradient(low="blue",high ="red") +
  theme(
    panel.grid.major=element_line(colour="grey"),
    panel.grid.minor=element_line(colour="grey"),
    legend.title = element_blank(),
) + theme_bw()

##自己画dotplot，全部双聚类放一起画

all_gene_enrich_result %>% as_tibble() -> all_gene_enrich_result

library(ggplot2)



ggplot(data = all_gene_enrich_result, mapping = aes(x=Description,y=group)) +
  geom_point(aes(size = GeneRatio,color=pvalue,shape = factor(ONTOLOGY))) +
  scale_color_gradient(low="blue",high ="red", trans = 'reverse') + 
  scale_y_continuous(breaks = seq(0, 10, 1)) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1,),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major=element_line(colour="grey"),
        panel.grid.minor=element_line(colour="grey"),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank()
        ) +
  labs(color = "p.value",size = "GeneRatio", shape = 'ontology', x="Term",y = "Bicluster")




##kegg富集分析

##先将entrez_id转化为kegg_id

kegg_id <- clusterProfiler::bitr_kegg(entrez_id[,1], fromType = 'ncbi-geneid', toType = 'kegg', organism = 'cel')

li <- c()
for(i in 1:104){
  m <- entrez_id[i,1]
  if(is.na(m)){
    li <- c(li, NA)
  }
  else{
    ch <- kegg_id[which(kegg_id[,1] %in% m), 2] 
    if(length(ch) > 0){
      li <- c(li, ch)
    }
    else{
      li <- c(li, NA)
    }
  }
} 

kegg_id <- data.frame('kegg_id' = li)
rownames(kegg_id) <- names_GENE

which(names_GENE %in% 'F58D2.1')
kegg_id[75,1] <- 'CELE_F58D2.1'

save(kegg_id, file = 'kegg_id.Rdata')

GENE_ID <- data.frame('name' = names_GENE, 'entrez_id' = entrez_id[,1], 'kegg_id' = kegg_id[,1])
save(GENE_ID, file = 'GENE_ID.Rdata')


##搞个函数，输入是基因名称，输出是基因转化为kegg_id
gene_kegg <- function(gene_name){
  l <- length(gene_name)
  gene_id_cl1 <- c()
  for(i in 1:l){
    num <- which(names_GENE %in% gene_name[i])
    gene_id_cl1 <- c(gene_id_cl1, kegg_id[num, 1])
  }
  names(gene_id_cl1) <- gene_name
  gene_id_cl1 <- data.frame(gene_id_cl1)
  return(gene_id_cl1)
}

gene_name <- All_B[[2]]$`gene names`
gene_id_cl1 <- gene_kegg(gene_name)

##kegg富集分析

gene_name <- All_B[[4]]$`gene names`
gene_id_cl1 <- gene_kegg(gene_name)


geneEnrichmentKEGG <- enrichKEGG(gene = gene_id_cl1[,1],
                                 organism = 'cel',
                                 keyType = "kegg",
                                 universe = kegg_id[,1],
                                 pAdjustMethod = 'BH',
                                 pvalueCutoff = 1,
                                 qvalueCutoff = 1,
                                 minGSSize = 1)

geneEnrichmentResultKEGG <- as.data.frame(geneEnrichmentKEGG@result)






dotplot(geneEnrichmentKEGG,title = 'geneEnrichmentKEGG',orderBy = "GeneRatio")
barplot(geneEnrichmentKEGG,title = 'geneEnrichmentKEGG')












#KEGG分析第十一类细胞基因聚类（没有通路）（不过全部基因存在通路）

uniprot_id <- clusterProfiler::bitr(geneID = names_GENE, fromType = "SYMBOL",toType = "UNIPROT", OrgDb = 'org.Ce.eg.db')
kegg_id <- clusterProfiler::bitr_kegg(geneID = uniprot_id[,2], fromType = "uniprot",toType = "kegg", organism = 'cel')

gene_name <- All_B[[2]]$`gene names`

geneIdUniprot <- clusterProfiler::bitr(geneID = gene_name, fromType = "SYMBOL",toType = "UNIPROT", OrgDb = 'org.Ce.eg.db')
geneIdkegg <-  clusterProfiler::bitr_kegg(geneID = geneIdUniprot[,2], fromType = "uniprot",toType = "kegg", organism = 'cel')







