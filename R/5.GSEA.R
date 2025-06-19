### Part -[1]-
#' 

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

batch.deg <- readRDS("./1.DEG/batch.large.small.deg.rds")
batch.deg <- tibble::rownames_to_column(batch.deg,"Symbol")
genelist <- batch.deg[,1:2] 
colnames(genelist) <- c("SYMBOL","logFC")

#' 

gene_id <- bitr(genelist$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_id <- merge(gene_id,genelist,by = 'SYMBOL')
gene_id <- gene_id[order(gene_id$logFC,decreasing = T),] 
gene_id <- distinct(gene_id,ENTREZID,.keep_all = T)

discovery_gene <- gene_id[,-c(1,2)]
names(discovery_gene) <- gene_id$ENTREZID


### Part -[2]-

library(DOSE)

kegg <- gseKEGG(discovery_gene,organism = 'hsa',pvalueCutoff = 0.05,pAdjustMethod = "none" )
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")  

discovery_kegg <- kegg
save(discovery_kegg,file = "./5.GSEA/discovery_large.small.RData")


### Part -[3]-
#'
load("./5.GSEA/discovery_large.small.RData")

gsea <-discovery_kegg
data <- gsea@result

paths <- c("hsa04658","hsa03030","hsa03410","hsa04659")
gseaplot2(gsea, paths, 
          color = "red", pvalue_table = F,title = "")

