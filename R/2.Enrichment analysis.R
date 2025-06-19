### Part -[1]-

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)

batch.deg <- readRDS("./DEG/batch.deg.rds")
dat <- batch.deg[batch.deg$change == "UP", ]
gene_id <- bitr(rownames(dat), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")


### Part -[2]-
#' GO

GO <- enrichGO(gene = gene_id[, 2],OrgDb = org.Hs.eg.db,pvalueCutoff = 1,qvalueCutoff = 1, ont = "all", readable = T) 
GO1.up <- as.data.frame(GO) 
GO1.up <- GO1.up[(GO1.up$pvalue < 0.05), ] 
write.xlsx(GO1.up, file = "./2.ENRICHMENT/GO1.up.xlsx")

#' KEGG

KEGG <- enrichKEGG(gene = gene_id[, 2],organism = "hsa", pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 1,use_internal_data = FALSE)
KEGG <- setReadable(KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
KEGG1.up <- as.data.frame(KEGG)
KEGG1.up <- KEGG1.up[(KEGG1.up$pvalue < 0.05), ]
write.xlsx(KEGG1.up, file = "./2.ENRICHMENT/KEGG1.up.xlsx")


### Part -[3]-
#' 

library(ggplot2)
GO1 <- read.xlsx('./2.ENRICHMENT/GO1.up.xlsx')
GO1$ID <- factor(GO1$ID, labels = GO1$ID) 

ggplot(GO1[1:10, ],aes(x=GeneRatio,y=Description,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("GO of D1 Lymph node")+
  xlab("GeneRatio")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(color="black",angle=0,vjust=1))+
    theme(axis.title=element_text(size=14,colour = 'black'), 
          axis.text=element_text(size=14,colour = 'black'), 
          axis.line = element_line(size=0.5, colour = 'black'), 
          panel.background = element_rect(color='black'),
          legend.key = element_blank() 
    )

library(ggplot2)
KEGG1 <- read.xlsx('./2.ENRICHMENT/KEGG1.up.xlsx')
KEGG1$ID <- factor(KEGG1$ID, labels = KEGG1$ID) 

ggplot(KEGG1[1:10, ],aes(x=GeneRatio,y=Description,colour=-1*log10(pvalue),size=Count))+
  geom_point()+
  scale_size(range=c(2, 8))+
  scale_colour_gradient(low = "blue",high = "red")+
  theme_bw()+
  ylab("KEGG of D1 Lymph node")+
  xlab("GeneRatio")+
  labs(color=expression(-log[10](PValue)))+
  theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
  theme(axis.title.y = element_text(margin = margin(r = 50)),axis.title.x = element_text(margin = margin(t = 20)))+
  theme(axis.text.x = element_text(color="black",angle=0,vjust=1))+
  theme(axis.title=element_text(size=14,colour = 'black'), 
        axis.text=element_text(size=14,colour = 'black'), 
        axis.line = element_line(size=0.5, colour = 'black'), 
        panel.background = element_rect(color='black'), 
        legend.key = element_blank() 
  )
