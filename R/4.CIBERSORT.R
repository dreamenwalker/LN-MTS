### Part --[1]----
#' 

library(tibble)
data.batch <- readRDS("data.batch.rds") 
mixture.file = tibble::rownames_to_column(data.batch, "Gene_symbol")
write.table(mixture.file, file = "mixture.file.txt", row.names = F, quote = F, sep = "\t") 


### Part --[2]----
#' CIBERSORT

source("CIBERSORT.R") 
TME.results <- CIBERSORT("mixture.file.txt","LM22.txt",perm = 500, QN = T) 
batch.cibersort <- TME.results
batch.cibersort <- batch.cibersort[, -(23:25)] 
batch.cibersort <- data.frame(batch.cibersort)
saveRDS(batch.cibersort, file = "./CIBERSORT/batch.cibersort.rds")


### Part --[3]----
#' 
re <- readRDS("./4.CIBERSORT/batch.cibersort.rds")
re2 <- as.data.frame(t(re))

library(openxlsx)
Clinical <- read.xlsx("Clinical.xlsx")
Pheno <- Clinical[Clinical$AnatomicSite != "Rectum", ]
Pheno <- Pheno[Pheno$SampleType == "Normal", ]
Pheno <- Pheno[Pheno$TNM.stage != "Ⅲ", ]
dim(Pheno)

a <- match(Pheno$SampleID, colnames(re2))
re2 <- re2[, a] 

Pheno <- Pheno[order(Pheno$AnatomicSite), ] 
metadata <- data.frame("Group" = Pheno$AnatomicSite) 
rownames(metadata) <- Pheno$SampleID

a <- match(rownames(metadata), colnames(re2))
re3 <- re2[, a]

groupcolor <- c("red", "#016D06")
names(groupcolor) <- c("D1 Lymph node", "D2 Lymph node") 

ann_colors <- list(Group = groupcolor) 

library(pheatmap)
pheatmap(re3,
  scale = "column",
  show_colnames = F, cluster_rows = T, cluster_cols = F,
  annotation_col = metadata, main = "Pheatmap of Immune Cells",
  annotation_colors = ann_colors,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(500)
)


### Part -[4]-
#' 

library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
mypalette <- colorRampPalette(brewer.pal(8, "Set1"))
dat <- data.frame(t(re3))
dat$Group <- metadata$Group 

dat <- dat %>%
  tibble::rownames_to_column("Sample") %>%
  gather(key = Cell_type, value = Proportion, -c(Sample, Group)) 

library(stringr)
library(ggpubr)
ggplot(dat, aes(Cell_type, Proportion, fill = Group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  ggtitle("Distance") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6)) +
  scale_fill_manual(values = c("#AA51A2", "#C0C000")) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), hide.ns = T, method = "kruskal.test")

