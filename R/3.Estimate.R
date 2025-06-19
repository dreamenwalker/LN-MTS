### Part --[1]----

data.batch <- readRDS("data.batch.rds")

write.table(data.batch, file = "data.batch.txt", quote = F, sep = "\t")
GSE5851 <- "data.batch.txt"


### Part --[2]----
#' estimate

library(utils)
library(estimate)

filterCommonGenes(input.f = GSE5851, output.f = "OV_10412genes.gct", id = "GeneSymbol") 

rt <- read.table("OV_10412genes.gct", skip = 2, header = TRUE, sep = "\t")
View(rt)


estimateScore(
  input.ds = "OV_10412genes.gct", 
  output.ds = "estimateScore.gct", 
  platform = "affymetrix"
) 
scores <- read.table("estimateScore.gct", skip = 2, header = T) 
View(scores)
rownames(scores) <- scores[, 1] 
batch.estimate <- t(scores[, 3:ncol(scores)]) 
batch.estimate <- data.frame(batch.estimate)
saveRDS(batch.estimate, file = "./estimate/batch.estimate.rds")


### Part --[3]----
#' 

re <- readRDS("./3.ESTIMATE/batch.estimate.rds")
re2 <- as.data.frame(t(re))

library(openxlsx)
Clinical <- read.xlsx("Clinical.xlsx")
Pheno <- Clinical[Clinical$AnatomicSite != "Rectum", ]
Pheno <- Pheno[Pheno$SampleType == "Normal", ]
Pheno <- Pheno[Pheno$TNM.stage != "Ⅲ", ]
dim(Pheno)
table(Pheno$TNM.stage)

a <- match(Pheno$SampleID, colnames(re2))
re2 <- re2[, a] 

Pheno <- Pheno[order(Pheno$AnatomicSite), ]
metadata <- data.frame("Group" = Pheno$AnatomicSite) 
rownames(metadata) <- Pheno$SampleID

a <- match(rownames(metadata), colnames(re2))
re3 <- re2[, a]

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8, "Set1"))
dat <- data.frame(t(re3))
dat$Group <- metadata$Group 

library(tidyr)
dat <- dat %>%
  tibble::rownames_to_column("Sample") %>%
  gather(key = Estimate, value = Proportion, -c(Sample, Group)) 

library(ggplot2)
library(ggpubr)
ggplot(dat, aes(Group, Proportion)) +
  geom_violin(aes(fill = Group), scale = "width") +
  geom_boxplot(width = 0.1) +
  ggtitle("Distance") +
  scale_fill_manual(values = c("red", "#016D06", "blue")) + 
  facet_wrap(~Estimate, scales = "free_y", ncol = 5, nrow = 1) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), hide.ns = T, method = "kruskal.test") +
  theme_bw()

