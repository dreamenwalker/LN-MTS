library(dplyr)
library(openxlsx)
library(tidyr)
library(stringr)
data.batch <- readRDS("data.batch.rds")
data.batch <- data.batch %>%
  filter(if_all(everything(), ~ . >= 0))
data <- log2(data.batch+1)
dim(data)

Clinical <- read.xlsx("Clinical.xlsx")
Pheno <- Clinical

a <- match(Pheno$SampleID, colnames(data))
data <- data[, a] #
dim(data)

#' Group
Pheno$TumorSize <- as.numeric(Pheno$TumorSize)
group_list <- ifelse(Pheno$TumorSize > median(Pheno$TumorSize),'Large','Small')
group_list <- factor(group_list, levels = c("Large", "Small"))

#' Matrix

library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(data)
contrast.matrix <- makeContrasts("Large-Small", levels = design)

#' 
fit <- lmFit(data, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
allDiff <- na.omit(allDiff)

test_p <- allDiff$P.Value <= 0.05 
test_up <- allDiff$logFC >= 0.5 
test_down <- allDiff$logFC <= -0.5 
allDiff$change <- ifelse(test_p & test_up, "UP", 
                         ifelse(test_p & test_down, "DOWN", "NOT")
)
table(allDiff$change)
batch.D1 <- allDiff
saveRDS(batch.D1, file = "./1.DEG/batch.D1.rds")

library(ggplot2)
library(dplyr)

ggplot(data = batch.D1, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.6, size = 1.5, aes(color = change)) +
  ylab("-log10(FDR)") +
  scale_color_manual(values = c("#34bfb5", "#828586", "#ff6633")) +
  geom_vline(xintercept = c(-0.5,0.5), lty = 4, col = "grey", lwd = 0.8) +
  theme_classic() +ggtitle("D1 Lymph node") +
  theme(plot.title = element_text(size = 12, hjust = 0.5),legend.title = element_blank(),
  )


#' pheatmap
dat <- batch.D1[batch.D1$change != "NOT", ]
dim(dat)
exp <- data[rownames(data) %in% rownames(dat), ] #
exp <- data[rownames(data) %in% rownames(dat), ]
dim(exp)
View(exp)

col_anno <- data.frame(names = colnames(data), Group_list = group_list) 
col_anno <- col_anno[order(col_anno$Group_list), ]
rownames(col_anno) <- col_anno[, 1]
a <- match(rownames(col_anno), colnames(exp))
exp <- exp[, a]
col_anno <- data.frame(col_anno[, -1])
rownames(col_anno) <- colnames(exp)
colnames(col_anno) <- "Group"
DEG <- dat

row_anno <- data.frame(Type = DEG$change, row.names = rownames(DEG)) 

library(pheatmap)
bk <- c(seq(-4, -0.1, by = 0.01), seq(0, 4, by = 0.01))
ann_colors <- list(
  Type = c("UP" = "#ff6633", "DOWN" = "#34bfb5"),
  Group = c("Large" = "#cf6bd6", "Small" = "#99cc00")
)
pheatmap(exp,
         scale = "row", cluster_row = T, cluster_col = F, show_rownames = T, annotation_row = row_anno, annotation_col = col_anno,
         clustering_column_method = "complete", show_colnames = F, main = "D1 Lymph node", annotation_colors = ann_colors,
         color = c(colorRampPalette(colors = c("#ff6633", "white"))(length(bk) / 2), colorRampPalette(colors = c("white", "#34bfb5"))(length(bk) / 2)),
         fontsize = 8, fontsize_row = 4, legend_breaks = seq(-4, 4, 2), breaks = bk
)
View(exp)