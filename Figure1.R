rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing/")

# count
library(stringr)
count <- read.table("sheep.count.txt", header = T, sep = "\t", row.names = 1, check.names = F)
colnames(count)
head(count)
colnames(count) <- gsub(".bam", "", colnames(count)) #字符替换

# count to fpkm
kb <- count$Length/1000;head(kb)
rpk <- count[,-1]/kb;head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000);head(tpm)
fpkm <- t(t(rpk)/colSums(count[,-1]) * 10^6);head(fpkm)
count <- count[,-1]
head(fpkm)
fpkm <- round(fpkm, 4)
save(count, fpkm, file = "sheep.exp.raw.rda")

# load data
load("/Users/mtdeng/Documents/Research/Splicing/sheep.exp.raw.rda")
head(fpkm)
fpkmE <- fpkm[,-c(1:4, 33:42)]
head(fpkmE)

fpkmE <- fpkmE[,c(26:28, 5:18, 1:4, 29:32, 19:25)] #re-arrange the order
colnames(fpkmE) <- c("MII1", "MII2", "MII3", "E2C1", "E2C2", "E2C3", "E2C4", "E2C5", "E4C1", "E4C2", "E4C3", "E4C4", 
                     "E8C1", "E8C2", "E8C3", "E8C4", "E8C5", "E16C1", "E16C2","E16C3",  
                     "E16C4", "Morula1", "Morula2", "Morula3", "Morula4", "EB1", "EB2", "EB3", "HB1", "HB2", "HB3", "HB4")

#boxplot
library(dplyr)
library(ggpubr)
library(reshape2)
library(stringr)
explot <- fpkmE[which(rowSums(fpkmE) > 0), ]
explot <- melt(explot) #长宽数据转换
colnames(explot) <- list("gene", "sample", "fpkm")
explot$group <- case_when(str_sub(explot$sample, 1, 2) == "MI"~"MII",
                          str_sub(explot$sample, 1, 2) == "E2"~"E2C",
                          str_sub(explot$sample, 1, 2) == "E4"~"E4C",
                          str_sub(explot$sample, 1, 2) == "E8"~"E8C",
                          str_sub(explot$sample, 1, 3) == "E16"~"E16",
                          str_sub(explot$sample, 1, 1) == "M"~"M",
                          str_sub(explot$sample, 1, 2) == "EB"~"EB",
                          str_sub(explot$sample, 1, 2) == "HB"~"HB")

explot$fpkm <- log10(explot$fpkm + 0.01)
ggboxplot(explot, x = "group", y = "fpkm", fill = "group", palette = "aass", title = "", outlier.shape = 26) + 
  labs(x = "", y = bquote(Log[10]~(FPKM + 0.01))) + #ylim(-1, 2) +
  stat_compare_means(label.x = 5) + theme_bw() +
  theme(legend.position = "none",
        axis.title =  element_text(size = 10, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 9, color = "black", vjust = 1, hjust = 1, angle = 45))

fpkmE <- fpkmE[, -c(which(colnames(fpkmE) == "MII2"),
                    which(colnames(fpkmE) == "E2C5"), 
                    which(colnames(fpkmE) == "E4C3"),
                    which(colnames(fpkmE) == "E8C3"),
                    which(colnames(fpkmE) == "E16C4"),
                    which(colnames(fpkmE) == "Morula1"))]

#re-define colnames
colnames(fpkmE) <- c("MII1", "MII2", "2C1", "2C2", "2C3", "2C4", "4C1", "4C2", "4C3", 
                     "8C1", "8C2", "8C3", "8C4", "16C1", "16C2", "16C3", 
                     "Morula1", "Morula2", "Morula3", "EB1", "EB2", "EB3", "HB1", "HB2", "HB3", "HB4")

# go to boxplot in line38
fpkmE <- as.data.frame(fpkmE)
save(fpkmE, file = "sheep.filiter.fpkm.rda")

# sample cluster
library(paletteer)
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
expSet <- t(log10(fpkmE + 1))
exp_tree <- hclust(dist(expSet), method = "average")
A2Rplot(exp_tree, k = 4, boxes = F, 
        col.up = paletteer_d("RColorBrewer::Paired", n = 9), 
        col.down = paletteer_d("RColorBrewer::Paired", n = 9), 
        type = "rec", knot.pos = "mean", show.labels = T,
        lty.up = 1, lty.down = 1, lwd.up = 2, lwd.down = 2)


library(ggtree)
ggtree(ape::as.phylo(exp_tree), linetype = "solid", layout = "circular") +
  geom_nodepoint(color = "green", size = 2, alpha=.8) +
  geom_tiplab(hjust = -0.3) + 
  geom_tippoint(size = 0.1) +
  #geom_text(aes(label = node), size=2) +
  theme(legend.title = element_text(face = "bold"), 
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.text = element_text(size = rel(0.8))) +
  geom_highlight(node = 30, fill = "#FF9999") +
  geom_highlight(node = 31, fill = "#99CC00") +
  geom_highlight(node = 36, fill = "blue") 

# PCA
library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
pca.data <- fpkmE
pca.data <- pca.data[which(rowSums(pca.data) > 0), ]
pca.data <- as.data.frame(t(pca.data))
pca.data$group <- case_when(str_sub(rownames(pca.data), 1, 3) == "MII"~"MII",
                            str_sub(rownames(pca.data), 1, 2) == "2C"~"2C",
                            str_sub(rownames(pca.data), 1, 2) == "4C"~"4C",
                            str_sub(rownames(pca.data), 1, 2) == "8C"~"8C",
                            str_sub(rownames(pca.data), 1, 2) == "16"~"16C",
                            str_sub(rownames(pca.data), 1, 1) == "M"~"Morula",
                            str_sub(rownames(pca.data), 1, 2) == "EB"~"EarlyB",
                            str_sub(rownames(pca.data), 1, 2) == "HB"~"ExpB")

pca <- prcomp(pca.data[,-ncol(pca.data)], center = T, scale. = T)
df <- pca$x
df <- as.data.frame(df)
summ <- summary(pca)
# 提取主成分的方差贡献率,生成坐标轴标题
xlab <- paste0("PC1 (", round(summ$importance[2,1]*100,2), "%)")
ylab <- paste0("PC2 (", round(summ$importance[2,2]*100,2), "%)")

# generate pca using ggplot2
ggplot(data = df, aes(x = PC1, y = PC2, color = pca.data$group)) +
  geom_point(size = 3) + labs(x = xlab, y = ylab, color = "group", title = "Plot of PCA score") +
  stat_ellipse(aes(fill = pca.data$group), type = "norm", geom = "polygon", alpha = 0.2, color = NA) +
  guides(fill = "none") + theme_minimal() +
  theme(legend.position = 'right', 
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4),'cm'),
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5))

# correlation
load("sheep.filiter.fpkm.rda")
library(reshape2)
library(tidyverse)
cor <- round(cor(fpkmE), 2)
head(cor)[,1:6]

library(ggcorrplot)
ggcorrplot(cor, method = "circle", lab = T, type = "lower",
           ggtheme = theme_classic()) 

library(corrplot)
cordat <- cor.mtest(fpkmE, conf.level = 0.95)
class(cordat)

library(RColorBrewer)
display.brewer.all()

color1 <- rev(brewer.pal(6, "RdYlGn"))
color2 <- colorRampPalette(c("Blue", "White", "tomato")) 
corrplot(cor, col = color1)

cor <- round(cor(fpkmE), 2)

corrplot(cor, col = color2(10), tl.col = "black",
         method = 'pie', order = 'original', type = "upper",
         cl.cex = 0.75, tl.cex = 0.75)

corrplot.mixed(cor, lower = '', upper = 'pie', order = 'hclust', cl.cex = 0.75, tl.cex = 0.75)

COL1(sequential = c("Oranges", "Purples", "Reds", "Blues", "Greens", 
                    "Greys", "OrRd", "YlOrRd", "YlOrBr", "YlGn"), n = 200)

COL2(diverging = c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", "RdYlBu"), n = 200)

corrplot(cor, p.mat = cordat$p, 
         method = 'circle', type = 'lower', insig = 'blank', col.lim = c(-1, 1), 
         addCoef.col ='black', number.cex = 0.75, cl.cex = 1, tl.cex = 0.75, 
         order = 'original', diag = F, col = color2(10), tl.col = "black")


# heatmap
library(pheatmap)
library(reshape2)
library(dplyr)

load("sheep.filiter.fpkm.rda")
head(fpkmE)
hm <- as.data.frame(fpkmE)

#heatmap using mean fpkm in each group
{ hm$MII <- rowMeans(hm[,c(1:2)])
  hm$E2 <- rowMeans(hm[,c(3:6)])
  hm$E4 <- rowMeans(hm[,c(7:9)])
  hm$E8 <- rowMeans(hm[,c(10:13)])
  hm$E16 <- rowMeans(hm[,c(14:16)])
  hm$M <- rowMeans(hm[,c(17:19)])
  hm$EB <- rowMeans(hm[,c(20:22)])
  hm$HB <- rowMeans(hm[,c(23:26)])
  hm <- hm[, -c(1:26)]
}

head(hm)
hm <- hm[rowSums(hm) > 0, ]
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))

ph <- 
  pheatmap(hm, legend = 10, fontsize = 8, show_rownames = F, 
           scale = "row", cluster_col = F,  cluster_rows = T, angle_col = 45,
           color = c(colorRampPalette(colors = c("royalblue", "white"))(length(bk)/2), 
                     colorRampPalette(colors = c("white", "red1"))(length(bk)/2)),
           legend_breaks = seq(-2, 2, 1), breaks = bk, cutree_rows = 6)

## clusters ----
row_cluster <- cutree(ph$tree_row, k = 6)
neworder <- hm[ph$tree_row$order,]
cluster <- row_cluster[match(rownames(neworder), names(row_cluster))]
data <- data.frame(neworder, cluster)

data$Cluster <- paste0("cluster", data$cluster)
data <- data[,-(ncol(data)-1)]
data$gene <- rownames(data)
data_new <- melt(data)
head(data_new)
table(data$Cluster)

ggplot(data_new, aes(variable, log10(value + 1), group = gene)) + 
  labs( x = "", y = bquote(Log[10]~(FPKM + 1))) +
  geom_line(color = "white", size = 0.4) + 
  #geom_hline(yintercept = 1, linetype = 2) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line", size = 0.6, color = "red1") + 
  facet_wrap(Cluster~., ncol = 1) + #分面，一行多列
  theme_classic() + scale_y_continuous(breaks = seq(0, 1, 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(color = "black", size = 9),
        panel.border = element_rect(colour = "black", fill = NA, size = 0),
        axis.text.x = element_text(size = 9, color = "black", angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black", angle = 90, vjust = 0.5))

# GO analysis ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

go <- data_new[data_new$Cluster == "cluster6", ]
go <- enrichGO(gene = go$gene, 
               OrgDb = org.Hs.eg.db, 
               keyType = 'SYMBOL', 
               ont = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.05)

dotplot(go, showCategory = 20) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 27)) + 
  scale_color_continuous(low = 'royalblue', high = 'red1') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
        text = element_text(color = "black", size = 14),
        axis.text.x = element_text(size = 14, color = "black", vjust = 1, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 14, color = "black"))
write.table(go, 'cluster6.go.txt', sep = '\t', row.names = F, quote = F)


# using ClusterGVis
library(ClusterGVis)

heatmap <- as.data.frame(fpkmE)
heatmap$MII <- rowMeans(heatmap[,1:2])
heatmap$'2C' <- rowMeans(heatmap[,3:6])
heatmap$'4C' <- rowMeans(heatmap[,7:9])
heatmap$'8C' <- rowMeans(heatmap[,10:13])
heatmap$'16C' <- rowMeans(heatmap[,14:16])
heatmap$Morula <- rowMeans(heatmap[,17:19])
heatmap$EarlyB <- rowMeans(heatmap[,20:22])
heatmap$ExpB <- rowMeans(heatmap[,23:26])
heatmap <- heatmap[, -c(1:26)]

heatmap <- heatmap[rowSums(heatmap) > 0.1, ]
getClusters(exp = heatmap)

# using mfuzz for clustering
ck <- clusterData(exp = heatmap, cluster.method = "kmeans", cluster.num = 3)

# add gene name
markGenes <- c("TEAD4", "TFAP2C", "RHOA", "BTG4", "KDM5B", "YAP1")

visCluster(object = ck, plot.type = "line")

visCluster(object = ck, border = F,
           plot.type = "both",
           column_names_rot = 45,
           add.box = T,
           add.line = F,
           line.side = "left",
           boxcol = ggsci::pal_npg()(8))

# enrich for clusters
library(ggsci)
library(org.Hs.eg.db)
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 10,
                        seed = 12345)

pdf('ClusterGVis.pdf', height = 7, width = 10)
visCluster(object = ck, border = F,
           plot.type = "both",
           column_names_rot = 45,
           add.box = F,
           add.line = T,
           annoTerm.data = enrich,
           line.side = "left",
           show_row_dend = F,
           textbox.pos = c(0.5, 1.5),
           boxcol = ggsci::pal_npg()(8),
           go.col = rep(pal_d3()(10),each = 3))
dev.off()

library(reshape2)
library(ggpubr)
c2 <- ck$wide.res[ck$wide.res$cluster == "2",]
head(c2)

explot <- heatmap[rownames(heatmap) %in% rownames(c2),]
explot <- explot[rowSums(explot) > 1, ]
explot <- melt(explot) #长宽数据转换
colnames(explot) <- list("sample", "value")
explot$lg <- log10(explot$value + 1)

ggboxplot(explot, x = "sample", y = "lg", fill = "sample", palette = "aass", outlier.shape = 26) + 
  labs(x = "", y = bquote(Log[10]~(FPKM + 1))) + ylim(0, 3) +
  #stat_compare_means(label.x = 5) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.title =  element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 11, color = "black", vjust = 1, hjust = 1, angle = 45))
save(ck, enrich, file = "ClusterGVis.rda")






