rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing/")

# DEGs-----
library(DESeq2)
load("sheep.filiter.count.rda")
load("sheep.filiter.fpkm.rda")
group <- factor(c(rep('MII', 2), rep('E2C', 4), rep('E4C', 3), rep('E8C', 4), 
                  rep('E16C', 3), rep('Morula', 3), rep('EarlyB', 3), rep('HB', 4)), 
                levels = c('MII', 'E2C', 'E4C', 'E8C', 'E16C', 'Morula', 'EarlyB', 'HB'))
colData <- data.frame(row.names = colnames(countE), group)
colData

dds <- DESeqDataSetFromMatrix(countE, colData, design = ~ group)
dds <- DESeq(dds)

E16vsMII <- results(dds, contrast = c("group", "E16C", "MII"))
table(E16vsMII$padj < 0.05)
E16vsMII <- cbind(E16vsMII, fpkmE)
E16vsMII <- as.data.frame(E16vsMII)
E16vsMII$change <- ifelse(E16vsMII$padj > 0.05, 'stable',
                          ifelse(E16vsMII$log2FoldChange > 1, 'up', 
                                 ifelse(E16vsMII$log2FoldChange < -1, 'down', 'stable')))
table(E16vsMII$change)
DEGs.E16vsMII <- subset(E16vsMII, E16vsMII$change != "stable")

MvsMII <- results(dds, contrast = c("group", "Morula", "MII"))
table(MvsMII$padj < 0.05)
MvsMII <- cbind(MvsMII, fpkmE)
MvsMII <- as.data.frame(MvsMII)
MvsMII$change <- ifelse(MvsMII$padj > 0.05, 'stable',
                        ifelse(MvsMII$log2FoldChange > 1, 'up', 
                               ifelse(MvsMII$log2FoldChange < -1, 'down', 'stable')))

table(MvsMII$change)
DEGs.MvsMII <- subset(MvsMII, MvsMII$change != "stable")

table(rownames(DEGs.MvsMII) %in% rownames(DEGs.E16vsMII))
table(rownames(DEGs.E16vsMII) %in% rownames(DEGs.MvsMII))
save(MvsMII, E16vsMII, file = "sheep.degs.rda")

# sheep maternal genes
library(pheatmap)
load("sheep.degs.rda")
sheep.degs <- subset(MvsMII, MvsMII$change != "stable")
mat <- subset(sheep.degs, sheep.degs$change == "down")
table(mat$change)

#dataframe using mean fpkm
{
  hm <- mat[, -c(1:6, ncol(mat))]
  hm <- hm[rowSums(hm) > 0, ]
  hm$MII <- rowMeans(hm[,c(1:2)])
  hm$E2 <- rowMeans(hm[,c(3:6)])
  hm$E4 <- rowMeans(hm[,c(7:9)])
  hm$E8 <- rowMeans(hm[,c(10:13)])
  hm$E16 <- rowMeans(hm[,c(14:16)])
  hm$M <- rowMeans(hm[,c(17:19)])
  hm$EB <- rowMeans(hm[,c(20:22)])
  hm$HB <- rowMeans(hm[,c(23:26)])
  hm <- hm[,-c(1:26)]
  head(hm)
  }

bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
ph <- pheatmap(hm, legend = 12, fontsize = 12, show_rownames = F, 
               scale = "row", cluster_col = F,  cluster_rows = T, angle_col = 45,
               color = c(colorRampPalette(colors = c("blue", "white"))(length(bk)/2), 
                         colorRampPalette(colors = c("white", "red1"))(length(bk)/2)),
               legend_breaks = seq(-2, 2, 1), breaks = bk, cutree_rows = 2)


## clusters ----
library(ggplot2)
library(reshape2)

row_cluster <- cutree(ph$tree_row, k = 2)
neworder <- hm[ph$tree_row$order,]
cluster <- row_cluster[match(rownames(neworder), names(row_cluster))]
data <- data.frame(neworder, cluster)
data$cluster <- paste0("cluster", data$cluster)
data <- data[,-(ncol(data)-1)]
data$gene <- rownames(data)
table(data$cluster)

data_new <- melt(data)
head(data_new)
ggplot(data_new, aes(variable, log2(value + 1), group = gene)) + 
  labs( x = "", y = bquote(Log[2]~(FPKM + 1))) +
  geom_line(color = "white", size = 0.6) + 
  geom_hline(yintercept = 1, linetype = 2) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line", size = 0.8, color = "red1") + 
  facet_wrap(cluster~., nrow = 3, ncol = 1) + #分面，一行多列
  theme_classic() + scale_y_continuous(breaks = seq(0, 4, 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(color = "black", size = 7),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
        axis.text.x = element_text(size = 7, color = "black", angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 7, color = "black", angle = 90, vjust = 0.5))


# GO analysis ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

go <- enrichGO(gene = rownames(mat), 
               OrgDb = org.Hs.eg.db, 
               keyType = 'SYMBOL', 
               ont = "BP",
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

##GSEA-----
gsea <- mat[mat$baseMean != 0,]
ENTREID <- bitr(rownames(gsea), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
rownames(ENTREID) <- ENTREID$SYMBOL

gesa.dat <- merge(ENTREID, gsea, by = 0)
gesa.sort <- gesa.dat[order(gesa.dat$log2FoldChange, decreasing = T),]
rownames(gesa.sort) <- gesa.sort$ENTREZID
gene_all <- gesa.sort$log2FoldChange
names(gene_all) <- gesa.sort$ENTREZID

library(enrichplot)
gseGO_all <- gseGO(gene_all, OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID", ont = "all", by = "fgsea",
                   minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 0.05)

gseGO_all@result[["Description"]]
gseaplot(gseGO_all, 9, color = "red", color.line = "red", pvalue_table = F)
#gsea term
#9, 25, 46, 47, 85
gseaplot2(gseGO_all, c(47), color = c("royalblue","steelblue"), rel_heights = c(0.6, 0.1, 0.3), 
          pvalue_table = F, ES_geom = "line", base_size = 7, title = "")


# specific genes expression -----
# TEAD4, TFAP2C, RHOA, BTG4, KDM5B, YAP1

library(ggsci)
library(ggpubr)
which(rownames(mat) == "KDM5B")
gene <- "KDM5B"
plot <- subset(fpkmE, rownames(fpkmE) == gene)
plot <- melt(plot)
plot$group <- factor(t(group))
head(plot)
plot$log <- log2(plot$value + 1)

ggviolin(plot, x = "group", y = "log", fill = "group", palette = "npg", size = 0) + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + theme_bw() + 
  theme(legend.position = "none", 
        panel.background = element_rect(fill = NA, colour = " black", size = 0.8),
        axis.title =  element_text(size = 9, color = "black"), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color = "black", vjust = 1, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 8, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  annotate("text", x = 4.5, y = Inf, label = gene, vjust = 1.5, size = 3)

#####################################WGCNA######################################
#using all genes instead of DEGs
#数据进行variance-stabilizing transformation标准化，也可以使用log2(FPKM+1)

library(WGCNA)
library(dplyr)
allowWGCNAThreads()  #打开多线程

load("/Users/mtdeng/Documents/Research/Splicing/ClusterGVis.rda")
load("/Users/mtdeng/Documents/Research/Splicing/sheep.exp.raw.rda")
head(fpkm)
fpkmE <- fpkm[,-c(1:4, 33:42)]
head(fpkmE)

fpkmE <- fpkmE[,c(26:28, 5:18, 1:4, 29:32, 19:25)] #re-arrange the order
colnames(fpkmE) <- c("MII1", "MII2", "MII3", "2C1", "2C2", "2C3", "2C4", "2C5", 
                     "4C1", "4C2", "4C3", "4C4", "8C1", "8C2", "8C3", "8C4", "8C5", 
                     "16C1", "16C2","16C3", "16C4", "Morula1", "Morula2", "Morula3", "Morula4", 
                     "EB1", "EB2", "EB3", "HB1", "HB2", "HB3", "HB4")

fpkmE <- fpkmE[, -c(which(colnames(fpkmE) == "MII2"),
                    which(colnames(fpkmE) == "2C5"), 
                    which(colnames(fpkmE) == "4C3"),
                    which(colnames(fpkmE) == "8C3"),
                    which(colnames(fpkmE) == "16C4"),
                    which(colnames(fpkmE) == "Morula1"))]

head(ck$wide.res)

library(ggpubr)
library(reshape2)
library(stringr)

c1 <- ck$wide.res[ck$wide.res$cluster == "1", ]

explot <- fpkmE[rownames(fpkmE) %in% rownames(c1),]
#explot <- explot[rowSums(explot) > 1, ]
explot <- melt(explot) #长宽数据转换
colnames(explot) <- list("gene", "sample", "value")
head(explot)

explot$group <- case_when(str_sub(explot$sample, 1, 3) == "MII"~"MII",
                          str_sub(explot$sample, 1, 2) == "2C"~"2C",
                          str_sub(explot$sample, 1, 2) == "4C"~"4C",
                          str_sub(explot$sample, 1, 2) == "8C"~"8C",
                          str_sub(explot$sample, 1, 2) == "16"~"16C",
                          str_sub(explot$sample, 1, 3) == "Mor"~"Morula",
                          str_sub(explot$sample, 1, 2) == "EB"~"EarlyB",
                          str_sub(explot$sample, 1, 2) == "HB"~"ExpB")

explot$lg <- log2(explot$value + 1)

ggboxplot(explot, 
          x = "group", y = "lg", fill = "group", 
          palette = "aass", notch = T, outlier.shape = 26) + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + 
  theme_bw() + ylim(0, 12) +
  theme(legend.position = "none",
        axis.title =  element_text(size = 11, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 10, color = "black", vjust = 1, hjust = 1, angle = 45))


#WGCNA
RNAexpr <- fpkmE[rownames(c1) %in% rownames(fpkmE), ]
RNAexpr <- RNAexpr[rowSums(RNAexpr) > 1,]
RNAexpr <- log2(RNAexpr + 1)
dim(RNAexpr);head(RNAexpr)

m.mad <- apply(RNAexpr, 1, mad)
datExpr <- RNAexpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.05))[2],0.01)),]
datExpr <- as.data.frame(t(datExpr)) #转换成列为基因名
head(datExpr)[, 1:6]

#样品聚类
par(mar = c(1,2,1,1), mfrow = c(1, 1), cex = 0.5)
datExpr.tree <- hclust(dist(datExpr), method = "ward.D2")
plot(datExpr.tree, 
     main = "Sample clustering", sub="", xlab="", ylab="Hight", 
     cex.lab = 2, cex.axis = 1, cex.main = 1, cex.lab = 1)

#添加样品信息
group <-  factor(c(rep('MII', 2), rep('2C', 4), rep('4C', 3), 
                   rep('8C', 4), rep('16C', 3), rep('Morula', 3), 
                   rep('EarlyB', 3), rep('ExpB', 4)), 
                 levels = c('MII','2C','4C','8C','16C','Morula','EarlyB','ExpB'))
datTraits <- data.frame(sample = rownames(datExpr), group = group)
rownames(datTraits) <- datTraits[,1]
head(datTraits)

# 批量计算软阈值
powers <- c(c(1:10), seq(from = 10, to = 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
best_beta <- sft$powerEstimate
best_beta

# 无标度网络的评估
par(cex = 2.0)
cutcex = 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     labels = powers, cex = cutcex, col = "green3")
abline(h = cutcex, col = "red")

# Soft threshold与平均连通性
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", type = "n", 
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cutcex, 
     col = "red")

# 一步法构建共表达矩阵
corType = "pearson"
cor <- WGCNA::cor
net <- blockwiseModules(datExpr, 
                        power = best_beta, 
                        maxBlockSize = ncol(datExpr), 
                        TOMType = "unsigned", minModuleSize = 50, 
                        reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = T, pamRespectsDendro = F, 
                        saveTOMs = F, saveTOMFileBase = "TOM", 
                        maxPOutliers = ifelse(corType == "pearson", 1, 0.05),
                        verbose = 3)

table(net$colors)
cor <- stats::cor

#层级聚类树展示模块
moduleColors <- labels2colors(net$colors)
table(moduleColors)
plotDendroAndColors(net$dendrograms[[1]], 
                    moduleColors[net$blockGenes[[1]]],
                    "Module", dendroLabels = F, hang = 0.03, 
                    addGuide = T, guideHang = 0.05)

#针对样品矩阵添加对应颜色
sample.colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                colors = c("grey","blue","red","green"), 
                                signed = F)

#系统聚类树及性状热图
plotDendroAndColors(datExpr.tree, sample.colors, 
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8, 
                    marAll = c(1,2,1,1), cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")

#模块和性状的关系
design <- model.matrix(~0 + datTraits$group)
colnames(design) <- levels(as.factor(datTraits$group))
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

head(design)
Trait <- as.data.frame(design[, 6])
names(Trait) <- "Morula"
MET <- orderMEs(cbind(MEs, Trait))  #add the weight to existing module eigengenes

# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", 
                      marDendro = c(1,5,1,2), 
                      marHeatmap = c(6, 6, 0, 2), 
                      cex.lab = 0.8, xLabelsAngle = 90)

#模块聚类
par(cex = 0.80)
plotEigengeneNetworks(MET, 
                      "Eigengene dendrogram", 
                      marDendro = c(0, 4, 2, 0), 
                      plotHeatmaps = T)

#性状与模块的热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", 
                      marHeatmap = c(3, 4, 2, 2),
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)

#模块与性状的相关性
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0) 
moduleTraitCor <- cor(MEs, design, use = "p")
nSamples <- nrow(datExpr)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

#相关性热图展示
par(mar = c(3,9,2,1))
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = colnames(design), yLabels = colnames(MEs), ySymbols = names(MEs), 
               colorLabels = T, colors = blueWhiteRed(50), zlim = c(-1, 1), 
               textMatrix = textMatrix, setStdMargins = F, cex.text = 0.83,
               main = paste("Gene expression and embryo stage relationships"))

#模块内基因分析
modnames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modnames, sep="")
names(MMPvalue) = paste("p.MM", modnames, sep="")

head(design)
stage <- as.data.frame(design[, 6])
names(stage) = "Morula"
geneTraitSignificance <- as.data.frame(cor(datExpr, stage, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(stage), sep="")
names(GSPvalue) <- paste("p.GS.", names(stage), sep="")

#select module
module <- "brown"
column <- match(module, modnames)
moduleGenes <- moduleColors==module

par(mar = c(5,5,3,2))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   cex.main = 1.1, cex.lab = 1.1, cex.axis = 1.1, 
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Genes significance for Morula", 
                   col = module, abline = T)
  
#Select module probes
probes <- colnames(datExpr)
modProbes <- probes[moduleGenes]
write.csv(modProbes, file = "Module_brown_m.csv")
  
#set color
library(gplots)
heatcol <- colorpanel(100, 'darkred', 'orange', 'lemonchiffon')
  
#set the random seed for reproducibility
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
nSelect = 400
set.seed(12345)
select <- sample(nGenes, size = nSelect)
dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = best_beta)
selectTOM <- dissTOM[select, select]

#re-cluster clustering tree to a subset of genes.
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]
plotDiss <- selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, 
        main = "Network heatmap plot, selected genes", 
        col = heatcol)
  
#Network heatmap plot of all genes, taking time!
dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = 6)
plotTOM <- dissTOM^7
diag(plotTOM) <- NA
geneTree <- net$dendrograms[[1]]
TOMplot(plotTOM, geneTree, moduleColors, 
        main = "Network heatmap plot, all DEGs",
        col = heatcol)
  
#导出网络用于Cytoscape
TOM <- TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)
module <- "brown"
probes <- colnames(datExpr)
inModule <- (moduleColors==module)
modProbes <- probes[inModule]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                weighted = TRUE, threshold = 0.02, nodeNames = modProbes, 
                                nodeAttr = moduleColors[inModule])

write.csv(cyt$edgeData, file = "CytoscapeInput-edge-brown.csv")
write.csv(cyt$nodeData, file = "CytoscapeInput-nodes-brown.csv")
  
# boxplot
library(ggpubr)
library(reshape2)
library(dplyr)

head(fpkmE)
explot <- melt(fpkmE) #长宽数据转换
colnames(explot) <- list("gene", "sample", "value")
head(explot)

explot$group <- case_when(str_sub(explot$sample, 1, 3) == "MII"~"MII",
                          str_sub(explot$sample, 1, 2) == "2C"~"2C",
                          str_sub(explot$sample, 1, 2) == "4C"~"4C",
                          str_sub(explot$sample, 1, 2) == "8C"~"8C",
                          str_sub(explot$sample, 1, 2) == "16"~"16C",
                          str_sub(explot$sample, 1, 3) == "Mor"~"Morula",
                          str_sub(explot$sample, 1, 2) == "EB"~"EarlyB",
                          str_sub(explot$sample, 1, 2) == "HB"~"ExpB")

explot$lg <- log2(explot$value + 1)

plot <- explot[explot$gene == "ZSCAN4", ]

ggboxplot(plot, x = "group", y = "lg", fill = "group", palette = "aass") + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1)), title = "ZSCAN4") + 
  theme_bw() + #ylim(0, 10) +
  theme(legend.position = "none",
        axis.title =  element_text(size = 13, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 0.8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 13, color = "black", vjust = 1, hjust = 1, angle = 45))

    
    
