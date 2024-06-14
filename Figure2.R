rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing/")

load("sheep_exp.raw.rda")
countE <- count[,-c(1:4,33:42)]

countE <- countE[,c(26:28, 5:18, 1:4, 29:32, 19:25)] #re-arrange the order
colnames(countE) <- c("MII1", "MII2", "MII3", "E2C1", "E2C2", "E2C3", "E2C4", "E2C5", "E4C1", "E4C2", "E4C3", 
                      "E4C4", "E8C1", "E8C2", "E8C3", "E8C4", "E8C5", "E16C1", "E16C2", "E16C3", "E16C4",
                      "M1", "M2", "M3", "M4", "EB1", "EB2", "EB3", "HB1", "HB2", "HB3", "HB4")

countE <- countE[, -c(which(colnames(countE) == "MII2"),
                      which(colnames(countE) == "E2C5"), 
                      which(colnames(countE) == "E4C3"),
                      which(colnames(countE) == "E8C3"),
                      which(colnames(countE) == "E16C4"),
                      which(colnames(countE) == "M1"))]

#re-define colnames
colnames(countE) <- c("MII1", "MII2", "E2C1", "E2C2", "E2C3", "E2C4", "E4C1", "E4C2", "E4C3", 
                      "E8C1", "E8C2", "E8C3", "E8C4", "E16C1", "E16C2", "E16C3", 
                      "M1", "M2", "M3", "EB1", "EB2", "EB3", "HB1", "HB2", "HB3", "HB4")
save(countE, file = "sheep.filiter.count.rda")

# DEG analysis ----
group <- factor(c(rep('MII', 2), rep('E2C', 4), rep('E4C', 3), rep('E8C', 4), 
                  rep('E16C', 3), rep('Morula', 3), rep('EarlyB', 3), rep('HB', 4)), 
                levels = c('MII', 'E2C', 'E4C', 'E8C', 'E16C', 'Morula', 'EarlyB', 'HB'))
colData <- data.frame(row.names = colnames(countE), group)
colData

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countE, colData, design = ~ group)
dds <- DESeq(dds)
dds

#E2vsMII ----
E2vsMII <- results(dds, contrast = c("group", "E2C", "MII"))
table(E2vsMII$padj < 0.05)
E2vsMII <- cbind(E2vsMII, fpkmE)
E2vsMII <- as.data.frame(E2vsMII)
E2vsMII$change <- ifelse(E2vsMII$padj > 0.05, 'stable',
                         ifelse(E2vsMII$log2FoldChange > 1, 'up', 
                                ifelse(E2vsMII$log2FoldChange < -1, 'down', 'stable')))
table(E2vsMII$change)

# plot of DEGs in E2vsMII
library(ggplot2)
E2vsMII_p <- E2vsMII
E2vsMII_p$MII <- rowMeans(E2vsMII_p[,c(7, 8)])
E2vsMII_p$E2 <- rowMeans(E2vsMII_p[,c(9, 12)])
head(E2vsMII_p)
dim(E2vsMII_p)

E2vsMII_p <- E2vsMII_p[,c(2, 34, 35)]
E2vsMII_p$group <- 
  ifelse(E2vsMII_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(E2vsMII_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(E2vsMII_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(E2vsMII_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(E2vsMII_p$group))
table(is.na(E2vsMII_p$log2FoldChange))
E2vsMII_p <- na.omit(E2vsMII_p)
E2vsMII_p <- E2vsMII_p[rowSums(E2vsMII_p[,c(2,3)]) > 0, ]
table(E2vsMII_p$group)

ggplot(E2vsMII_p, aes(x = log2(MII + 1), y = log2(E2 + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~MII), y = bquote(Log[2]~(FPKM + 1)~of~E2)) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 11948"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 2431"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 2803"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 1006"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 1492"), color = "red3")

save(E2vsMII, E2vsMII_p, file = "DEGs.E2_MII.rda")

#E4vsE2 ----
E4vsE2 <- results(dds, contrast = c("group", "E4C", "E2C"))
table(E4vsE2$padj < 0.05)
E4vsE2 <- cbind(E4vsE2, fpkmE)
E4vsE2 <- as.data.frame(E4vsE2)
E4vsE2$change <- ifelse(E4vsE2$padj > 0.05, 'stable',
                         ifelse(E4vsE2$log2FoldChange > 1, 'up', 
                                ifelse(E4vsE2$log2FoldChange < -1, 'down', 'stable')))
table(E4vsE2$change)

# plot of DEGs in E4vsE2
E4vsE2_p <- E4vsE2
E4vsE2_p$E4 <- rowMeans(E4vsE2_p[,c(13, 15)])
E4vsE2_p$E2 <- rowMeans(E4vsE2_p[,c(9, 12)])
head(E4vsE2_p)
dim(E4vsE2_p)

E4vsE2_p <- E4vsE2_p[,c(2, 34, 35)]
E4vsE2_p$group <- 
  ifelse(E4vsE2_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(E4vsE2_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(E4vsE2_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(E4vsE2_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(E4vsE2_p$group))
table(is.na(E4vsE2_p$log2FoldChange))
E4vsE2_p <- na.omit(E4vsE2_p)
E4vsE2_p <- E4vsE2_p[rowSums(E4vsE2_p[,c(2,3)]) > 0, ]
table(E4vsE2_p$group)

ggplot(E4vsE2_p, aes(x = log2(E2 + 1), y = log2(E4 + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~E2), y = bquote(Log[2]~(FPKM + 1)~of~E4)) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 14271"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 1360"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 2179"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 907"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 1503"), color = "red3")

save(E4vsE2, E4vsE2_p, file = "DEGs.E4_E2.rda")

#E8vsE4 ----
E8vsE4 <- results(dds, contrast = c("group", "E8C", "E4C"))
table(E8vsE4$padj < 0.05)
E8vsE4 <- cbind(E8vsE4, fpkmE)
E8vsE4 <- as.data.frame(E8vsE4)
E8vsE4$change <- ifelse(E8vsE4$padj > 0.05, 'stable',
                        ifelse(E8vsE4$log2FoldChange > 1, 'up', 
                               ifelse(E8vsE4$log2FoldChange < -1, 'down', 'stable')))
table(E8vsE4$change)

# plot of DEGs in E8vsE4
E8vsE4_p <- E8vsE4
E8vsE4_p$E4 <- rowMeans(E8vsE4_p[,c(13, 15)])
E8vsE4_p$E8 <- rowMeans(E8vsE4_p[,c(16, 19)])
head(E8vsE4_p)
dim(E8vsE4_p)

E8vsE4_p <- E8vsE4_p[,c(2, 34, 35)]
E8vsE4_p$group <- 
  ifelse(E8vsE4_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(E8vsE4_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(E8vsE4_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(E8vsE4_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(E8vsE4_p$group))
table(is.na(E8vsE4_p$log2FoldChange))
E8vsE4_p <- na.omit(E8vsE4_p)
E8vsE4_p <- E8vsE4_p[rowSums(E8vsE4_p[,c(2,3)]) > 0, ]
table(E8vsE4_p$group)

ggplot(E8vsE4_p, aes(x = log2(E4 + 1), y = log2(E8 + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~E4), y = bquote(Log[2]~(FPKM + 1)~of~E8)) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 13696"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 2275"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 1987"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 1215"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 1389"), color = "red3")
save(E8vsE4, E8vsE4_p, file = "DEGs.E8_E4.rda")


#E16vsE8----
E16vsE8 <- results(dds, contrast = c("group", "E16C", "E8C"))
table(E16vsE8$padj < 0.05)
E16vsE8 <- cbind(E16vsE8, fpkmE)
E16vsE8 <- as.data.frame(E16vsE8)
E16vsE8$change <- ifelse(E16vsE8$padj > 0.05, 'stable',
                        ifelse(E16vsE8$log2FoldChange > 1, 'up', 
                               ifelse(E16vsE8$log2FoldChange < -1, 'down', 'stable')))
table(E16vsE8$change)

# plot of DEGs in E16vsE8
E16vsE8_p <- E16vsE8
E16vsE8_p$E16 <- rowMeans(E16vsE8_p[,c(20, 22)])
E16vsE8_p$E8 <- rowMeans(E16vsE8_p[,c(16, 19)])
head(E16vsE8_p)
dim(E16vsE8_p)

E16vsE8_p <- E16vsE8_p[,c(2, 34, 35)]
E16vsE8_p$group <- 
  ifelse(E16vsE8_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(E16vsE8_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(E16vsE8_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(E16vsE8_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(E16vsE8_p$group))
table(is.na(E16vsE8_p$log2FoldChange))
E16vsE8_p <- na.omit(E16vsE8_p)
E16vsE8_p <- E16vsE8_p[rowSums(E16vsE8_p[,c(2,3)]) > 0, ]
table(E16vsE8_p$group)

load("/Users/mtdeng/Documents/Research/Splicing/degs/DEGs.E16_E8.rda")
ggplot(E16vsE8_p, aes(x = log2(E8 + 1), y = log2(E16 + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~'8C'), y = bquote(Log[2]~(FPKM + 1)~of~'16C')) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 11906"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 2192"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 2621"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 1676"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 2085"), color = "red3")

save(E16vsE8, E16vsE8_p, file = "DEGs.E16_E8.rda")

#MvsE16 ----
MvsE16 <- results(dds, contrast = c("group", "Morula", "E16C"))
table(MvsE16$padj < 0.05)
MvsE16 <- cbind(MvsE16, fpkmE)
MvsE16 <- as.data.frame(MvsE16)
MvsE16$change <- ifelse(MvsE16$padj > 0.05, 'stable',
                         ifelse(MvsE16$log2FoldChange > 1, 'up', 
                                ifelse(MvsE16$log2FoldChange < -1, 'down', 'stable')))
table(MvsE16$change)

# plot of DEGs in MvsE16
MvsE16_p <- MvsE16
MvsE16_p$E16 <- rowMeans(MvsE16_p[,c(20, 22)])
MvsE16_p$M <- rowMeans(MvsE16_p[,c(23, 25)])
head(MvsE16_p)
dim(MvsE16_p)

MvsE16_p <- MvsE16_p[,c(2, 34, 35)]
MvsE16_p$group <- 
  ifelse(MvsE16_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(MvsE16_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(MvsE16_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(MvsE16_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(MvsE16_p$group))
table(is.na(MvsE16_p$log2FoldChange))
MvsE16_p <- na.omit(MvsE16_p)
MvsE16_p <- MvsE16_p[rowSums(MvsE16_p[,c(2,3)]) > 0, ]
table(MvsE16_p$group)

load("/Users/mtdeng/Documents/Research/Splicing/degs/DEGs.M_E16.rda")
ggplot(MvsE16_p, aes(x = log2(E16 + 1), y = log2(M + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~'16C'), y = bquote(Log[2]~(FPKM + 1)~of~Morula)) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 6786"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 3125"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 3851"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 2555"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 6106"), color = "red3")
save(MvsE16, MvsE16_p, file = "DEGs.M_E16.rda")

#EBvsM ----
EBvsM <- results(dds, contrast = c("group", "EarlyB", "Morula"))
table(EBvsM$padj < 0.05)
EBvsM <- cbind(EBvsM, fpkmE)
EBvsM <- as.data.frame(EBvsM)
EBvsM$change <- ifelse(EBvsM$padj > 0.05, 'stable',
                        ifelse(EBvsM$log2FoldChange > 1, 'up', 
                               ifelse(EBvsM$log2FoldChange < -1, 'down', 'stable')))
table(EBvsM$change)

# plot of DEGs in EBvsM
EBvsM_p <- EBvsM
EBvsM_p$EB <- rowMeans(EBvsM_p[,c(26, 28)])
EBvsM_p$M <- rowMeans(EBvsM_p[,c(23, 25)])
head(EBvsM_p)
dim(EBvsM_p)

EBvsM_p <- EBvsM_p[,c(2, 34, 35)]
EBvsM_p$group <- 
  ifelse(EBvsM_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(EBvsM_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(EBvsM_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(EBvsM_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(EBvsM_p$group))
table(is.na(EBvsM_p$log2FoldChange))
EBvsM_p <- na.omit(EBvsM_p)
EBvsM_p <- EBvsM_p[rowSums(EBvsM_p[,c(2,3)]) > 0, ]
table(EBvsM_p$group)

ggplot(EBvsM_p, aes(x = log2(M + 1), y = log2(EB + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~Morula), y = bquote(Log[2]~(FPKM + 1)~of~Early~blastocyst)) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 10115"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 4037"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 3978"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 2665"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 2770"), color = "red3")
save(EBvsM, EBvsM_p, file = "DEGs.EB_M.rda")

#HBvsEB ----
HBvsEB <- results(dds, contrast = c("group", "HB", "EarlyB"))
table(HBvsEB$padj < 0.05)
HBvsEB <- cbind(HBvsEB, fpkmE)
HBvsEB <- as.data.frame(HBvsEB)
HBvsEB$change <- ifelse(HBvsEB$padj > 0.05, 'stable',
                       ifelse(HBvsEB$log2FoldChange > 1, 'up', 
                              ifelse(HBvsEB$log2FoldChange < -1, 'down', 'stable')))
table(HBvsEB$change)


# plot of DEGs in HBvsEB
HBvsEB_p <- HBvsEB
HBvsEB_p$EB <- rowMeans(HBvsEB_p[,c(26, 28)])
HBvsEB_p$HB <- rowMeans(HBvsEB_p[,c(29, 32)])
head(HBvsEB_p)
dim(HBvsEB_p)

HBvsEB_p <- HBvsEB_p[,c(2, 34, 35)]
HBvsEB_p$group <- 
  ifelse(HBvsEB_p$log2FoldChange >= 2.321928, "FC>=5",
         ifelse(HBvsEB_p$log2FoldChange >= 1, '5>FC>=2', 
                ifelse(HBvsEB_p$log2FoldChange >= -1, "1/2<= FC<2",
                       ifelse(HBvsEB_p$log2FoldChange >= -2.321928, "1/5<=FC<1/2", 'FC<=1/5'))))

table(is.na(HBvsEB_p$group))
table(is.na(HBvsEB_p$log2FoldChange))
HBvsEB_p <- na.omit(HBvsEB_p)
HBvsEB_p <- HBvsEB_p[rowSums(HBvsEB_p[,c(2,3)]) > 0, ]
table(HBvsEB_p$group)

ggplot(HBvsEB_p, aes(x = log2(EB + 1), y = log2(HB + 1), color = group)) + 
  labs(x = bquote(Log[2]~(FPKM + 1)~of~Early~blastocyst), 
       y = bquote(Log[2]~(FPKM + 1)~of~Expended~blastocyst)) +
  geom_point(size = 0.8) + theme_bw() + 
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.1),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5)) +
  scale_color_manual(values = c("grey80", "steelblue1", "springgreen3", "royalblue", "red3")) +
  annotate("text", x = 4.0, y = 14, label = paste("1/2<=FC<2: 12501"), color = "grey80") +
  annotate("text", x = 4.1, y = 13, label = paste("1/5<=FC<1/2: 3975"), color = "steelblue1") +
  annotate("text", x = 3.7, y = 12, label = paste("5>FC>=2: 2652"), color = "springgreen3") +
  annotate("text", x = 3.5, y = 11, label = paste("FC<=1/5: 2768"), color = "royalblue") +
  annotate("text", x = 3.3, y = 10, label = paste("FC>=5: 1161"), color = "red3")
save(HBvsEB, HBvsEB_p, file = "DEGs.HB_EB.rda")

# heatmap for all DEGs (by fc)----
library(pheatmap)
library(reshape2)
library(dplyr)

#load data
{
  load("DEGs.E2_MII.rda")
  load("DEGs.E4_E2.rda")
  load("DEGs.E8_E4.rda")
  load("DEGs.E16_E8.rda")
  load("DEGs.M_E16.rda")
  load("DEGs.EB_M.rda")
  load("DEGs.HB_EB.rda")
  }

# generate a data frame with only group and gene name
{
  E2vsMII_c <- E2vsMII_p[E2vsMII_p$group != "1/2<= FC<2", ]
  E2vsMII_c$id <- rownames(E2vsMII_c)
  E4vsE2_c <- E4vsE2_p[E4vsE2_p$group != "1/2<= FC<2", ]
  E4vsE2_c$id <- rownames(E4vsE2_c)
  E8vsE4_c <- E8vsE4_p[E8vsE4_p$group != "1/2<= FC<2", ]
  E8vsE4_c$id <- rownames(E8vsE4_c)
  
  E16vsE8_c <- E16vsE8_p[E16vsE8_p$group != "1/2<= FC<2", ]
  E16vsE8_c$id <- rownames(E16vsE8_c)
  MvsE16_c <- MvsE16_p[MvsE16_p$group != "1/2<= FC<2", ]
  MvsE16_c$id <- rownames(MvsE16_c)
  EBvsM_c <- EBvsM_p[EBvsM_p$group != "1/2<= FC<2", ]
  EBvsM_c$id <- rownames(EBvsM_c)
  HBvsEB_c <- HBvsEB_p[HBvsEB_p$group != "1/2<= FC<2", ]
  HBvsEB_c$id <- rownames(HBvsEB_c)

  heatdat <- full_join(E2vsMII_c, E4vsE2_c, by = "id")
  heatdat <- full_join(heatdat, E8vsE4_c, by = "id")
  heatdat <- full_join(heatdat, E16vsE8_c, by = "id")
  heatdat <- full_join(heatdat, MvsE16_c, by = "id")
  heatdat <- full_join(heatdat, EBvsM_c, by = "id")
  heatdat <- full_join(heatdat, HBvsEB_c, by = "id")
  head(heatdat)
  }

heatdat <- heatdat[, -c(8, 11, 16, 19, 24, 27)]
colnames(heatdat) <- 
  c("log2FCE2_MII", "MII", "E2", "FCE2_MII", "id",
    "log2FCE4_E2", "E4", "FCE4_E2", "log2FCE8_E4", "E8", "FCE8_E4",
    "log2FCE16_E8", "E16", "FCE16_E8", "log2FCM_E16", "M", "FCM_E16",
    "log2FCEB_M", "EB", "FCEB_M", "log2FCHB_EB", "HB", "FCHB_EB")

rownames(heatdat) <- heatdat$id
heatdat <- heatdat[, -5]

load("sheep.filiter.fpkm.rda")
hm <- fpkmE[rownames(fpkmE) %in% rownames(heatdat), ]
hm <- as.data.frame(hm)

#heatmap using mean fpkm in each group
{
  hm$MII <- rowMeans(hm[,c(1:2)])
  hm$E2 <- rowMeans(hm[,c(3:6)])
  hm$E4 <- rowMeans(hm[,c(7:9)])
  hm$E8 <- rowMeans(hm[,c(10:13)])
  hm$E16 <- rowMeans(hm[,c(14:16)])
  hm$M <- rowMeans(hm[,c(17:19)])
  hm$EB <- rowMeans(hm[,c(20:22)])
  hm$HB <- rowMeans(hm[,c(23:26)])
  hm <- hm[, -c(1:26)]
  }

bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
ph <- 
  pheatmap(hm, legend = 10, fontsize = 8, show_rownames = F, 
           scale = "row", cluster_col = F,  cluster_rows = T, angle_col = 45,
           color = c(colorRampPalette(colors = c("royalblue", "white"))(length(bk)/2), 
                     colorRampPalette(colors = c("white", "red1"))(length(bk)/2)),
           legend_breaks = seq(-2, 2, 1), breaks = bk, cutree_rows = 3)

## clusters ----
row_cluster <- cutree(ph$tree_row, k = 3)
neworder <- hm[ph$tree_row$order,]
cluster <- row_cluster[match(rownames(neworder), names(row_cluster))]
data <- data.frame(neworder, cluster)

data$Cluster <- paste0("cluster", data$cluster)
data <- data[,-(ncol(data)-1)]
data$gene <- rownames(data)
data_new <- melt(data)
head(data_new)
table(data$Cluster)

ggplot(data_new, aes(variable, log2(value + 1), group = gene)) + 
  labs( x = "", y = bquote(Log[2]~(FPKM + 1))) +
  geom_line(color = "white", size = 0.6) + 
  geom_hline(yintercept = 1, linetype = 2) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line", size = 0.8, color = "red1") + 
  facet_wrap(Cluster~., nrow = 3, ncol = 1) + #分面，一行多列
  theme_classic() + scale_y_continuous(breaks = seq(0, 4, 1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(color = "black", size = 8),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black", angle = 90, vjust = 0.5))

#sanky-----
library(reshape2)
library(tidyverse)
{ sankydat <- heatdat[, c(4,7,10,13,16,19,22)]
  dim(sankydat)
  head(sankydat)
  
  sankydat$FCE2_MII <- replace_na(sankydat$FCE2_MII, "stable")
  sankydat$FCE4_E2 <- replace_na(sankydat$FCE4_E2, "stable")
  sankydat$FCE8_E4 <- replace_na(sankydat$FCE8_E4, "stable")
  sankydat$FCE16_E8 <- replace_na(sankydat$FCE16_E8, "stable")
  sankydat$FCM_E16 <- replace_na(sankydat$FCM_E16, "stable")
  sankydat$FCEB_M <- replace_na(sankydat$FCEB_M, "stable")
  sankydat$FCHB_EB <- replace_na(sankydat$FCHB_EB, "stable")
  }

sankydat <- sankydat %>% 
  select(colnames(sankydat)) %>% 
  mutate(FCE2_MII = case_when(FCE2_MII == "FC>=5"~"up", 
                              FCE2_MII == "5>FC>=2"~"up",
                              FCE2_MII == "stable"~"stable",
                              FCE2_MII == "1/5<=FC<1/2"~"down", 
                              FCE2_MII == "FC<=1/5"~"down"),
         FCE4_E2 = case_when(FCE4_E2 == "FC>=5"~"up", 
                             FCE4_E2 == "5>FC>=2"~"up",
                             FCE4_E2 == "stable"~"stable",
                             FCE4_E2 == "1/5<=FC<1/2"~"down", 
                             FCE4_E2 == "FC<=1/5"~"down"),
         FCE8_E4 = case_when(FCE8_E4 == "FC>=5"~"up", 
                             FCE8_E4 == "5>FC>=2"~"up",
                             FCE8_E4 == "stable"~"stable",
                             FCE8_E4 == "1/5<=FC<1/2"~"down", 
                             FCE8_E4 == "FC<=1/5"~"down"),
         FCE16_E8 = case_when(FCE16_E8 == "FC>=5"~"up", 
                              FCE16_E8 == "5>FC>=2"~"up",
                              FCE16_E8 == "stable"~"stable",
                              FCE16_E8 == "1/5<=FC<1/2"~"down", 
                              FCE16_E8 == "FC<=1/5"~"down"),
         FCM_E16 = case_when(FCM_E16 == "FC>=5"~"up", 
                             FCM_E16 == "5>FC>=2"~"up",
                             FCM_E16 == "stable"~"stable",
                             FCM_E16 == "1/5<=FC<1/2"~"down", 
                             FCM_E16 == "FC<=1/5"~"down"),
         FCEB_M = case_when(FCEB_M == "FC>=5"~"up", 
                            FCEB_M == "5>FC>=2"~"up",
                            FCEB_M == "stable"~"stable",
                            FCEB_M == "1/5<=FC<1/2"~"down", 
                            FCEB_M == "FC<=1/5"~"down"),
         FCHB_EB = case_when(FCHB_EB == "FC>=5"~"up", 
                             FCHB_EB == "5>FC>=2"~"up",
                             FCHB_EB == "stable"~"stable",
                             FCHB_EB == "1/5<=FC<1/2"~"down", 
                             FCHB_EB == "FC<=1/5"~"down"))

library(ggalluvial)
library(dplyr)

sankydat <- sankydat %>%
  group_by(FCE2_MII, FCE4_E2, FCE8_E4, FCE16_E8, FCM_E16, FCEB_M, FCHB_EB) %>% 
  summarise(count = n())
head(sankydat)

sankydat <- sankydat[sankydat$count >= 5,]
ggplot(as.data.frame(sankydat),
       aes(axis1 = FCE2_MII, axis2 = FCE4_E2, axis3 = FCE8_E4,
           axis4 = FCE16_E8, axis5 = FCM_E16, axis6 = FCEB_M,
           axis7 = FCHB_EB, y= count)) +
  scale_x_discrete(limits = c("FCE2_MII", "FCE4_E2", "FCE8_E4", "FCE16_E8", 
                              "FCM_E16", "FCEB_M", "FCHB_EB"),
                   expand = c(0.2, 0.1)) + geom_flow() +
  geom_alluvium(aes(fill = FCE2_MII)) + 
  geom_stratum(alpha = 0.5) + 
  theme_minimal() 


#signal intensity-----
library(ggpubr)
intensity <- readxl::read_excel("fig2.xlsx", sheet = 5)
head(intensity)
intensity <- as.data.frame(intensity)
intensity <- intensity %>%
  mutate(
    group = case_when(group == "RNAPII-S2-8cell" ~ "8C",
                      group == "RNAPII-S2-16cell" ~ "16C",
                      group == "RNAPII-S2-morula" ~ "Morula"))
table(intensity$group)

mycomp <- list(c("8C", "E16C"), c("8C", "Morula"), c("E16C", "Morula"))

dat <- intensity[intensity$group == "RNAPII_Ser2P", ]
ggboxplot(dat, x = "sample", y = "mean", fill = "sample", palette = "aass", outlier.shape = 26) + 
  labs(x = "", y = "Signal intensity") + ylim(0, 50) +
  stat_compare_means(comparisons = mycomp, label = "p.signif", 
                     label.y = c(40, 45, 40)) + 
  theme_bw() +
  theme(legend.position = "",
        axis.title =  element_text(size = 10, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1.2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 9, color = "black", vjust = 1, hjust = 1, angle = 45))



