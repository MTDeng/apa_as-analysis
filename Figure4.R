rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing/apa/sheep")

####################expression of apa genes in sheep###########################
load("sheep.filiter.fpkm.rda")
paf <- read.csv("/Users/mtdeng/Documents/Research/Splicing/apa/APA.factors.csv", header = F)
paf_exp <- fpkmE[rownames(fpkmE) %in% paf$V1,]
head(paf_exp)

library(pheatmap)
bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
pheatmap(paf_exp, legend = 10, fontsize = 12, show_rownames = T, 
         scale = "row", cluster_col = F,  cluster_rows = F, angle_col = 45,
         color = c(colorRampPalette(colors = c("royalblue", "white"))(length(bk)/2),
                   colorRampPalette(colors = c("white", "red"))(length(bk)/2)),
         legend_breaks = seq(-2, 2, 1), breaks = bk)

paf_exp$gene <- rownames(paf_exp)

library(reshape2)
library(stringr)
library(dplyr)
paf_exp <- melt(paf_exp)
colnames(paf_exp) <- c("gene", "sample", "fpkm")
head(paf_exp)
table(paf_exp$sample)
paf_exp$group <- case_when(str_sub(paf_exp$sample, 1, 3) == "MII"~"1MII",
                           str_sub(paf_exp$sample, 1, 2) == "2C"~"2C",
                           str_sub(paf_exp$sample, 1, 2) == "4C"~"4C",
                           str_sub(paf_exp$sample, 1, 2) == "8C"~"8C",
                           str_sub(paf_exp$sample, 1, 2) == "16"~"a16C",
                           str_sub(paf_exp$sample, 1, 3) == "Mor"~"bMorula",
                           str_sub(paf_exp$sample, 1, 2) == "EB"~"cEarlyB",
                           str_sub(paf_exp$sample, 1, 2) == "HB"~"dExpB")

library(ggpubr)

paf_exp$exp <- log2(paf_exp$fpkm + 1)
ggboxplot(paf_exp, x = "gene", y = "exp", fill = "group", palette = "aass", title = "") + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + 
  theme_bw() +
  theme(legend.position = "right",
        axis.title =  element_text(size = 10, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 9, color = "black", vjust = 1, hjust = 1, angle = 45))

###########################PAS in sheep#######################################
#read data
library(dplyr)
dat_E8_16 <- read.table("E8&16_All_Prediction_Results.txt", header = T, row.names = 1, sep = '\t', quote = "")
dim(dat_E8_16)
head(dat_E8_16[, (ncol(dat_E8_16)-5): ncol(dat_E8_16)])
dat_E8_16 <- dat_E8_16[, c(3,(ncol(dat_E8_16)-5): ncol(dat_E8_16))]
dat_E8_16 <- dat_E8_16[, c(1:5)]
head(dat_E8_16)
dat_E8_16 <- dat_E8_16 %>% dplyr::rename("E8" = "Group_A_Mean_PDUI", 
                                         "E16" = "Group_B_Mean_PDUI",
                                         "E8_E16" = "PDUI_Group_diff")
head(dat_E8_16)
write.table(dat_E8_16, "PDUI-Pval.E8_E16.txt", row.names = T, col.names = T, sep = '\t', quote = F)
system(paste("sed -i '' '1s/^/Gene\t/'", "PDUI-Pval.E8_E16.txt"))
system("less -S PDUI-Pval.E8_E16.txt | sed -e '1s/^Gene\t/gene_id\tgene_name\tchr\tstrand\t/' -e 's/|/\t/g' | 
       cut -f1,2,5,6,7,8,9 > PDUI-Pval.E8_E16.new.txt") #拆分row name中的｜


dat_E16_M <- read.table("E16&M_All_Prediction_Results.txt", header = T, row.names = 1, sep = '\t', quote = "")
dim(dat_E16_M)
head(dat_E16_M[, (ncol(dat_E16_M)-5): ncol(dat_E16_M)])
dat_E16_M <- dat_E16_M[, c(3, (ncol(dat_E16_M)-5): ncol(dat_E16_M))]
dat_E16_M <- dat_E16_M[, c(1:5)]
head(dat_E16_M)
dat_E16_M <- dat_E16_M %>% dplyr::rename("E16" = "Group_A_Mean_PDUI", 
                                         "M" = "Group_B_Mean_PDUI",
                                         "E16_M" = "PDUI_Group_diff")
head(dat_E16_M)

table(rownames(dat_E8_16) %in% rownames(dat_E16_M))

write.table(dat_E16_M, "PDUI-Pval.E16_M.txt", row.names = T, col.names = T, sep = '\t', quote = F)

system(paste("sed -i '' '1s/^/Gene\t/'", "PDUI-Pval.E16_M.txt"))
system("less -S PDUI-Pval.E16_M.txt | sed -e '1s/^Gene\t/gene_id\tgene_name\tchr\tstrand\t/' -e 's/|/\t/g' | 
       cut -f1,2,5,6,7,8,9 > PDUI-Pval.E16_M.new.txt") #拆分row name中的｜

#################################diff PAS in E16_E8#############################
PDUI_E8_16 <- read.table("PDUI-Pval.E8_E16.new.txt", row.names = 1, header = T)
head(PDUI_E8_16)

PDUI_E8_16$change <- ifelse(PDUI_E8_16$P_val > 0.05, 'stable',
                         ifelse(PDUI_E8_16$P_val < 0.05 & PDUI_E8_16$E8_E16 <= -0.05, 'shorten', 
                                ifelse(PDUI_E8_16$P_val < 0.05 & PDUI_E8_16$E8_E16 >= 0.05, 'lengthen', 'stable')))

table(PDUI_E8_16$change)
write.table(PDUI_E8_16, "PDUI-diff.E8_E16.txt", row.names = T, col.names = T, sep = '\t', quote = F)


dat <- PDUI_E8_16[PDUI_E8_16$change == "shorten" | PDUI_E8_16$change == "lengthen",]
dat <- na.omit(dat)
write.csv(dat, file = "dat.csv")

#plot of diff PAS
library(ggplot2)
table(is.na(PDUI_E8_16$change))
PDUI_E8_16 <- na.omit(PDUI_E8_16)
head(PDUI_E8_16)

ggplot(PDUI_E8_16, aes(x = E8, y = E16, color = change)) + 
  geom_point(size = 0.5) + 
  labs(x = bquote(Lengthened~PAS~of~E8C), y = bquote(Lengthened~PAS~of~E16C)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.2),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_color_manual(values = c("#0000FF90", "#FF000090", "grey90")) +
  annotate("text", x = 0.7, y = 0.01, label = paste("n = 615"), color = "#0000FF90") +
  annotate("text", x = 0.2, y = 0.99, label = paste("n = 227"), color = "#FF000090")

ggplot(PDUI_E8_16, aes(x = -E8_E16, y = -log10(P_val), color = change)) + 
  geom_point(size = 0.5) + theme_bw() + xlim(-1, 0.6) +
  labs(x = bquote(PAS~of~E16C/E8C), y = bquote(-log[10](P_val))) +
  geom_hline(yintercept = 1.3, linetype = 2) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = 2) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.2),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_color_manual(values = c("tomato", "royalblue", "grey90")) +
  annotate("text", x = -0.75, y = 160, label = paste("n = 615"), color = "tomato2") +
  annotate("text", x = 0.35, y = 160, label = paste("n = 227"), color = "royalblue")

##gene expression with PAS
load("/Users/mtdeng/Documents/Research/Splicing/sheep.filiter.fpkm.rda")

PDUI <- PDUI_E8_16[!duplicated(PDUI_E8_16$gene_name),]
rownames(PDUI) <- PDUI$gene_name

expr <- merge(PDUI, fpkmE, by = 0)
head(expr)

expr <- expr[, c(2:3, 8, 18:24)]
expr$E8 <- rowMeans(expr[,c(4:7)])
expr$E16 <- rowMeans(expr[,c(8:10)])
head(expr)

library(ggpubr)
ggboxplot(expr, x = "change", y = "log2(E16+1)", notch = T,
          fill = "change", palette = "aass", outlier.shape = 26) + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + #ylim(0, 40) +
  stat_compare_means(label = "p.signif", method = 't.test', 
                     ref.group = "lengthen", label.y = 11.5) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title =  element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45)) +
  scale_x_discrete(labels = c("stable", "shorten", "lengthen"))

#############################diff PAS in M_E16##################################
PDUI_E16_M <- read.table("PDUI-Pval.E16_M.new.txt", row.names = 1, header = T)
head(PDUI_E16_M)

PDUI_E16_M$change <- ifelse(PDUI_E16_M$P_val > 0.05, 'stable',
                            ifelse(PDUI_E16_M$P_val < 0.05 & PDUI_E16_M$E16_M <= -0.05, 'shorten', 
                                   ifelse(PDUI_E16_M$P_val < 0.05 & PDUI_E16_M$E16_M >= 0.05, 'lengthen', 'stable')))

table(PDUI_E16_M$change)
write.table(PDUI_E16_M, "PDUI-diff.E16_M.txt", row.names = T, col.names = T, sep = '\t', quote = F)

#plot of diff PAS
library(ggplot2)
table(is.na(PDUI_E16_M$change))
PDUI_E16_M <- na.omit(PDUI_E16_M)
head(PDUI_E16_M)

table(PDUI_E16_M$change)
ggplot(PDUI_E16_M, aes(x = E16, y = M, color = change)) + 
  geom_point(size = 0.5) + 
  labs(x = bquote(Lengthened~PAS~of~E16C), y = bquote(Lengthened~PAS~of~M)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.2),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_color_manual(values = c("#0000FF90", "#FF000090", "grey90")) +
  annotate("text", x = 0.7, y = 0.01, label = paste("n = 294"), color = "#0000FF90") +
  annotate("text", x = 0.2, y = 0.99, label = paste("n = 1058"), color = "#FF000090")

ggplot(PDUI_E16_M, aes(x = -E16_M, y = -log10(P_val), color = change)) + 
  geom_point(size = 0.5) + theme_bw() + xlim(-1, 1) +
  labs(x = bquote(PAS~of~M/E16C), y = bquote(-log[10](P_val))) +
  geom_hline(yintercept = 1.3, linetype = 2) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = 2) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.2),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_color_manual(values = c("royalblue", "tomato", "grey90")) +
  annotate("text", x = -0.55, y = 230, label = paste("n = 294"), color = "royalblue") +
  annotate("text", x = 0.55, y = 230, label = paste("n = 1058"), color = "tomato2")

#PAS with gene expression
PDUI <- PDUI_E16_M[!duplicated(PDUI_E16_M$gene_name),]
rownames(PDUI) <- PDUI$gene_name

expr <- merge(PDUI, fpkmE, by = 0)
head(expr)

expr <- expr[, c(2:3, 8, 22:27)]
expr$E16 <- rowMeans(expr[,c(4:6)])
expr$M<- rowMeans(expr[,c(7:9)])
head(expr)

library(ggpubr)
ggboxplot(expr, x = "change", y = "log2(M+1)", fill = "change", 
          palette = "aass", outlier.shape = 26, notch = T) + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + #ylim(0, 40) +
  stat_compare_means(label = "p.signif", method = 't.test', 
                     ref.group = "stable", label.y = 11.5) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title =  element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1.05),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))


###############################PAS in mice######################################
rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing/apa/mice")

dat_E2_L2 <- read.table ("E2_L2_All_Prediction_Results.txt", header = T, row.names = 1, sep = '\t', quote = "")
dim(dat_E2_L2)
head(dat_E2_L2[, (ncol(dat_E2_L2)-5): ncol(dat_E2_L2)])
table(dat_E2_L2$Pass_Filter)

dat_E2_L2 <- dat_E2_L2[, (ncol(dat_E2_L2)-5): ncol(dat_E2_L2)]
dat_E2_L2 <- dat_E2_L2[, c(1:4)]
head(dat_E2_L2)
dat_E2_L2 <- dat_E2_L2 %>% dplyr::rename("E2" = "Group_A_Mean_PDUI", 
                                         "L2" = "Group_B_Mean_PDUI",
                                         "E2_L2" = "PDUI_Group_diff")
head(dat_E2_L2)
write.table(dat_E2_L2, "PDUI-Pval.E2_L2.txt", sep = "\t", quote = F)

system(paste("sed -i '' '1s/^/Gene\t/'", "PDUI-Pval.E2_L2.txt"))
system("less -S PDUI-Pval.E2_L2.txt | sed -e '1s/^Gene\t/gene_id\tgene_name\tchr\tstrand\t/' -e 's/|/\t/g' | 
       cut -f1,2,5,6,7,8 > PDUI-Pval.E2_L2.new.txt")

dat_PN5_E2 <- read.table("PN5_E2_All_Prediction_Results.txt", header = T, row.names = 1, sep = '\t', quote = "")
dim(dat_PN5_E2)
head(dat_PN5_E2[, (ncol(dat_PN5_E2)-5): ncol(dat_PN5_E2)])
dat_PN5_E2 <- dat_PN5_E2[, c(3, (ncol(dat_PN5_E2)-5): ncol(dat_PN5_E2))]
dat_PN5_E2 <- dat_PN5_E2[, c(1:5)]
head(dat_PN5_E2)
dat_PN5_E2 <- dat_PN5_E2 %>% dplyr::rename("PN5" = "Group_A_Mean_PDUI",
                                           "E2" = "Group_B_Mean_PDUI",
                                           "PN5_E2" = "PDUI_Group_diff")
head(dat_PN5_E2)

table(rownames(dat_PN5_E2) %in% rownames(dat_E2_L2))

write.table(dat_PN5_E2, "PDUI-Pval.PN5_E2.txt", row.names = T, col.names = T, sep = '\t', quote = F)

system(paste("sed -i '' '1s/^/Gene\t/'", "PDUI-Pval.PN5_E2.txt"))
system("less -S PDUI-Pval.PN5_E2.txt | sed -e '1s/^Gene\t/gene_id\tgene_name\tchr\tstrand\t/' -e 's/|/\t/g' | 
       cut -f1,2,5,6,7,8,9 > PDUI-Pval.PN5_E2.new.txt") #拆分row name中的｜

#################################diff PAS in PN5_E2#############################
PDUI_PN5_E2 <- read.table("PDUI-Pval.PN5_E2.new.txt", row.names = 1, header = T)
head(PDUI_PN5_E2)

PDUI_PN5_E2$change <- ifelse(PDUI_PN5_E2$P_val > 0.05, 'stable',
                            ifelse(PDUI_PN5_E2$P_val < 0.05 & PDUI_PN5_E2$PN5_E2 <= -0.05, 'shorten', 
                                   ifelse(PDUI_PN5_E2$P_val < 0.05 & PDUI_PN5_E2$PN5_E2 >= 0.05, 'lengthen', 'stable')))

table(PDUI_PN5_E2$change)
write.table(PDUI_PN5_E2, "PDUI-diff.PN5_E2.txt", row.names = T, col.names = T, sep = '\t', quote = F)

#plot of diff PAS
library(ggplot2)
table(is.na(PDUI_PN5_E2$change))
PDUI_PN5_E2 <- na.omit(PDUI_PN5_E2)
head(PDUI_PN5_E2)
table(PDUI_PN5_E2$change)

ggplot(PDUI_PN5_E2, aes(x = -PN5_E2, y = -log10(P_val), color = change)) + 
  geom_point(size = 0.5) + theme_bw() + 
  labs(x = bquote(PAS~of~E2C/PN5), y = bquote(-log[10](P_val))) +
  geom_hline(yintercept = 1.3, linetype = 2) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = 2) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.2),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_color_manual(values = c( "royalblue", "tomato", "grey90")) +
  annotate("text", x = -0.2, y = 100, label = paste("n = 248"), color ="royalblue") +
  annotate("text", x = 0.35, y = 100, label = paste("n = 897"), color = "tomato2")


####################expression of apa genes in mice############################
library(stringr)
count <- read.table("/Users/mtdeng/Documents/Research/Splicing/mammal/mice.count.txt", header = T)
head(count)
count$Geneid <- str_sub(count$Geneid, 1, 18)
table(duplicated(count$Geneid))
rownames(count) <- count[, 1]
count <- count[, -1]

library(org.Mm.eg.db)
columns(org.Mm.eg.db)
symbol <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(count), keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID"))
symbol[1:6,]
symbol <- symbol[!duplicated(symbol$ENSEMBL),]
rownames(symbol) <- symbol$ENSEMBL
count <- merge(symbol, count, by = 0)
count <- count[, -c(1,2,4)]
count <- count %>%
  group_by(SYMBOL) %>%
  summarise_all(mean)
count <- as.data.frame(count)
count <- na.omit(count)
rownames(count) <- count$SYMBOL
head(count)
count <- count[, -1]

# count to fpkm
kb <- count$Length/1000;head(kb)
rpk <- count[,-1]/kb;head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
fpkm <- t(t(rpk)/colSums(count[,-1]) * 10^6);head(fpkm)
head(fpkm)
fpkm <- round(fpkm, 4)

count <- count[, -1]

save(count, file = "mice.exp.count.rda")
save(fpkm, file = "mice.exp.fpkm.rda")


library(DESeq2)
group <- factor(c(rep('GO', 2), rep('MII', 2), rep('PN5', 2), rep('E2', 2), 
                  rep('L2', 2), rep('E4', 2), rep('E8', 2), rep('ICM', 2)), 
                levels = c('GO', 'MII', 'PN5', 'E2', 'L2', 'E4', 'E8', 'ICM'))
colData <- data.frame(row.names = colnames(count), group)
colData

count <- round(count)
dds <- DESeqDataSetFromMatrix(count, colData, design = ~ group)
dds <- DESeq(dds)

E2_PN5 <- results(dds, contrast = c("group", "E2", "PN5"))
table(E2_PN5$padj < 0.05)
E2_PN5 <- cbind(E2_PN5, fpkm)
E2_PN5 <- as.data.frame(E2_PN5)
E2_PN5$change <- ifelse(E2_PN5$padj > 0.05, 'stable',
                        ifelse(E2_PN5$log2FoldChange > 1, 'up', 
                               ifelse(E2_PN5$log2FoldChange < -1, 'down', 'stable')))
table(E2_PN5$change)

L2_E2 <- results(dds, contrast = c("group", "L2", "E2"))
table(L2_E2$padj < 0.05)
L2_E2 <- cbind(L2_E2, fpkm)
L2_E2 <- as.data.frame(L2_E2)
L2_E2$change <- ifelse(L2_E2$padj > 0.05, 'stable',
                       ifelse(L2_E2$log2FoldChange > 1, 'up', 
                              ifelse(L2_E2$log2FoldChange < -1, 'down', 'stable')))
table(L2_E2$change)

##################gene expression with PAS###################
load("/Users/mtdeng/Documents/Research/Splicing/mice.exp.fpkm.rda")

PDUI <- PDUI_PN5_E2[!duplicated(PDUI_PN5_E2$gene_name),]
rownames(PDUI) <- PDUI$gene_name

expr <- merge(PDUI, fpkm, by = 0)
head(expr)

expr <- expr[, c(2, 8, 13:16)]
expr$PN <- rowMeans(expr[,c(3:4)])
expr$E2 <- rowMeans(expr[,c(5:6)])
head(expr)
table(expr$change)

library(ggpubr)
ggboxplot(expr, x = "change", y = "log2(E2+1)", notch = T,
          fill = "change", palette = "aass", outlier.shape = 26) + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + ylim(0, 11) +
  stat_compare_means(label = "p.signif", method = 't.test', 
                     ref.group = "stable", label.y = 10) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title =  element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45)) +
  scale_x_discrete(labels = c("stable", "shorten", "lengthen"))


PDUI_E2_L2 <- read.table("PDUI-Pval.E2_L2.new.txt", row.names = 1, header = T)
head(PDUI_E2_L2)

PDUI_E2_L2$change <- ifelse(PDUI_E2_L2$P_val > 0.05, 'stable',
                             ifelse(PDUI_E2_L2$P_val < 0.05 & PDUI_E2_L2$E2_L2 <= -0.05, 'shorten', 
                                    ifelse(PDUI_E2_L2$P_val < 0.05 & PDUI_E2_L2$E2_L2 >= 0.05, 'lengthen', 'stable')))

table(PDUI_E2_L2$change)
write.table(PDUI_E2_L2, "PDUI-diff.E2_L2.txt", row.names = T, col.names = T, sep = '\t', quote = F)

ggplot(PDUI_E2_L2, aes(x = -E2_L2, y = -log10(P_val), color = change)) + 
  geom_point(size = 0.5) + theme_bw() + 
  labs(x = bquote(PAS~of~L2C/E2C), y = bquote(-log[10](P_val))) +
  geom_hline(yintercept = 1.3, linetype = 2) + 
  geom_vline(xintercept = c(-0.05, 0.05), linetype = 2) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.2),
        text = element_text(color = "black", size = 12),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5)) +
  scale_color_manual(values = c( "royalblue", "tomato", "grey90")) +
  annotate("text", x = -0.5, y = 190, label = paste("n = 560"), color ="royalblue") +
  annotate("text", x = 0.4, y = 190, label = paste("n = 933"), color = "tomato2")


#plot of diff PAS
table(is.na(PDUI_E2_L2$change))
PDUI_E2_L2 <- na.omit(PDUI_E2_L2)
head(PDUI_E2_L2)
table(PDUI_E2_L2$change)
PDUI <- PDUI_E2_L2[!duplicated(PDUI_E2_L2$gene_name),]
rownames(PDUI) <- PDUI$gene_name

expr <- merge(PDUI, fpkm, by = 0)
head(expr)

expr <- expr[, c(2, 7, 14:17)]
expr$E2 <- rowMeans(expr[,c(3:4)])
expr$L2 <- rowMeans(expr[,c(5:6)])
head(expr)
table(expr$change)

library(ggpubr)
ggboxplot(expr, x = "change", y = "log2(L2+1)", notch = T,
          fill = "change", palette = "aass", outlier.shape = 26) + 
  labs(x = "", y = bquote(Log[2]~(FPKM + 1))) + ylim(0, 13) +
  stat_compare_means(label = "p.signif", method = 't.test', 
                     ref.group = "stable", label.y = 11) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title =  element_text(size = 12, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1.05),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45)) +
  scale_x_discrete(labels = c("stable", "shorten", "lengthen"))







