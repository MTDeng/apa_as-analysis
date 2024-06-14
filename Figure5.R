rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing")

##################sheep as analysis#####################
as.dat <- read.table("splicing/sheep/summary_sheep.txt", 
                     header = T, sep = "\t", check.names = F)
head(as.dat)

#boxplot
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)

plot <- melt(as.dat) #长宽数据转换
head(plot)
colnames(plot) <- list("Group", "EventType", "Change", "Num")
plot$Change <- case_when(plot$Change == "SignificantEventsJC"~"Total",
                         plot$Change == "SigEventsJCSample1HigherInclusion"~"Increased",
                         plot$Change == "SigEventsJCSample2HigherInclusion"~"Decreased")

plot.all <- plot[plot$Change != "Total",]
plot.all$Group <- factor(plot.all$Group, 
                         levels = c("2C_MII", "4C_2C", "8C_4C", "16C_8C", "M_16C", "EB_M"))

ggplot(plot.all, aes(Group, Num, fill = Change)) + 
  labs(x = "", y = "Number of events") + 
  geom_bar(stat = "identity", width = 0.75) + 
  scale_fill_manual(values=c("royalblue", "tomato")) + 
  stat_summary(aes(group = 1), fun = "mean", geom = "line", size = 0.6, color = "black") + 
  facet_wrap(EventType~., nrow = 1, ncol = 5) + #分面，一行多列
  theme_classic() + scale_y_continuous(breaks = seq(0, 4000, 2000)) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(color = "black", size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black", angle = 90, vjust = 0.5, hjust = 0.5))

#################mice as analysis#####################
as.dat <- read.table("splicing/mice/summary_mice.txt", 
                     header = T, sep = "\t", check.names = F)
head(as.dat)
plot <- melt(as.dat) #长宽数据转换
head(plot)
colnames(plot) <- list("Group", "EventType", "Change", "Num")
table(plot$Change)
plot$Change <- case_when(plot$Change == "SignificantEventsJC"~"Total",
                         plot$Change == "SigEventsJCSample1HigherInclusion"~"Increased",
                         plot$Change == "SigEventsJCSample2HigherInclusion"~"Decreased")

plot.all <- plot[plot$Change != "Total",]
table(plot.all$Group)

plot.all$Group <- factor(plot.all$Group, 
                         levels = c("PN5_MII", "E2_PN5", "L2_E2", "4C_L2", "8C_4C"))

ggplot(plot.all, aes(Group, Num, fill = Change)) + 
  labs(x = "", y = "Number of events") + geom_bar(stat = "identity", width = 0.75) + 
  scale_fill_manual(values=c("royalblue1", "tomato2")) + 
  #geom_text(aes(label = Num), size = 3) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line", size = 0.6, color = "black") + 
  facet_wrap(EventType~., nrow = 1, ncol = 5) + #分面，一行多列
  theme_classic() + scale_y_continuous(breaks = seq(0, 2000, 500)) +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(color = "black", size = 10),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
        axis.text.x = element_text(size = 8, color = "black", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, color = "black", angle = 90, vjust = 0.5, hjust = 0.5))

###########expression of as genes##################

as <- readxl::read_xlsx("splicing/AS_pathway.xlsx", sheet = 1)
as <- as.data.frame(as)
load("degs/DEGs.M_E16.rda")
plot.as <- MvsE16[rownames(MvsE16) %in% as$gene_name,]
head(plot.as)
s.as <- plot.as[, c("log2FoldChange", "change")]
head(s.as)

#mice as
library(stringr)
load("mice.exp.count.rda")
head(count)
# DEG analysis ----
group <- factor(c(rep('Growing', 2), rep('MII', 2), rep('Zygote', 2), rep('E2C', 2),
                  rep('L2C', 2), rep('E4C', 2), rep('E8C', 2), rep('ICM', 2)), 
                levels = c('Growing', 'MII', 'Zygote', 'E2C', 'L2C', 'E4C', 'E8C', 'ICM'))
colData <- data.frame(row.names = colnames(count), group)
colData

library(DESeq2)
count <- round(count)
dds <- DESeqDataSetFromMatrix(count, colData, design = ~ group)
dds <- DESeq(dds)
dds

L2vsE2 <- results(dds, contrast = c("group", "L2C", "E2C"))
table(L2vsE2$padj < 0.05)
L2vsE2 <- as.data.frame(L2vsE2)
save(L2vsE2, file = "mice.DEGs.L2vsE2.rda")

head(L2vsE2)
rownames(L2vsE2) <- toupper(rownames(L2vsE2))

m.as <- L2vsE2[rownames(L2vsE2) %in% as$gene_name, ]
head(m.as)
plot.as <- merge(m.as, s.as, by = 0)
head(plot.as)
rownames(plot.as) <- plot.as[, 1]
plot.as <- plot.as[, c(3, 8)]
head(plot.as)
colnames(plot.as) <- c("m.log2fc", "s.log2fc")

plot.as$label <- rownames(plot.as)
plot.as$change <-  
  ifelse(plot.as$m.log2fc >= 1 & plot.as$s.log2fc >= 1, "red",
         ifelse(plot.as$m.log2fc <= -1 & plot.as$s.log2fc <= -1, "blue", "grey"))

change <- c(red = 'red', grey = "grey8", blue = 'royalblue')

library(ggplot2)

ggplot(plot.as, aes(m.log2fc, s.log2fc, fill = change)) +
  geom_point(size = 0.8, aes(color = change)) + theme_bw() +
  geom_hline(yintercept = c(-1, 1), lty = 2, col = "black", lwd = 0.3) + 
  geom_vline(xintercept = c(-1, 1), lty = 2, col = "black", lwd = 0.3) +
  geom_text_repel(data = plot.as, aes(m.log2fc, s.log2fc, label = label, color = change), size = 2.6) +
  scale_color_manual(values = change) + 
  labs(y = bquote(Log[2]~(M~vs.16C)~of~sheep), x = bquote(Log[2]~(L2~vs.E2)~of~mice)) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.05),
        axis.title = element_text(size = 12, color = "black", vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black", angle = 90, hjust = 0.5, vjust = 0.5))

##GSEA-----

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
load("degs/DEGs.M_E16.rda")
gsea <- MvsE16

ENTREID <- bitr(rownames(gsea), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
ENTREID <- ENTREID[!duplicated(ENTREID$SYMBOL),]
rownames(ENTREID) <- ENTREID$SYMBOL

gsea.dat <- merge(ENTREID, gsea, by = 0)
gsea.sort <- gsea.dat[order(gsea.dat$log2FoldChange, decreasing = T),]
head(gsea.sort)

dat <- readxl::read_xlsx("/Users/mtdeng/Documents/Research/Splicing/splicing/venn_result14993.xlsx", 
                         sheet = 3)
s.gene <- gsea.sort[gsea.sort$SYMBOL %in% dat$gene,]
gene_s <- s.gene$log2FoldChange
names(gene_s) <- s.gene$ENTREZID
gene_s <- na.omit(gene_s)

library(enrichplot)
gseGO_all <- gseGO(gene_s, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "all", by = "fgsea",
                   minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 0.05)
gseGO_all@result[["Description"]]
gseaplot2(gseGO_all, c(52, 62, 72), color = rainbow(3), rel_heights = c(0.6, 0.2, 0), 
          pvalue_table = F, ES_geom = "line", base_size = 10, title = "")

library(ReactomePA)
gseKEGG_all <- gsePathway(gene_s, organism = "hsa", by = "fgsea", 
                          minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 0.05)
gseKEGG_all@result[["Description"]]
gseaplot2(gseGO_all, c(47), color = c("royalblue","steelblue"), rel_heights = c(0.6, 0.1, 0.3), 
          pvalue_table = F, ES_geom = "line", base_size = 7, title = "")


load("mice.DEGs.L2vsE2.rda")

library(org.Mm.eg.db)

gsea <- L2vsE2
ENTREID <- bitr(rownames(gsea), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db") 
ENTREID <- ENTREID[!duplicated(ENTREID$SYMBOL),]
rownames(ENTREID) <- ENTREID$SYMBOL

gsea.dat <- merge(ENTREID, gsea, by = 0)
gsea.sort <- gsea.dat[order(gsea.dat$log2FoldChange, decreasing = T),]
head(gsea.sort)

dat <- readxl::read_xlsx("/Users/mtdeng/Documents/Research/Splicing/splicing/venn_result14993.xlsx", 
                         sheet = 2)
s.gene <- gsea.sort[gsea.sort$SYMBOL %in% dat$name,]
gene_s <- s.gene$log2FoldChange
names(gene_s) <- s.gene$ENTREZID
gene_s <- na.omit(gene_s)

library(enrichplot)
gseGO_all <- gseGO(gene_s, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "all", by = "fgsea",
                   minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 0.05)
gseGO_all@result[["Description"]]
gseaplot2(gseGO_all, c(2, 18, 55, 57), color = rainbow(4), rel_heights = c(0.6, 0.2, 0), 
          pvalue_table = F, ES_geom = "line", base_size = 12, title = "")

library(ReactomePA)
gseKEGG_all <- gsePathway(gene_s, organism = "hsa", by = "fgsea", 
                          minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 0.05)
gseKEGG_all@result[["Description"]]
gseaplot2(gseGO_all, c(47), color = c("royalblue","steelblue"), rel_heights = c(0.6, 0.1, 0.3), 
          pvalue_table = F, ES_geom = "line", base_size = 7, title = "")




rm(list = ls())
setwd("/Users/mtdeng/Documents/Research/Splicing")

dat <- readxl::read_xlsx("fig6.xlsx")
head(dat)
dat <- as.data.frame(dat)

library(ggpubr)
dat$dev <- dat$dev*100
ggboxplot(dat, x = "group", y = "dev", fill = "sample", palette = "npj") + 
  labs(x = "", y = "2-cell embryos (%)") + ylim(40, 100) +
  theme_bw() +  scale_fill_manual(values=c("royalblue1", "tomato2")) +
  theme(legend.position = "top",
        axis.title =  element_text(size = 14, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))

dat <- readxl::read_xlsx("fig6.xlsx", sheet = 2)
head(dat)
dat <- as.data.frame(dat)

library(ggpubr)
dat$dev <- dat$dev*100
ggboxplot(dat, x = "group", y = "dev", fill = "group") + 
  labs(x = "", y = "Percentage of blastocysts") + ylim(0, 100) + 
  stat_compare_means(method = 't.test', label = "p.signif", ref.group = "Control", label.y = 90) +
  theme_bw() + scale_fill_manual(values=c("royalblue1", "tomato2")) +
  theme(legend.position = "",
        axis.title =  element_text(size = 14, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))

dat <- readxl::read_xlsx("fig6.xlsx", sheet = 3)
head(dat)
dat <- as.data.frame(dat)

library(ggpubr)
ggboxplot(dat, x = "group", y = "signal", fill = "group", add = 'jitter') + 
  labs(x = "", y = "5-EU signal intensity") + ylim(0, 40) + 
  stat_compare_means(method = 't.test', label.y = 36, label = "p.signif", ref.group = "Control") +
  theme_bw() + scale_fill_manual(values=c("royalblue1", "tomato2")) +
  theme(legend.position = "",
        axis.title =  element_text(size = 14, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))

dat <- readxl::read_xlsx("fig6.xlsx", sheet = 4)
head(dat)
dat <- as.data.frame(dat)

library(ggpubr)
ggboxplot(dat, x = "group", y = "signal", fill = "group", add = 'jitter') + 
  labs(x = "", y = "RNAPII signal intensity") + ylim(0, 30) + 
  stat_compare_means(method = 't.test', label.y = 28,  label = "p.signif", ref.group = "Control") +
  theme_bw() + scale_fill_manual(values=c("royalblue1", "tomato2")) +
  theme(legend.position = "",
        axis.title =  element_text(size = 14, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))

dat <- readxl::read_xlsx("fig6.xlsx", sheet = 5)
head(dat)
dat <- as.data.frame(dat)

library(ggpubr)
ggboxplot(dat, x = "group", y = "signal", fill = "group", add = 'jitter') + 
  labs(x = "", y = "RNAPII-Ser2P signal intensity") + ylim(0, 30) + 
  stat_compare_means(method = 't.test', label.y = 28, label = "p.signif", ref.group = "Control") +
  theme_bw() + scale_fill_manual(values=c("royalblue1", "tomato2")) +
  theme(legend.position = "",
        axis.title =  element_text(size = 14, color = "black"), 
        panel.background = element_rect(fill = NA, colour = " black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, color = "black", angle = 90, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))
