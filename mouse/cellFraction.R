
rm(list=ls())
setwd("~/Documents/work/20240816_LiangZhu/20251111_scRNAseq_MouseModel_YJW/")
source("../../public/WFMyFunction.R")

## ------- 1. major cell types -------

library(readxl)
cellTypeDF <- read_excel("./cellFraction.xlsx", sheet="cellFraction1")
cellTypeDF <- as.data.frame(cellTypeDF)
rownames(cellTypeDF) <- cellTypeDF$Samples

# cellTypeDF <- cellTypeDF[1:9,]s

library(dplyr)
library(tidyr)

df_long <- cellTypeDF %>%
  pivot_longer(
    cols = -c(Samples),
    names_to = "CellType", 
    values_to = "Proportion"
  )

df <- as.data.frame(df_long)
cellColLabels <- c("#5B3178",
                   "#005AA5",
                   "#D8A7C4",
                   "#EA8F22",
                   "#7EB44C", 
                   "#8C564B",
                   "#E15032",
                   "#98C6E5",
                   "#80807F",
                   "#F0BB75",
                   "#EB365A",
                   "#2E5A6F",
                   "#ae88c3")

cellColLabels <- c("#1F77B4", "#D62728", "#AA40FC","#EFD358", "#FF7F0E","#279E68")

df$Samples <- factor(df$Samples, levels = c('WT','TKO'))
df$CellType <- factor(df$CellType, levels = unique(df$CellType))


pdf("./allCellFraction.1.pdf",10,7)
ggplot(df, aes(x = Samples, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cellColLabels) +
  labs(title = "", x = "", y = "Percentage") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black")) + 
  theme(plot.title = element_text(size = 25, colour = "black", hjust=0.5)) +
  theme(legend.text = element_text(size = 25, colour = "black")) +
  theme(legend.title = element_text(size = 25, colour = "black")) + 
  theme(axis.ticks.x = element_line(linewidth = 1, colour = "black"),
        axis.ticks.y = element_line(linewidth = 1, colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 18, colour = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18, colour = "black")) +
  theme(panel.border = element_rect(linewidth = 2.5, colour = "black"))
dev.off()



## 
cellTypeDF <- read_excel("./cellFraction.xlsx", sheet="cellFraction2")
cellTypeDF <- as.data.frame(cellTypeDF)
rownames(cellTypeDF) <- cellTypeDF$Samples

# cellTypeDF <- cellTypeDF[1:9,]s

library(dplyr)
library(tidyr)

df_long <- cellTypeDF %>%
  pivot_longer(
    cols = -c(Samples),
    names_to = "CellType", 
    values_to = "Proportion"
  )

df <- as.data.frame(df_long)
cellColLabels <- c("#5B3178",
                   "#005AA5",
                   "#D8A7C4",
                   "#EA8F22",
                   "#7EB44C", 
                   "#8C564B",
                   "#E15032",
                   "#98C6E5",
                   "#80807F",
                   "#F0BB75",
                   "#EB365A",
                   "#2E5A6F",
                   "#ae88c3")

# cellColLabels <- c("#1F77B4", "#D62728", "#AA40FC","#EFD358", "#FF7F0E",,"#279E68")

df$Samples <- factor(df$Samples, levels = c('WT','TKO'))
df$CellType <- factor(df$CellType, levels = unique(df$CellType))


pdf("./allCellFraction.2.pdf",10,7)
ggplot(df, aes(x = Samples, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cellColLabels) +
  labs(title = "", x = "", y = "Percentage") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black")) + 
  theme(plot.title = element_text(size = 25, colour = "black", hjust=0.5)) +
  theme(legend.text = element_text(size = 25, colour = "black")) +
  theme(legend.title = element_text(size = 25, colour = "black")) + 
  theme(axis.ticks.x = element_line(linewidth = 1, colour = "black"),
        axis.ticks.y = element_line(linewidth = 1, colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 18, colour = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18, colour = "black")) +
  theme(panel.border = element_rect(linewidth = 2.5, colour = "black"))
dev.off()




## ------- 2. subset cell types -------

cellTypeDF <- read_excel("./cellFraction.xlsx", sheet="subMyeloid")
cellTypeDF <- as.data.frame(cellTypeDF)
rownames(cellTypeDF) <- cellTypeDF$Samples

df_long <- cellTypeDF %>%
  pivot_longer(
    cols = -c(Samples),
    names_to = "CellType", 
    values_to = "Proportion"
  )


df <- as.data.frame(df_long)
df$Samples <- factor(df$Samples, levels = c('WT','TKO'))
df$CellType <- factor(df$CellType, levels = unique(df$CellType))


pdf("./allCellFraction.subMyeloid.pdf",10,7)
ggplot(df, aes(x = Samples, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cellColLabels) +
  labs(title = "", x = "", y = "Percentage") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.title.x = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black")) + 
  theme(plot.title = element_text(size = 25, colour = "black", hjust=0.5)) +
  theme(legend.text = element_text(size = 25, colour = "black")) +
  theme(legend.title = element_text(size = 25, colour = "black")) + 
  theme(axis.ticks.x = element_line(linewidth = 1, colour = "black"),
        axis.ticks.y = element_line(linewidth = 1, colour = "black"),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(size = 18, colour = "black", angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18, colour = "black")) +
  theme(panel.border = element_rect(linewidth = 2.5, colour = "black"))
dev.off()

