#code for figure 4C in R, using ggplot2

library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)


#load the table with normalised counts from DESeq2
Deseq_GFPp <- read.delim("~/RNA-seq/DeSeq2/mm39/rnaseq_rscripts/outputs/mergedTables/Normalised_reads_table.txt")

View(Deseq_GFPp)

#keep data columns and rows of interest

head(Deseq_GFPp[,c(2,4:17)])

df1 <-Deseq_GFPp[,c(2,4:17)]

#genes of interest list

list <- c("Shox2", "Pitx1", "Prrx1", "Jund", "Sox9",
          "Krt14", "Cldn5", "Mrc1", "Ttn", "Runx2", "Dlx2", "Tbx15", "Egr1", "Lhx9", "Junb")
GOI_df <- filter(df1, gene_short_name %in% list)
GOI_df <- data.frame(GOI_df)
head(GOI_df)




colnames(GOI_df) <- c("genes", "FL_wt_1", "FL_wt_2", "Inv1_1", "Inv1_2", "Inv2_1", "Inv2_2","Rel1_1", "Rel1_2", "Rel2_1", "Rel2_2", "Rel3_1", "Rel3_2", "HL_wt_1", "HL_wt_2")
row.names(GOI_df) <- GOI_df$genes

#log2 ratio transformation

GOI_df$FL_FL <- log2(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))
GOI_df$I1_FL <- log2(rowMeans(GOI_df[,c('Inv1_1', 'Inv1_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))
GOI_df$I2_FL <- log2(rowMeans(GOI_df[,c('Inv2_1', 'Inv2_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))
GOI_df$R1_FL <- log2(rowMeans(GOI_df[,c('Rel1_1', 'Rel1_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))
GOI_df$R2_FL <- log2(rowMeans(GOI_df[,c('Rel2_1', 'Rel2_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))
GOI_df$R3_FL <- log2(rowMeans(GOI_df[,c('Rel3_1', 'Rel3_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))
GOI_df$HL_FL <- log2(rowMeans(GOI_df[,c('HL_wt_1', 'HL_wt_2')])/(rowMeans(GOI_df[,c('FL_wt_1', 'FL_wt_2')])))


GOI_df_3 <- subset(GOI_df, select = c("genes", "FL_FL", "I1_FL", "I2_FL","I1_FL", "R2_FL", "R3_FL", "HL_FL"))

gene_order <- c("Pitx1",  "Junb", "Jund", "Egr1", "Runx2", "Sox9", "Shox2","Tbx15", "Lhx9","Prrx1",
                "Cldn5", "Dlx2", "Krt14", "Ttn", "Mrc1")

GOI_df_3$genes <- factor(GOI_df_3$genes, levels = gene_order)

m_GOI_df_3 <- melt(GOI_df_3)


m_GOI_df_3$genes <- factor(m_GOI_df_3$genes, levels = gene_order)

e <- ggplot(m_GOI_df_3, aes(x = value, y = genes, fill = variable))
e + geom_bar(stat = "identity", position = "dodge", colour="black") + theme_bw()
