# Libraries
library(tidyverse)
library(ggplot2)
#library(hrbrthemes)
library(tm)
library(proustr)
library(readr)
library(scales)
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
#library(ggVennDiagram)
library(eulerr)
#install.packages("BiocManager")
library(BiocManager)
#install() # Install BioConductor core packages
install("maEndToEnd")
install("ArrayExpress")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)e


Disorders <- read_csv("DEGDataSample .csv")


#Within Study Pooled Disorders
#UltraPool <- read_csv("BIO363UltraPool.csv")
#view(PooledDisorders)
#PooledDisorders$ADHD_MDD <- PooledDisorders$ADHD %in% PooledDisorders$MDD
#PooledDisorders$ADHD_MDD <- PooledDisorders[PooledDisorders$MDD %in% PooledDisorders$ADHD]

#Depression x Bipolar Disorder
genes_to_test <- (Disorders$BD)
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
plot(barplot(GO_results,showCategory = 10))
png("MDDBD.png", res = 250, width = 1200, height = 1000)
#png(fit)
#fit
#PooledDisorders$ALL <- df$ADHD == df$MDD & 
# work on GO stuff beloq
#library(clusterProfiler)

OldDisorders <- read_csv("BIO363IPPool.csv")
reference <- (OldDisorders$MDD)
GO_results <- enrichGO(gene = reference, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results)
plot(barplot(GO_results,showCategory = 10))
#png("MDDBD.png", res = 250, width = 1200, height = 1000)
#print(fit)
#fit


# ignore below
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("AnnotationDbi")
# BiocManager::install(version = "3.21")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install(version = "3.21")

