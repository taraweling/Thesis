# Libraries
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(readr)
library(scales)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(eulerr)

df <- read_csv("DEGData.csv") 

gene_dict <- split(df[, c("DISORDER", "STUDY", "YEAR", "TISSUE", "LOG2FC", "PVAL")],
                   df$GENEID)
# Each element is a data frame of all records for that gene.

## disorder types


# following Akula, significance of overlapping genes were calcuated
# using SuperExactTest in R.