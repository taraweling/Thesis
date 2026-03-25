# Libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
#library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(readr)
library(scales)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library(eulerr)

df <- read_csv("DEGData.csv", stringsAsFactors=FALSE)
#Within Study Pooled Disorders

## INDIVIDUAL
df$GENEID   <- as.character(df$GENEID)
df$DISORDER <- as.character(df$DISORDER)
df$STUDY    <- as.character(df$STUDY)
df$TISSUE   <- as.character(df$TISSUE)
df$YEAR     <- as.integer(df$YEAR)
df$LOG2FC   <- as.numeric(df$LOG2FC)
df$PVAL     <- as.numeric(df$PVAL)

genes_per_disorder <- split(df$GENEID, df$DISORDER)

#genes_per_disorder <- lapply(split(df$GENEID, df$DISORDER), unique)

# DOESN'T WORK BELOW
genes_per_disorder <- df %>%
    group_by(DISORDER) %>%
    summarise(GENES = list(unique(GENEID)), .groups = "drop") %>%
    deframe()

#AD  <- genes_per_disorder[["AD"]]
#ADHD  <- genes_per_disorder[["ADHD"]]
#ASD <- genes_per_disorder[["ASD"]]
#MDD  <- genes_per_disorder[["MDD"]]
BD  <- genes_per_disorder[["BD"]]
#OCD  <- genes_per_disorder[["OCD"]]
SZ  <- genes_per_disorder[["SZ"]]

# volcano plots
bd <- subset(df, DISORDER == "BD")
plot(bd$LOG2FC,
     -log10(bd$PVAL),
     pch = 16,
     xlab = "Log2 Fold Change",
     ylab = "-log(P-value)",
     main = "BD Differential Expression")


# Find Overlap

BDSZ = list(BD,SZ) 
names(BDSZ)=c("Bipolar","Schizophrenia")
comp2 <- ggVennDiagram(BDSZ, label = "both") 
  # comp2+ ggplot2::scale_fill_gradient(low = "yellow", high = "red")
comp2 + 
  #scale_fill_viridis_c(begin = 0.6, end = 0.9, option = "D", direction = -1) +
  scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

MDDBD=list(MDD,BD) 
names(MDDBD)=c("Depression","Bipolar")
comp3 <- ggVennDiagram(MDDBD, label = "both") 
comp3 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 


ADHDMDD=list(ADHD,MDD) 
names(ADHDMDD)=c("ADHD","Depression")
comp4 <- ggVennDiagram(ADHDMDD, label = "both") 
comp4 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

BDADHD=list(BD,ADHD) 
names(BDADHD)=c("Bipolar","ADHD")
comp5 <- ggVennDiagram(BDADHD, label = "both") 
comp5 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

OCDAD=list(OCD,AD) 
names(OCDAD)=c("OCD","Anxiety")
comp6 <- ggVennDiagram(OCDAD, label = "both") 
comp6 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

ASDADHD=list(ASD,ADHD) 
names(ASDADHD)=c("Autism","ADHD")
comp7 <- ggVennDiagram(ASDADHD, label = "both") 
comp7 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

ADHDSZ=list(ADHD,SZ) 
names(ADHDSZ)=c("ADHD","Schizophrenia")
comp8 <- ggVennDiagram(ADHDSZ, label = "both") 
comp8 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

SZASD=list(SZ,ASD) 
names(SZASD)=c("Schizophrenia", "Autism")
comp9 <- ggVennDiagram(SZASD, label = "both") 
comp9 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

OCDASD=list(OCD,ASD) 
names(OCDASD)=c("OCD", "Autism")
comp10 <- ggVennDiagram(OCDASD, label = "both") 
comp10 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 


# Graph of overlap by genetic correlation (linkage value from Brainstorm Consortium et al. (2018))
disorders <- c("AD-MDD","BD-SZ","BD-MDD","ADHD-MDD","BD-ADHD","AD-OCD","ADHD-ASD","ADHD-SZ","ASD-SZ","ASD-OCD")
rval <- c(0.80,0.70,0.47,0.42,0.40,0.40,0.36,0.36,0.20,0.20)
overlap <- c(0.03,0.74,0.1,0.06,0.1,0.02,0.03,0.1,0.06,0.04)
combo_df <- data.frame(X = disorders, Y = rval, Z = overlap)
print(combo_df)

ggplot(combo_df, aes(x = rval, y = overlap, label = disorders)) +
    geom_point(size = 4) +
    geom_text(nudge_y = 0.02) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        x = "Genetic Correlation (r)",
        y = "Percent Gene Overlap"
    )


ggplot(combo_df, aes(x = disorders, y = rval, fill=overlap)) +
    geom_col(position = "dodge") +
    scale_color_gradient(
        low = "red",
        high = "blue",
        labels = scales::percent,    # Format the labels as percentages
        limits = c(0, 1),             # Optional: ensure the scale ranges from 0% to 100%
        name = "Percentage Overlap", # Set the legend title
    )
    labs(title = "Comparison of genetic correlation and DEG overlap", y = "rval", x = "Disorder Pairs")

# Triplets
SZBDADHD=list(SZ,BD,ADHD) 
names(SZBDADHD)=c("Schizophrenia", "Bipolar","ADHD")
comp11 <- ggVennDiagram(SZBDADHD, label = "both") 
comp11 + scale_fill_distiller(palette = "Blues", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

BDMDDADHD=list(BD,MDD,ADHD) 
names(BDMDDADHD)=c("Bipolar","Depression","ADHD")
comp12 <- ggVennDiagram(BDMDDADHD, label = "both") 
comp12 + scale_fill_distiller(palette = "Greens", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

ASDADHDSZ=list(ASD,ADHD,SZ) 
names(ASDADHDSZ)=c("Autism","ADHD","Schizophrenia")
comp13 <- ggVennDiagram(ASDADHDSZ, label = "both", digits = 3) 
comp13 + scale_fill_distiller(palette = "Purples", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

ASDBDSZMDD=list(ASD,BD,SZ,MDD,AD,ADHD,OCD) 
names(ASDBDSZMDD)=c("Autism","Bipolar","Schizophrenia","Depression","Anxiety","ADHD","OCD")
comp14 <- ggVennDiagram(ASDBDSZMDD, label = "count", digits = 3) 
comp14 + scale_fill_distiller(palette = "Greys", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

# Akula et al. Rep
BDMDDSZ=list(BD,MDD,SZ) 
names(BDMDDSZ)=c("Bipolar","Depression","Schizophrenia")
comp15 <- ggVennDiagram(BDMDDSZ, label = "both", digits = 3) 
comp15 + scale_fill_distiller(palette = "Reds", direction = 1) + 
  scale_x_continuous(expand = expansion(mult = .3)) +
  guides(fill = FALSE) 

#
# disorder pairs
d <- c("MDD-AD","BD-SZ","MDD-BD","ADHD-MDD","BD-ADHD","OCD-AD","ASD-ADHD","ADHD-SZ","SZ-ASD","OCD-ASD")
# list of number of overlapping symptoms
#s <- c("")
# r value from Brainstorm Consortium et al
r <- c()
# number of overlapping DEGs

# percent of overlapping DEGs

#plot(venn(PooledDisorders))
#plot(euler(BDSZ), quantities = TRUE)
# label_alpha = 0, label = "both",
#set_size = 0.2) + ggplot2::scale_fill_gradient(low = "yellow", high = "red")
####################

ASD_ADHD_AD <- PooledDisorders[c("ASD", "ADHD", "AD")]
MDD_AD <- PooledDisorders[c("MDD", "AD")]
MDD_AD <- list(MDD = PooledDisorders[c("MDD")], AD = PooledDisorders[c("AD")])
#view(MDD_AD)

ggVennDiagram(MDD_AD, title = "test",
  #data(list = "ADHD, BD"),
  category.names = c("MDD", "AD"),
  label_alpha = 0, label = "both",
  set_size = 0.2) +
  ggplot2::scale_fill_gradient(low = "yellow", high = "red")

  
#title(main="Major Depressive Disorder Overlap With Anxiety Disorder")

