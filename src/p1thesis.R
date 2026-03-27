install.packages("data.table")
install.packages("tidyverse")
install.packages("dplyr")
install.packages("BiocManager")
install.packages("readr")
install.packages("scales")
install.packages("org.Hs.eg.db")
install.packages("clusterProfiler")
install.packages("stats")
install.packages("data.table")
install.packages("igraph")
install.packages("clusterProfiler")
install.packages("ggraph")
install.packages("ComplexHeatmap")
install.packages("circlize")
install.packages('pheatmap')
install.packages("DOSE")
install.packages("enrichplot")
install.packages("ggupset")
install.packages("jsonlite")

BiocManager::install(version = "3.22") # keep at 3.22 for consistency
BiocManager::install("clusterProfiler", lib="~/R/library")
BiocManager::install("ComplexHeatmap", lib="~/R/library")
BiocManager::install("biomaRt", lib="~/R/library")
BiocManager::install("org.Hs.eg.db", lib="~/R/library")
BiocManager::install("enrichplot", lib="~/R/library")

# Optional cleaning environment before starting :D
#rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
#gc() # free up memory and report the memory usage
#options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation


# Core
library(tidyverse)
library(dplyr)
library(BiocManager)

# Annotation
library(biomaRt)
library(readr)
library(scales)
library(org.Hs.eg.db)
library(clusterProfiler) # PEA analysis

# Statistics
library(stats)
library(data.table)

# Networks
library(igraph)
library(ggraph)

# Visualization
library(ComplexHeatmap)
library(circlize)
library(DOSE)
library(pheatmap)
library(ggupset)
library(RColorBrewer)

# Enrichment
library(enrichplot)

# if not using RStudio, getwd() - if incorrect
# setwd("absolute/path/to/your/project")
# Create output directories
if(!dir.exists("results")) dir.create("results")
if(!dir.exists("results/figures")) dir.create("results/figures")
if(!dir.exists("results/tables")) dir.create("results/tables")
if(!dir.exists("results/enrichment")) dir.create("results/enrichment")
if(!dir.exists("results/networks")) dir.create("results/networks")

# 1. Data Import and Cleaning
#DEGDataSample <- read_csv("DEGDataSample.csv") # change to actual file later
deg <- fread("DEGData.csv") 

deg <- deg %>%
    mutate(
        DISORDER = as.factor(DISORDER),
        TISSUE = as.factor(TISSUE),
        logP = -log10(PVAL),
        Direction = ifelse(LOG2FC > 0 & PVAL < 0.05, "UP", "DOWN")) 

# Remove duplicated gene entries within same disorder-study-tissue
deg <- deg %>%
    group_by(GENEID, DISORDER, STUDY, TISSUE) %>%
    slice_min(PVAL, n = 1) %>%
    ungroup()

# 2. Collapse Across Studies Within Each Disorder
meta_disorder <- deg %>%
    group_by(GENEID, DISORDER) %>%
    summarise(
        mean_log2FC = mean(LOG2FC),
        min_p = min(PVAL),
        combined_p = p.adjust(min(PVAL), method = "BH"),
        Direction = ifelse(mean_log2FC > 0, "UP", "DOWN"),
        .groups = "drop") # . references 

#meta_disorder = core aggregated gene-level dataset 
#collapses study-level differential expression results 
#into a single estimate per gene per disorder.

meta_disorder$GENEID <- sub("\\..*", "", meta_disorder$GENEID) # strip version nums from ensembl id

write.csv(meta_disorder,
          "results/tables/meta_disorder_collapsed.csv",
          row.names = FALSE)

# 3. Gene Annotation
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org") # if SQL not working, specify a different host

annot <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
    filters = "ensembl_gene_id",
    values = unique(meta_disorder$GENEID),
    mart = mart)

annot <- annot %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)

meta_disorder <- left_join(
    meta_disorder,
    annot,
    by = c("GENEID" = "ensembl_gene_id"))

# 4. DEG Summary Statistics

disorder_summary <- meta_disorder %>%
    group_by(DISORDER) %>%
    summarise(
        n_genes = n(),
        n_up = sum(Direction == "UP"),
        n_down = sum(Direction == "DOWN"),
        mean_effect = mean(abs(mean_log2FC)))  # what does mean effect show
# mean effect concerns the average up or downreg

write.csv(disorder_summary,"results/tables/disorder_summary.csv",row.names = FALSE)

# Interpretation dimensions: Regulatory asymmetry (up vs down), Effect size magnitude, Gene count 
#########
# 5. Functional Enrichment Per Disorder
# Evaluate: Synaptic signaling, Immune activation, Chromatin remodeling, Metabolic shifts

run_enrichment <- function(df, disorder_name){
    
    genes <- df %>%
        filter(DISORDER == disorder_name) %>%
        pull(entrezgene_id) %>%
        na.omit()
    
    ego <- enrichGO(
        gene = genes,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        readable = TRUE)
    
    ekegg <- enrichKEGG(gene = genes,organism = "hsa")
    
    return(list(GO = ego, KEGG = ekegg))}

results_list <- lapply(levels(meta_disorder$DISORDER),
                       function(d) run_enrichment(meta_disorder, d))

names(results_list) <- levels(meta_disorder$DISORDER)
#head(results_list)

    # Save GO and KEGG tables
for(d in names(results_list)){
    
    go_res <- results_list[[d]]$GO
    if(is.null(go_res) || nrow(as.data.frame(go_res)) == 0){
        print(paste("Skipping", d, "- no enrichment results"))
        next}
        
    #go_res <- as.data.frame(results_list[[d]]$GO)
    #kegg_res <- as.data.frame(results_list[[d]]$KEGG)
    
    write.csv(go_res,
              paste0("results/enrichment/", d, "_GO.csv"),
              row.names = FALSE)
    
    write.csv(kegg_res,
              paste0("results/enrichment/", d, "_KEGG.csv"),
              row.names = FALSE)
    
    # Save dotplot
    #p_go <- dotplot(results_list[[d]]$GO, showCategory = 10)+ scale_size_area(max_size = 10)
    p_go <- dotplot(go_res, showCategory = 10)+ scale_size_area(max_size = 10)
    ggsave(
        filename = paste0("results/figures/", d, "_GO_dotplot.png"),
        plot = p_go,
        width = 8,
        height = 6,
        dpi = 300)}

#clusterprofiler done already?
#res <- qs::qread(file.path(out_path, paste0(filename_prefix, '_resclusterp_enrichres.qs')))

    # 23. Pathway Enrichment Dotplot Comparison
dotplot(results_list[[1]]$GO, showCategory = 10)+ scale_size_area(max_size = 10)

    # 14. Pathway-Level Cross-Disorder Comparison
go_lists <- lapply(results_list, function(x){x$GO@result$Description})

########

# 6. GSEA Using Ranked Log2FC

run_gsea <- function(disorder_name){
    
    ranked <- meta_disorder %>%
        filter(DISORDER == disorder_name) %>%
        arrange(desc(mean_log2FC)) %>%
        select(entrezgene_id, mean_log2FC) %>%
        na.omit()
    
    geneList <- ranked$mean_log2FC
    names(geneList) <- ranked$entrezgene_id
    
    gsea <- gseGO(
        geneList = geneList,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        keyType = "ENTREZID",
        verbose = FALSE
    )
    
    return(gsea)
}

# gsea viz
gsea_results <- lapply(levels(meta_disorder$DISORDER), run_gsea)
names(gsea_results) <- levels(meta_disorder$DISORDER)

# x-axis: ranked genes by log2FC

# vertical ticks: genes belonging to the pathway

#running score: enrichment signal across the ranked list

# 67Peak deviation indicates strongest enrichment.


p_gsea1 <- gseaplot2(
    gsea_results[["AD"]],
    geneSetID = 1)

ggsave(
    "results/figures/AD_GSEA_enrichment.png",
    p_gsea1,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea2 <- gseaplot2(
    gsea_results[["ADHD"]],
    geneSetID = 1)

ggsave(
    "results/figures/ADHD_GSEA_enrichment.png",
    p_gsea2,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea3 <- gseaplot2(
    gsea_results[["ASD"]],
    geneSetID = 1)

ggsave(
    "results/figures/ASD_GSEA_enrichment.png",
    p_gsea3,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea4 <- gseaplot2(
    gsea_results[["BD"]],
    geneSetID = 1)

ggsave(
    "results/figures/BD_GSEA_enrichment.png",
    p_gsea4,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea5 <- gseaplot2(
    gsea_results[["MDD"]],
    geneSetID = 1)

ggsave(
    "results/figures/MDD_GSEA_enrichment.png",
    p_gsea5,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea6 <- gseaplot2(
    gsea_results[["OCD"]],
    geneSetID = 1)

ggsave(
    "results/figures/OCD_GSEA_enrichment.png",
    p_gsea6,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea7 <- gseaplot2(
    gsea_results[["SZ"]],
    geneSetID = 1)

ggsave(
    "results/figures/SZ_GSEA_enrichment.png",
    p_gsea7,
    width = 8,
    height = 6,
    dpi = 300)

######
# dotplot shows top enriched pathways simultaneously.
p_gsea_dot <- dotplot(
    gsea_results[["AD"]],
    showCategory = 10)

ggsave(
    "results/figures/AD_GSEA_dotplot.png",
    p_gsea_dot,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea_dot2 <- dotplot(
    gsea_results[["ADHD"]],
    showCategory = 10)

ggsave(
    "results/figures/ADHD_GSEA_dotplot.png",
    p_gsea_dot2,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea_dot2 <- dotplot(
    gsea_results[["ASD"]],
    showCategory = 10)

ggsave(
    "results/figures/ASD_GSEA_dotplot.png",
    p_gsea_dot2,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea_dot3 <- dotplot(
    gsea_results[["BD"]],
    showCategory = 10)

ggsave(
    "results/figures/BD_GSEA_dotplot.png",
    p_gsea_dot3,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea_dot4 <- dotplot(
    gsea_results[["MDD"]],
    showCategory = 10)

ggsave(
    "results/figures/MDD_GSEA_dotplot.png",
    p_gsea_dot4,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea_dot5 <- dotplot(
    gsea_results[["OCD"]],
    showCategory = 10)

ggsave(
    "results/figures/OCD_GSEA_dotplot.png",
    p_gsea_dot5,
    width = 8,
    height = 6,
    dpi = 300)
##

p_gsea_dot6 <- dotplot(
    gsea_results[["SZ"]],
    showCategory = 10)

ggsave(
    "results/figures/SZ_GSEA_dotplot.png",
    p_gsea_dot6,
    width = 8,
    height = 6,
    dpi = 300)
##

#####
# ridge plot shows ranked gene distribution contributing to each pathway.
p_gsea_ridge <- ridgeplot(
    gsea_results[["AD"]],
    showCategory = 10)

ggsave(
    "results/figures/AD_GSEA_ridgeplot.png",
    p_gsea_ridge,
    width = 8,
    height = 6,
    dpi = 300)

# 
p_gsea_ridge1 <- ridgeplot(
    gsea_results[["ADHD"]],
    showCategory = 10)

ggsave(
    "results/figures/ADHD_GSEA_ridgeplot.png",
    p_gsea_ridge1,
    width = 8,
    height = 6,
    dpi = 300)
# 
p_gsea_ridge2 <- ridgeplot(
    gsea_results[["ASD"]],
    showCategory = 10)

ggsave(
    "results/figures/ASD_GSEA_ridgeplot.png",
    p_gsea_ridge2,
    width = 8,
    height = 6,
    dpi = 300)

# 
p_gsea_ridge3 <- ridgeplot(
    gsea_results[["BD"]],
    showCategory = 10)

ggsave(
    "results/figures/BD_GSEA_ridgeplot.png",
    p_gsea_ridge3,
    width = 8,
    height = 6,
    dpi = 300)

# 
p_gsea_ridge4 <- ridgeplot(
    gsea_results[["MDD"]],
    showCategory = 10)

ggsave(
    "results/figures/MDD_GSEA_ridgeplot.png",
    p_gsea_ridge4,
    width = 8,
    height = 6,
    dpi = 300)

# 
p_gsea_ridge5 <- ridgeplot(
    gsea_results[["OCD"]],
    showCategory = 10)

ggsave(
    "results/figures/OCD_GSEA_ridgeplot.png",
    p_gsea_ridge,
    width = 8,
    height = 6,
    dpi = 300)

#
p_gsea_ridge6 <- ridgeplot(
    gsea_results[["SZ"]],
    showCategory = 10)

ggsave(
    "results/figures/SZ_GSEA_ridgeplot.png",
    p_gsea_ridge6,
    width = 8,
    height = 6,
    dpi = 300)


# PART 2: PAIRWISE DISORDER ANALYSIS
# 7. Overlap Statistics
presence_matrix <- meta_disorder %>%
    select(GENEID, DISORDER) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = DISORDER,
                values_from = present,
                values_fill = 0)

# all pairwise comparisons
calc_jaccard <- function(a, b){
    intersect = sum(a == 1 & b == 1)
    union = sum(a == 1 | b == 1)
    return(intersect / union)}


# 8. Fisher test NOT WORKING

fisher_pair <- function(d1, d2){
    
    g1 <- meta_disorder %>% filter(DISORDER == d1) %>% pull(GENEID)
    g2 <- meta_disorder %>% filter(DISORDER == d2) %>% pull(GENEID)
    
    # universe is formed from all protein coding genes
    universe <- getBM(
        attributes = "ensembl_gene_id",
        filters = "biotype",
        values = "protein_coding",
        mart = mart
    )$ensembl_gene_id
    
    #a = genes DE in both disorders
    #b = genes DE only in disorder1
    #c = genes DE only in disorder2
    #d = genes not DE in either disorder
        # testing matrix via print
        a <- length(intersect(g1, g2))
        b <- length(setdiff(g1, g2))
        c <- length(setdiff(g2, g1))
        d <- length(setdiff(universe, union(g1, g2)))
        
        mat <- matrix(c(a, b, c, d), nrow = 2)
        
        print(paste("Testing", d1, "vs", d2))
        print(c(a = a, b = b, c = c, d = d))
        print(mat)
    
    ft <- fisher.test(mat)
    
    data.frame(
        overlap = a,
        odds_ratio = ft$estimate,
        log10_p = -log10(ft$p.value))} # take log to prevent extremely small p value

# run fisher for every disorder pair
disorders <- unique(meta_disorder$DISORDER)

pairs <- combn(disorders, 2, simplify = FALSE)

fisher_results <- lapply(pairs, function(pair){
    
    res <- fisher_pair(pair[1], pair[2])
    
    data.frame(
        disorder1 = pair[1],
        disorder2 = pair[2],
        overlap = res$overlap,
        odds_ratio = res$odds_ratio,
        log10_p = res$log10_p)})

fisher_df <- do.call(rbind, fisher_results) # builds fisher df
# Convert back to normal p-values and apply FDR correction.
fisher_df$p_value <- 10^(-fisher_df$log10_p)
fisher_df$FDR <- p.adjust(fisher_df$p_value, method="BH")

write.csv(fisher_df,"results/tables/fisher_pairwise.csv",row.names = FALSE)

# VISUALIZE FISHER (log p value should be 100-150)

# create matrix
disorders <- unique(meta_disorder$DISORDER)

mat <- matrix(NA, length(disorders), length(disorders))
rownames(mat) <- disorders
colnames(mat) <- disorders

for(i in 1:nrow(fisher_df)){
    
    d1 <- fisher_df$disorder1[i]
    d2 <- fisher_df$disorder2[i]
    
    mat[d1, d2] <- fisher_df$log10_p[i]
    mat[d2, d1] <- fisher_df$log10_p[i]}

diag(mat) <- 0

png("results/figures/fisher_overlap_heatmap.png", width=2000, height=2000, res=300)

Heatmap(
    mat,
    name = "-log10(p)",
    cluster_rows = TRUE,
    cluster_columns = TRUE)

dev.off()

# visualize overlap as a network 
edges <- fisher_df

g <- graph_from_data_frame(
    edges[,c("disorder1","disorder2")],
    directed = FALSE)

E(g)$weight <- fisher_df$log10_p

plot(
    g,
    edge.width = E(g)$weight / 20,
    vertex.size = 30)

# 9. Concordance of Directionality
directionality <- meta_disorder %>%
    select(GENEID, DISORDER, Direction) %>%
    pivot_wider(names_from = DISORDER,
                values_from = Direction)

# Compare direction consistency per gene

# 10. Cross-Disorder Network
# Filter genes shared across multiple disorders
shared_genes <- meta_disorder %>%
    group_by(GENEID) %>%
    filter(n_distinct(DISORDER) >= 2) %>%
    ungroup()

edges <- shared_genes %>%
    select(DISORDER, GENEID)

g <- graph_from_data_frame(edges, directed = FALSE)

V(g)$type <- V(g)$name %in% unique(meta_disorder$DISORDER)

V(g)$color <- ifelse(V(g)$type, "red", "gold")
V(g)$size <- ifelse(V(g)$type, 15, 4)

plot(
    g,
    vertex.size = 4,
    vertex.label = NA,
    edge.arrow.size = 0.2)

# Save the plot
png("results/networks/shared_gene_network.png",
    width = 2400,
    height = 2400,
    res = 300)

plot(
    g,
    layout = layout_with_fr(g),
    vertex.label = NA,
    edge.color = "gray70")

dev.off()

disorders <- unique(meta_disorder$DISORDER)

shared_genes <- meta_disorder %>%
    group_by(GENEID) %>%
    filter(n_distinct(DISORDER) >= 2) %>%
    ungroup()

edges <- shared_genes %>% select(DISORDER, GENEID)

g <- graph_from_data_frame(edges, directed = FALSE)

V(g)$type <- V(g)$name %in% disorders
#table(V(g)$type)
stopifnot(is_bipartite(g))

proj <- bipartite_projection(g)

gene_network <- proj$proj2

#optionally remove nodes with no edges
deg <- igraph::degree(gene_network)
gene_network <- delete_vertices(
    gene_network,
    V(gene_network)[deg == 0])
#sum(igraph::degree(gene_network) == 0)

deg <- as.numeric(igraph::degree(gene_network))
wdeg <- as.numeric(igraph::strength(gene_network))
btw <- as.numeric(igraph::betweenness(
    gene_network,
    weights = 1 / E(gene_network)$weight))

centrality <- data.frame(
    gene = V(gene_network)$name,
    degree = deg,
    weighted_degree = wdeg,
    betweenness = btw)

"centrality <- data.frame(
    gene = V(gene_network)$name,
    betweenness = betweenness(gene_network))"

centrality2 <- data.frame(
    gene = V(gene_network)$name,
    degree = degree(gene_network))

# Find hub genes
hub_genes <- centrality[order(-centrality$weighted_degree), ]
head(hub_genes, 20)

# Global Cross-Disorder Integration

# 11. Identify Pan-Disorder Genes

gene_counts <- meta_disorder %>%
    group_by(GENEID) %>%
    summarise(n_disorders = n_distinct(DISORDER))

pan_genes <- gene_counts %>%
    filter(n_disorders == max(n_disorders))

# 12. Global Meta-Effect
global_effect <- meta_disorder %>%
    group_by(GENEID) %>%
    summarise(
        avg_effect = mean(mean_log2FC),
        consistency = sd(mean_log2FC),
        n_disorders = n_distinct(DISORDER)
    )

# 13. Hierarchical Clustering of Disorders
# Reconstructs Transcriptional similarity structure among disorders.
effect_matrix <- meta_disorder %>%
    select(GENEID, DISORDER, mean_log2FC) %>%
    pivot_wider(names_from = DISORDER,
                values_from = mean_log2FC,
                values_fill = 0)

mat <- as.matrix(effect_matrix[,-1])
rownames(mat) <- effect_matrix$GENEID

dist_mat <- dist(t(mat))
hc <- hclust(dist_mat)

plot(hc)


Reduce(intersect, go_lists)

# PART IV: Visualizations
# 15. Volcano Plot Per Disorder
plot_volcano <- function(disorder_name){
    
    df <- meta_disorder %>%
        filter(DISORDER == disorder_name)
    
    p <- ggplot(df, aes(x = mean_log2FC, y = -log10(min_p))) +
        geom_point(alpha = 0.6) +
        geom_vline(xintercept = 0, linetype = "dashed") +
        theme_minimal() +
        labs(
            title = paste("Volcano Plot:", disorder_name),
            x = "Mean log2FC",
            y = "-log10(p)")
    
    return(p)}

# 16. Effect Size Density Comparison Across Disorders
# Detect: Shifted distributions, Broad vs narrow dysregulation patterns
ggplot(meta_disorder, aes(x = mean_log2FC, fill = DISORDER)) +
    geom_density(alpha = 0.4) +
    theme_minimal() +
    labs(
        title = "Effect Size Distributions",
        x = "Mean log2FC",
        y = "Density")

p_density <- ggplot(meta_disorder, aes(x = mean_log2FC, fill = DISORDER)) +
    geom_density(alpha = 0.4) +
    theme_minimal()

ggsave("results/figures/effect_density.png",
       p_density,
       width = 8,
       height = 6,
       dpi = 300)

# 17. Heatmap of Cross-Disorder Gene Effects
# Detects Genetic proximity structure + Disorder grouping patterns

top_genes <- meta_disorder %>%
    group_by(GENEID) %>%
    summarise(total_abs = sum(abs(mean_log2FC))) %>%
    arrange(desc(total_abs)) %>%
    slice(1:100) %>%
    pull(GENEID)

heat_df <- meta_disorder %>%
    filter(GENEID %in% top_genes) %>%
    select(GENEID, DISORDER, mean_log2FC) %>%
    pivot_wider(names_from = DISORDER,
                values_from = mean_log2FC,
                values_fill = 0)

mat <- as.matrix(heat_df[,-1])
rownames(mat) <- heat_df$GENEID

Heatmap(
    mat,
    name = "log2FC",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    cluster_rows = TRUE,
    cluster_columns = TRUE)

# 18. Jaccard Similarity Heatmap Between Disorders

disorders <- levels(meta_disorder$DISORDER)

jmat <- matrix(0, length(disorders), length(disorders))
rownames(jmat) <- disorders
colnames(jmat) <- disorders

for(i in 1:length(disorders)){
    for(j in 1:length(disorders)){
        g1 <- presence_matrix[[disorders[i]]]
        g2 <- presence_matrix[[disorders[j]]]
        jmat[i,j] <- calc_jaccard(g1, g2)}}

Heatmap(
    jmat,
    name = "Jaccard",
    cluster_rows = TRUE,
    cluster_columns = TRUE)

write.csv(jmat,
          "results/tables/jaccard_matrix.csv")

# save heatmap
png("results/figures/jaccard_heatmap.png",
    width = 2000, height = 2000, res = 300)

Heatmap(jmat)

dev.off()


# 19. Network Visualization of Cross-Disorder Hubs
library(ggraph)

centrality_df <- data.frame(
    gene = V(gene_network)$name,
    degree = degree(gene_network)
)

hub_genes <- centrality_df %>%
    arrange(desc(degree)) %>%
    slice(1:50) %>%
    pull(gene)

sub_net <- induced_subgraph(gene_network, hub_genes)

ggraph(sub_net, layout = "fr") +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(size = degree(sub_net))) +
    theme_void()
# 20. Up vs Down Regulation Bar Plot Per Disorder
meta_disorder %>%
    group_by(DISORDER, Direction) %>%
    summarise(n = n()) %>%
    ggplot(aes(x = DISORDER, y = n, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
        title = "Regulatory Polarity by Disorder",
        y = "Number of Genes"
    )

p_bar <- meta_disorder %>%
    group_by(DISORDER, Direction) %>%
    summarise(n = n()) %>%
    ggplot(aes(x = DISORDER, y = n, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal()

ggsave("results/figures/up_down_barplot.png",
       p_bar,
       width = 8,
       height = 6,
       dpi = 600)

# 21. PCA of Disorder Effect Profiles
pca <- prcomp(t(mat), scale. = TRUE)

pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Disorder = rownames(pca$x)
)

ggplot(pca_df, aes(x = PC1, y = PC2, label = Disorder)) +
    geom_point(size = 4) +
    geom_text(vjust = -1) +
    theme_minimal() +
    labs(title = "PCA of Disorder Gene Profiles")


# 22. Circos Plot of Shared Genes
shared_counts <- combn(disorders, 2, function(pair){
    g1 <- meta_disorder %>% filter(DISORDER == pair[1]) %>% pull(GENEID)
    g2 <- meta_disorder %>% filter(DISORDER == pair[2]) %>% pull(GENEID)
    data.frame(
        from = pair[1],
        to = pair[2],
        value = length(intersect(g1, g2))
    )
}, simplify = FALSE)

shared_df <- bind_rows(shared_counts)

chordDiagram(shared_df)

# SAVE PAIRWASE EFFECT SCATTER PLOTS
for(i in 1:(length(disorders)-1)){
    for(j in (i+1):length(disorders)){
        
        d1 <- disorders[i]
        d2 <- disorders[j]
        
        p_pair <- compare_pair(d1, d2)
        
        ggsave(
            filename = paste0("results/figures/", d1, "_vs_", d2, ".png"),
            plot = p_pair,
            width = 6,
            height = 6,
            dpi = 300
        )
    }
}



# 24. Gene Effect Consistency Scatter Between Pairs
compare_pair <- function(d1, d2){
    
    df1 <- meta_disorder %>%
        filter(DISORDER == d1) %>%
        select(GENEID, mean_log2FC)
    
    df2 <- meta_disorder %>%
        filter(DISORDER == d2) %>%
        select(GENEID, mean_log2FC)
    
    merged <- inner_join(df1, df2, by = "GENEID")
    
    ggplot(merged, aes(x = mean_log2FC.x, y = mean_log2FC.y)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", se = FALSE) +
        theme_minimal() +
        labs(
            title = paste(d1, "vs", d2),
            x = d1,
            y = d2
        )
}

compare_pair("BD", "SCZ")
compare_pair("BD", "ASD")
compare_pair("BD", "AD")
compare_pair("BD", "MDD")

##############
#SAVE BLOCKS BELOW 
# assumed the following objects exist: meta_disorder, disorder_summary, 
#jmat, hc, pca_df, gene_network

# Automatically Save All Disorders Volcano Plots
for(d in levels(meta_disorder$DISORDER)){
    p <- plot_volcano(d)
    
    ggsave(
        filename = paste0("results/figures/", d, "_volcano.png"),
        plot = p,
        width = 8,
        height = 6,
        dpi = 600
    )
}

# Save Heatmap of Top Genes
png("results/figures/global_effect_heatmap.png",
    width = 2500, height = 2500, res = 300)

Heatmap(mat)

dev.off()

png("heatmap.png", width = 2000, height = 2000, res = 300)
Heatmap(mat)
dev.off()

# Save PCA Plot

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Disorder)) +
    geom_point(size = 4) +
    geom_text(vjust = -1) +
    theme_minimal()

ggsave("results/figures/pca_disorders.png",
       p_pca,
       width = 6,
       height = 5,
       dpi = 300)

# Save Hierarchical Clustering
png("results/figures/disorder_dendrogram.png",
    width = 2000, height = 2000, res = 300)

plot(hc)

dev.off()

# Save Gene Network as Image

p_net <- ggraph(sub_net, layout = "fr") +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(size = degree(sub_net))) +
    theme_void()

ggsave("results/networks/gene_network.png",
       p_net,
       width = 8,
       height = 6,
       dpi = 300)


# Export Network Edge List
# Allows downstream use in Cytoscape.
edge_list <- as_data_frame(gene_network)

write.csv(edge_list,
          "results/networks/gene_network_edges.csv",
          row.names = FALSE)

# Save Circos Plot
png("results/figures/circos_overlap.png",
    width = 2000, height = 2000, res = 300)

chordDiagram(shared_df)

dev.off()

# Save Pan-Disorder Genes
write.csv(pan_genes,
          "results/tables/pan_disorder_genes.csv",
          row.names = FALSE)

# Save Global Effect Table
write.csv(global_effect,
          "results/tables/global_meta_effects.csv",
          row.names = FALSE)
dim(global_effect)

