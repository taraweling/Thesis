# packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"), ask = FALSE)

library(GEOquery)
library(limma)
library(dplyr)

# load data
gset <- getGEO("GSE61672", GSEMatrix = TRUE)[[1]]

expr <- exprs(gset)
pheno <- pData(gset)

# enforce matrix
expr <- as.matrix(expr)
mode(expr) <- "numeric"

# ensure alignment
pheno <- pheno[colnames(expr), ]

# inspect phenotype structure (run once, then adapt parsing if needed)
colnames(pheno)

# ---- PARSE phenotype fields ----
# adjust patterns if your GEO fields differ

pheno$anxiety_status <- ifelse(
    grepl("anxiety", pheno$characteristics_ch1, ignore.case = TRUE),
    1, 0
)

pheno$age <- as.numeric(sub(".*age[:= ]", "", pheno$characteristics_ch1))

pheno$sex <- factor(
    ifelse(grepl("female", pheno$characteristics_ch1, ignore.case = TRUE),
           "F", "M")
)

pheno$race <- factor(
    sub(".*race[:= ]", "", pheno$characteristics_ch1)
)

# ---- FILTER before modeling ----
keep <- rowMeans(expr) > 5
expr <- expr[keep, ]

# ---- DESIGN MATRIX ----
design_df <- data.frame(
    group = factor(pheno$anxiety_status),
    age   = pheno$age,
    sex   = pheno$sex,
    race  = pheno$race
)

design <- model.matrix(~ group + age + sex + race, data = design_df)

# ---- DIFFERENTIAL EXPRESSION ----
fit <- lmFit(expr, design)
fit <- eBayes(fit)

deg <- topTable(
    fit,
    coef = "group1",
    number = Inf,
    adjust.method = "fdr"
)

sig <- deg %>%
    filter(adj.P.Val < 0.1)

# ---- MAP PROBES TO GENES ----
fdat <- fData(gset)

deg$gene <- fdat$GeneSymbol[
    match(rownames(deg), rownames(fdat))
]

sig_genes <- deg %>%
    filter(adj.P.Val < 0.1, !is.na(gene), gene != "") %>%
    group_by(gene) %>%
    slice_max(order_by = abs(logFC), n = 1)

# ---- OUTPUT CHECKS ----
nrow(sig)         # number of significant probes
nrow(sig_genes)   # number of unique genes
head(sig_genes)

"# use identical regression but without subsetting by sex
group <- factor(pheno$anxiety_status)

design <- model.matrix(
    ~ group + pheno$age + pheno$sex + pheno$race
)

# run differential expression
fit <- lmFit(expr, design)
fit <- eBayes(fit)

deg <- topTable(
    fit,
    coef = 'group1',
    number = Inf,
    adjust.method = 'fdr'
)

# apply the study threshold
sig <- deg %>%
    filter(adj.P.Val < 0.1)

nrow(sig)

# convert probe IDs to genes
deg$gene <- fData(gset)$GeneSymbol[match(rownames(deg), rownames(fData(gset)))]

sig_genes <- sig %>%
    group_by(gene) %>%
    slice_max(order_by = abs(logFC), n = 1)

# remove low expression probes
keep <- rowMeans(expr) > 5
expr <- expr[keep, ]

design <- model.matrix(~ pheno$anxiety_status)
fit <- lmFit(expr, design)
fit <- eBayes(fit)

deg <- topTable(fit, coef=2, number=Inf, adjust="fdr")
sig <- deg[deg$adj.P.Val < 0.1, ]

"
deg$significant <- deg$adj.P.Val < 0.1


ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(
        title = "Volcano Plot",
        x = "log2 Fold Change",
        y = "-log10(FDR)"
    )

library(pheatmap)

# match significant genes back to expression matrix
sig_ids <- rownames(sig)
expr_sig <- expr[sig_ids, ]

# scale genes (row-wise z-score)
expr_scaled <- t(scale(t(expr_sig)))

annotation <- data.frame(
    Anxiety = pheno$anxiety_status
)
rownames(annotation) <- colnames(expr_scaled)

pheatmap(
    expr_scaled,
    annotation_col = annotation,
    show_rownames = FALSE,
    show_colnames = FALSE
)

pca <- prcomp(t(expr), scale. = TRUE)

pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Anxiety = pheno$anxiety_status
)

ggplot(pca_df, aes(PC1, PC2, color = Anxiety)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "PCA of Expression Data")

top_genes <- sig_genes %>%
    arrange(desc(abs(logFC))) %>%
    head(20)

ggplot(top_genes, aes(x = reorder(gene, logFC), y = logFC)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(
        title = "Top DEGs",
        x = "Gene",
        y = "log2 Fold Change"
    )

