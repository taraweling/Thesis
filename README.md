# Cross-Disorder GRN Thesis Pipeline

This repository analyzes cross-disorder differential expression in the context of brain gene regulatory networks (GRNs).  
It combines:

1. Tissue-specific GRN matrices (PANDA/GRAND-derived; bipartite TF -> target edges)
2. Cross-study DEG tables for psychiatric/neurologic disorders
3. Python graph filtering/scoring/visualization
4. R-based enrichment, overlap, and meta-disorder analyses

## Repository Layout

- `src/run.py`: main Python pipeline for GRN + DEG integration.
- `src/graph_utils.py`: data loading, Ensembl mapping, GRN filtering/merging helpers.
- `src/graph_algos.py`: regulatory scoring, summary stats, overlap metrics, CSV writers.
- `src/graph_viz.py`: PyVis/Matplotlib graph visualizations.
- `src/degensembl.py`: helper script to fill missing Ensembl IDs in DEG tables.
- `src/thesisp1.R`: primary R analysis pipeline (enrichment, GSEA, overlap tests, networks).
- `src/p1thesis.R`: older/alternate thesis pipeline version.
- `src/p1permutationtest.R`: standalone permutation overlap test script.
- `src/p1deganalysis.R`, `src/p1geneontology.R`, `src/p1venndiagrams.R`, `src/wingodatasetanalysis.R`: exploratory/auxiliary analysis scripts.
- `src/data/`: input data files.
- `src/results/`: generated figures, tables, and HTML network visualizations.
- `src/lib/`: vendored JS/CSS assets used by HTML network outputs.
- `src/run_merge.sbatch`: SLURM wrapper to run `run.py`.

## Data Inputs

### GRN Matrices (CSV)
Used by `run.py` with edge thresholding:

- `src/data/Brain_Other.csv`
- `src/data/Brain_Basal_Ganglia.csv`
- `src/data/Brain_Cerebellum.csv`

Format: first column = TF, remaining columns = target genes, cells = edge weights.

### DEG Tables (CSV)
Main Python default:

- `src/data/DEGDataStrictLFC.csv`

Other variants included:

- `DEGData.csv`, `DEGDataStrictestLFC.csv`, `DEGDataSample.csv`

Expected columns:

- `DISORDER, STUDY, YEAR, TISSUE, GENEID, LOG2FC, PVAL` (some files also include `-LOGPVAL`)

## Python Pipeline (`src/run.py`)

Run from inside `src/`:

```bash
cd src
python run.py
```

What it does:

1. Reads brain GRN matrices and filters edges by weight (`threshold = 1.7`).
2. Merges tissue GRNs with averaged duplicate TF->gene edge weights.
3. Converts IDs to Ensembl where needed (Ensembl REST API).
4. Loads disorder DEG lists for:
   - `AD, ADHD, ASD, BD, MDD, OCD, SZ`
5. Builds per-disorder network variants:
   - DEG TF only (`tf_grn`)
   - DEG targets only (`deg_grn`)
   - DEG TF + DEG targets (`detf_deggrn`, when distinct)
6. Computes regulator scores, edge/log2FC summaries, and pairwise Jaccard overlaps.
7. Writes per-disorder and cross-disorder outputs to `src/results/`.

Key Python outputs:

- `src/results/deggrn_disorder_summary.csv`
- `src/results/deggrn_overlap_summary.csv`
- `src/results/deggrn_jaccard_*.csv`
- `src/results/overlay_tf_grn_disorders.html`
- `src/results/pyviz_*_deggrn.html` and `src/results/deggrn_*.html`

## R Pipeline (`src/thesisp1.R`)

Run from inside `src/`:

```bash
cd src
Rscript thesisp1.R
```

High-level steps:

1. Load and clean DEG table (`DEGData.csv`)
2. Collapse to per-gene/per-disorder meta effects
3. Annotate with HGNC/Entrez (biomaRt)
4. GO/KEGG enrichment and GSEA
5. Pairwise overlap statistics (Fisher/permutation variants)
6. Cross-disorder visualizations and network hub analyses
7. Export tables/figures under `src/results/tables`, `src/results/figures`, `src/results/networks`, `src/results/enrichment`

## Environment Requirements

### Python

Minimum packages used by core pipeline:

- `requests`
- `numpy`
- `matplotlib`
- `networkx`
- `pyvis`
- `urllib3`

### R

Scripts use many packages, including:

- CRAN: `tidyverse`, `dplyr`, `data.table`, `ggplot2`, `igraph`, `ggraph`, `pheatmap`, `circlize`, `RColorBrewer`, `ggupset`
- Bioconductor: `biomaRt`, `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `ComplexHeatmap`, `DOSE`, `GEOquery`, `limma`

## Notes

- `run.py` assumes execution from `src/` (relative paths like `data/...` and `results/...`).
- Ensembl ID conversion in Python requires internet access to Ensembl REST.
- R annotation/enrichment steps require internet access for some biomaRt/Bioc calls.
- `src/p1thesis.R` and other `p1*.R` scripts are older/experimental variants; `thesisp1.R` is the most complete end-to-end R workflow in this repo.
