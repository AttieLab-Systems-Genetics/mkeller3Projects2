# R Analysis Scripts

This folder contains a collection of R scripts and RMarkdown files for processing, analyzing, and visualizing multi-omic data (lipids, metabolites, RNA) in Diversity Outbred (DO) mice. The scripts are organized into functional pipelines below.

## How I created this file and repo

This repo mirrors files on the [ResearchDrive](https://it.wisc.edu/services/researchdrive/)
in folder `mkeller3/General/Projects2/R scripts`.
This repo is private, as is the ResearchDrive (RD) area.
I created this `README` document using a simple prompt to
"create a `README.md` document that concisely outlines the contents of this folder".
This section was added afterward.

I connected the `R scripts` folder on RD with this GitHub repository (or "repo")
using help from
[Get started with GitHub](https://happygitwithr.com/usage-intro#usage-intro).
Specifically, I used instructions about
[Existing project, GitHub last](https://happygitwithr.com/existing-github-last#existing-github-last)
to set up the connection.

## Pipelines

- **[Data Processing & Normalization](#1-data-processing--normalization)**: Tools for cleaning, batch correction, and Rank-Z transformations.
- **[QTL Mapping & Hotspot Analysis](#2-qtl-mapping--hotspot-analysis)**: Pipelines for genetic association mapping and hotspot detection.
- **[Network Analysis (WGCNA)](#3-network-analysis-wgcna)**: Core scripts for co-expression network construction and module analysis.
- **[Mediation Analysis](#4-mediation-analysis)**: Bayesian causal modeling tools for multi-omic data.
- **[Visualization](#5-visualization)**: Specialized plotting functions for summaries and heatmap visualization.
- **[Utilities & Development](#6-utilities--development)**: Shared helper functions and prototyping scripts.

### 1. Data Processing & Normalization

Scripts for cleaning, batch correction, and statistical transformations.

- `Assessment of batch effects and ComBat correction.R`: Comprehensive tool for assessing batch influence and applying ComBat normalization.
- `AUC_Batch_Churchill_loop.R`: Batch adjustment for Area Under the Curve (AUC) data using a random effects model (LMM).
- `AUC_Batch_Churchill_loop_13C.R`: Specialized version of the batch adjustment for 13C isotope data.
- `combat_correction_FAHFA_lipids_MPK.R`: Targetted ComBat correction for FAHFA lipids.
- `Script to rankz FAHFA lipids in DO mice.R`: Rank-Z transformation for lipid data to ensure normality.
- `Script to rankz metabolite AUC values.R`: Rank-Z transformation for metabolite AUC values.
- `Script to rankz transform RNA data for WGCNA.R`: Rank-Z transformation prepared specifically for downstream WGCNA analysis.

### 2. QTL Mapping & Hotspot Analysis

Tools for mapping quantitative trait loci (QTL) and identifying enrichment hotspots.

- `Standalone script for analyzing diet and sex specific liver metabolite QTL.R`: An integrated pipeline for metabolite QTL mapping with diet and sex considerations.
- `Mapping FAHFA lipids DO adipose rz transformed ComBat corrected values.Rmd`: Rmd notebook for mapping lipids from specialized transformed/corrected data.
- `Manhattan plot and identification of hotspots FAHFA QTL DO500 adipose.R`: Identification of QTL hotspots and generation of Manhattan plots for FAHFA traits.
- `Trans eQTL hotspot analysis for sex and diet split eQTL summary files.R`: Analysis of trans-eQTL hotspots across different sex/diet subsets.
- `Extracting QTL scans from RDS files for plotting SNP association FAHFA DO500 adipose.Rmd`: Utilities for extracting and visualizing SNP-level associations from saved QTL objects.

### 3. Network Analysis (WGCNA)

Core scripts for Weighted Gene Co-expression Network Analysis.

- `Script to run modified WGCNA modules.R`: A robust, modified version of the WGCNA pipeline for identifying modules in metabolomic or transcriptomic data.
- `module_gsea_enrichment.R`: Tools for functional enrichment and GSEA on identified network modules.
- `corrplot for ME ME trait correlation and pvalues.R`: Post-WGCNA analysis to correlate Module Eigengenes (MEs) with clinical traits.
- `corrplot for trait trait correlations with Pval and ready input file for Cytoscape network construction.R`: Generates correlation networks and prepares export files for Cytoscape.

### 4. Mediation Analysis

Bayesian modeling and interactive tools for causal inference.

- `Additive_mediator_script_for_shiny_v3.1_NewDO.R`: Core wrapper for Bayesian mediation (bmediatR), designed to support interactive Shiny applications.
- `Additive_mediator_script_for_shiny_v3_NewDO.R`: Previous version of the mediation wrapper.
- `source_data_loading_and_functionset_v2.R`: Essential helper functions and data loading logic for the mediation tools.

### 5. Visualization

Specialized scripts for generating high-quality analysis plots.

- `Standalone retangular Manhattan plot Grcm39.R`: Generates standard Manhattan plots using GRCm39 coordinates.
- `qtl_heatmap_v6_SM_ME_QTL.R`: Visualizes multi-trait QTL summaries via heatmaps.
- `script to draw heatmap for ME Add QTL 500DO SM sex specific MEs.R`: Targetted heatmap for sex-specific module eigengene QTL.
- `script to draw heatmap for diet dependent clinical traits.R`: Heatmap visualization for diet effects.
- `ggplots_boxwhisker_gridwork_RNA_modules_Ex1_Lamming.R` & `...Ex3_Lamming.R`: Grid-based boxplots for RNA modules.

### 6. Utilities & Development

- `prototype_running_new_fxn.R`: Incremental prototypes for testing new analysis functions.
- `Shiny app prototype for Mafa SNP DEG analysis.R`: Integrated genomic explorer for MafA SNPs and DEGs, featuring Manhattan plots and locus-specific schematics.
- `source_data_loading_and_functionset_v2.R`: Shared utility functions used across multiple pipelines.
