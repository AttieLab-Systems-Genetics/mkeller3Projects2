\# Run these lines in R console

\# 1. Load all required libraries

library(shiny)

library(data.table)

library(ggplot2)

library(plotly)

library(GenomicScores)

library(phastCons35way.UCSC.mm39)

library(GenomicRanges)

\# 2. Source the QTL Mapper function from the shared drive

source(\"W:/General/Projects2/R scripts/QTL_Mapper_Functions.R\")

\# 3. Launch the Shiny App

launch_qtl_app()
