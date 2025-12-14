# Coral-Skeleton-Endoliths-Project
R scripts for analyzing endolithic microbial communities in coral skeletons (Complex vs. Robust clades).
# Analysis of Endolithic Microbial Communities in Coral Skeletons

This repository contains the R code used for statistical analysis and visualization in the paper:
**"Host evolutionary lineage shapes the assembly and network topology of endolithic bacterial and archaeal communities in scleractinian coral skeletons"**

## Contents
- `01_Alpha_Beta_Diversity.R`: Code for rarefaction curves, PCoA, and boxplots (Figures 2 & S1-S3).
- `02_Composition_LEfSe.R`: Code for alluvial plots and preparing data for LEfSe (Figure 3 & S4).
- `03_Network_Analysis.R`: Code for co-occurrence network construction and topological role analysis (Figure 4 & S5).
- `04_Functional_Prediction.R`: Code for analyzing PICRUSt2 output (Figure 5).

## Dependencies
The analysis was performed using R (v4.3.0). Key packages include:
- phyloseq
- vegan
- ggplot2
- igraph
- ggtree

## Data Availability
The raw sequencing data and metadata associated with this code are available on Figshare: [Insert your Figshare Link Here]
