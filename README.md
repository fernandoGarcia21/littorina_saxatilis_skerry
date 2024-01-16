# Input data and scripts for the analysis in Garcia et al. 2024

This directory contains the datasets and scripts needed to conduct the analyses in our MS, currently available on BioRxiv as a preprint:

**Predicting rapid adaptation in time from adaptation in space: a 30-year field experiment in marine snails**

doi: https://doi.org/10.1101/2023.09.27.559715

Clarifying notes:

- The Data folder contains all input datasets that are needed to run the scripts.

- The scripts were organized in directories according to the arrangement of the figures/analysis in the manuscript:
  - The figure_1, figure_2, and figure_3 directories contain the scripts for the three figures in the main text. 
  - The supplementary_figures direcotory contains the scripts to generate the figures in the supplementary materials and methods: Figure S3, S4, S5, S6, S11, S12, S14, S15, S16, and S17. 
  - The demographic_inference directory contains the scripts for the downstream analysis of the demographic inference, the estimation of the expected range of allele frequency change without selection, as well as the Mathematica workbook of the likelihood surface of demographic parameters (Skerry interpolation 9.23 v2.nb).
  - The generals directory has one script with general functions that are invoked from multiple other scripts.
  - The prediction_tests directory has the scripts for Chi-square tests, Student t-tests, and Fisher's exact test that are requiered to generate Tables S9, S10, S11, and S12.
  - The selection_phenotypic_analysis directory contain the scripts to estimate the Strength of selection based on phenotypes.
  