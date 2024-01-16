# Scripts to test the hipothesis predictions at phenotypes, SNPs, and inversion level in Garcia et al. 2024

This directory contains the scripts needed to test the accuracy of the predictions in our MS, currently available on BioRxiv as a preprint:

**Predicting rapid adaptation in time from adaptation in space: a 30-year field experiment in marine snails**

doi: https://doi.org/10.1101/2023.09.27.559715

Clarifying notes:
- 2_Main_Predictions_Tests.R is the main script that performs Chi-square tests, Student t-tests, and Fisher's exact test as shown in Tables S9 to S12. This script imports the GeneralsSkerry.R script and call their functions. To run the hypothesis tests, just the 2_Main_Predictions_Tests.R needs to be run and make sure the GeneralsSkerry.R is available in the specified path.
- To be able to test the hypothesis for inversions, the input files with karyotypes ids for complex inversions and karyotypes of all individuals for all inversions must be previously generated through the script 1_Inversions_Individuals_Karyotypes.R and must be available in the Data directory.