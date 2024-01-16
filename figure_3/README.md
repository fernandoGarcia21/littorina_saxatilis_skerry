# Scripts to generate the Figure 3 in Garcia et al. 2024

This directory contains the scripts needed to generate the Figure 3 in our MS, currently available on BioRxiv as a preprint:

**Predicting rapid adaptation in time from adaptation in space: a 30-year field experiment in marine snails**

doi: https://doi.org/10.1101/2023.09.27.559715

Clarifying notes:
- 2_Main_Fig3.R is the main script that generates the Figure 3. This script imports other two scripts and call their functions. To generate Figure 3 just the 2_Main_Fig3.R needs to be run and make sure all other scripts are available in the specified paths.
- The input files with the trajectories of arrangement frequencies must be previously generated through the script 1_InversionsGenerateTrajectoriesData.R and must be available in the Data directory. 
- Likewise, the 500K random draws from the likelihood surface must be available in the Data folder too: Interpolation demographic parameters_random values 23 Sept 2023.csv.