# Runs seqsnptm003_demographic inference.R for a grid of parameter combinations.
#
#author: "Anja Westram"
#
######################################################################################################

#!/bin/bash -ve

# request memory for job (default 6G, max 72G)
#$ -l mem=6G
#$ -l rmem=6G
# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)
#$ -l h_rt=39:59:00
# current environment settings are used for the job???
#$ -V
#$ -t 1-1875
##$ -tc 20

taskid=${SGE_TASK_ID}

a=`head -n $taskid SEQSNPTM002_GRID.txt |tail -n 1`

N0=`echo $a | cut -f 1 -d' '`
r=`echo $a | cut -f 2 -d' '`
K=`echo $a | cut -f 3 -d' '`
M=`echo $a | cut -f 4 -d' '`
gen_factor=`echo $a | cut -f 5 -d' '`

module add apps/R/3.3.2/gcc-4.8.5
Rscript --vanilla ./seqsnptm003_TM_new.R $N0 $r $K $M $gen_factor