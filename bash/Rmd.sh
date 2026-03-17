#!/bin/bash -l
#SBATCH --time=48:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=150g
#SBATCH --tmp=150g
#SBATCH --mail-type=END
#SBATCH --mail-user=cbrault@umn.edu
#SBATCH --output=GP.pred.RKHS.CV00
#SBATCH --job-name=GP.pred.RKHS.CV00


conda activate r_env
export RSTUDIO_PANDOC=/users/7/cbrault/miniconda/envs/r_env/bin/
#export R_LIBS_USER="/users/7/cbrault/R/x86_64-pc-linux-gnu-library/4.4/"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

#R -e "rmarkdown::render('/users/7/cbrault/Predictathon/analysis/combine_GRM.Rmd')"
R -e "rmarkdown::render('/users/7/cbrault/Predictathon/analysis/genom_pred.Rmd')"

