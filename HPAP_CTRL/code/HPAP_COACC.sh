#!/bin/bash

# set the number of nodes (each with several cpus/cores)
#SBATCH --export=ALL
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --mem=320G

# Output files
#SBATCH -o COACC_HPAP_CTRL.out # STDOUT
#SBATCH -e COACC_HPAP_CTRL.err # STDERR

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=georgina.fuentes@upf.edu

# load necessary modules with `module load XXX`
module load Singularity/3.8.7-GCCcore-13.2.0

# write the code you will run in the clusster
singularity exec /scratch/lab_lpasquali/shared_data/rstudio-singularity/bioconductor_R-4.4.2_RELEASE_3_20_MACS2.sif Rscript /homes/users/gfuentes/scratch/projects/spicey_old/HPAP_CTRL/code/COACC_HPAP.R
