
# Script that calls all other scripts 
# Scripts must be called in this order

source('R/required_funs_libs.R')

## Dropseq
# Merge data into single DGE and get qc plots
source('R/drop_merge_qc.R')

# Get cell-types and single cell trajectories 
source('R/drop_celltypes_trajectory.R')

# Depending on your RAM, you may have to run rm(list=ls()) before starting the DroNc-seq analysis
# none of the results for above are needed below

## DroNcseq
# Merge data into single DGE and get qc plots
source('R/dronc_merge_qc.R')

# Run Seurat on merged Drop DGE
source('R/dronc_celltypes_trajectory.R')


