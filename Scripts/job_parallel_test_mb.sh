#!/bin/bash

# ---------------------------------------------------------
# TUM - Technichal University of Munich
#
# Authors:  Aydin Uzun
# Date: 2022
# Purpose:  Run the R script on SLURM calculating the ppb errors according to Taylor's framwork 
# ---------------------------------------------------------

# change this part
# set the name of batch request depending on the instrument name
#SBATCH -J parallel_mb
# set the directory depending on the date
#SBATCH -D /dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_20180905/by-id
# set the output logs depending on the date
#SBATCH -o /dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_20180905/by-id/logs/myjob.%j.%N.out
# change this part

#SBATCH --export=NONE
#SBATCH --get-user-env	
# use cm2_tiny and only 1 node
#SBATCH --clusters=cm2_tiny 
#SBATCH --partition=cm2_tiny
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=aydin.uzun@tum.de
#SBATCH --time=03:00:00

# load correct dependencies to compile R libraries
module unload intel-mpi 
module unload intel
module load gcc/11.2.0
# tools for required packages
# The required packages are raster, geosphere, MASS, ncdf4
module load geos/3.9.1-gcc11
module load proj/8.1.0-gcc11
module load sqlite/3.36.0-gcc11
module load netcdf-hdf5-all/4.7_hdf5-1.10-intel21-serial
module load gdal/3.3.3
module load pkgconf/1.8.0
module load slurm_setup
# load the R module
module load r/4.1.2-gcc11-mkl

# R script to run - change this part depending on the instrument name
Rscript parallel_test_mb.R