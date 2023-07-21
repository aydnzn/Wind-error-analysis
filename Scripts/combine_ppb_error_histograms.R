# ---------------------------------------------------------
# TUM - Technichal University of Munich
#
# Authors:  Aydin Uzun
# Date: 2022
# Purpose:  Create daily ppb error histogram - This script could be run on the cluster. 
# Run: Run this script on the cluster - or you can do it on your local PC if you have the histogram matrices
#      and scaling factors. 
#      For cluster:
#                 1. Load correct dependencies to compile R libraries
#                           module unload intel-mpi
#                           module unload intel
#                           module load gcc
#                 2. Load your desired R module 
#                           module load r/3.6.3-gcc11-mkl  
#                 3. Run your script
#                           Rscript combine_ppb_error_histograms.R
# ---------------------------------------------------------
rm(list = ls())
############################ Day under investigation HAM_XXXXXXX 
# The working directory is "/dss/dsstumfs01/pn69ki/pn69ki-dss-0000/STILT_Ham/HAM_large/HAM_XXXXXXXX/by-id
# Upload this script to the working directory and submit the job 
# XXXXXXXX being e.g. 20180905
Date <- "20210905"
# Please note that you need to include the scaling factors for the corresponding day in your working directory
# The scaling factors for each day can be found on '/Volumes/esm/11-Thesis/
# 03-Scientific-Internship/2021 FP Aydin Uzun/Data/Scaling_factors'
# scales <- read.csv(file = 'HAM_XXXXXXXX_scaling.csv') e.g. 20210905
scales <- read.csv(file = paste("HAM_",Date,"_scaling.csv",sep=""))
############################ Day under investigation HAM_XXXXXXX 



# Subscripts 1,2,3 correspond to the different inventories used
# e.g. mb_hist_ppb_errors_1 contains the histogram counts of the bins created 
# after the error simulations are done using the particle trajectory files for this specific day 
# and for this specific instrument location
# So there are different distributions corresponding to different instruments, different datetimes and different particle release heights
# The goal is to combine them to create a daily distribution
mb_hist_ppb_errors_1 <- readRDS("parallel_mb_hist_ppb_errors_1_counts.rds") # int [1:20,1:507] 
# 507 is the number of different combinations of timepoints and particle release heights
# there are 39 timepoints and 13 particle release heights
mc_hist_ppb_errors_1 <- readRDS("parallel_mc_hist_ppb_errors_1_counts.rds") # int [1:20,1:507]
md_hist_ppb_errors_1 <- readRDS("parallel_md_hist_ppb_errors_1_counts.rds") # int [1:20,1:507]
me_hist_ppb_errors_1 <- readRDS("parallel_me_hist_ppb_errors_1_counts.rds") # int [1:20,1:507]

mb_hist_ppb_errors_2 <- readRDS("parallel_mb_hist_ppb_errors_2_counts.rds") # int [1:20,1:507]
mc_hist_ppb_errors_2 <- readRDS("parallel_mc_hist_ppb_errors_2_counts.rds") # int [1:20,1:507]
md_hist_ppb_errors_2 <- readRDS("parallel_md_hist_ppb_errors_2_counts.rds") # int [1:20,1:507]
me_hist_ppb_errors_2 <- readRDS("parallel_me_hist_ppb_errors_2_counts.rds") # int [1:20,1:507]

mb_hist_ppb_errors_3 <- readRDS("parallel_mb_hist_ppb_errors_3_counts.rds") # int [1:20,1:507]
mc_hist_ppb_errors_3 <- readRDS("parallel_mc_hist_ppb_errors_3_counts.rds") # int [1:20,1:507]
md_hist_ppb_errors_3 <- readRDS("parallel_md_hist_ppb_errors_3_counts.rds") # int [1:20,1:507]
me_hist_ppb_errors_3 <- readRDS("parallel_me_hist_ppb_errors_3_counts.rds") # int [1:20,1:507]

vert_pres <-scales$x # vertical scaling factors for the 13 atm. layers
# please remember that these are in the increasing order
# first element corresponds to first atm. layer [20m], last element corresponds to the 13.th atm layer [2220m]

# Read the names of the trajectories related to "mb"
file_names_mb <- readRDS("parallel_mb_particle_files.rds") # chr [1:507]
# e.g. mb_202109050630_1060m" "mb_202109050630_1244m" "mb_202109050630_1432m" .....
file_names_mb_list <- strsplit(file_names_mb, split = "_") # list of 507
# e.g. chr[1:3] "mb" "202109050630" "1060m"
# chr[1:3] "mb" "202109050630" "1244m" .....
file_names_mb_matrix <- matrix(unlist(file_names_mb_list),nrow=length(file_names_mb_list),ncol=3,byrow=TRUE)
# e.g.
#   V1 V2 V3
# 1 mb 202109050630 1060m
# 2 mb 202109050630 1244m
# 3 mb 202109050630 1432m
# .....
# 507 .............
file_dates_mb_matrix<-file_names_mb_matrix[,2:3]
# e.g.
#      [,1]           [,2]   
#[1,] "202109050630" "1060m"
#[2,] "202109050630" "1244m"
#[3,] "202109050630" "1432m"
#[4,] "202109050630" "1623m"
#[5,] "202109050630" "1818m"
#[6,] "202109050630" "186m" 
#[7,] "202109050630" "2017m"
#[8,] "202109050630" "20m"  
#[9,] "202109050630" "2220m"
#[10,] "202109050630" "355m" 
#[11,] "202109050630" "527m" 
#[12,] "202109050630" "701m" 
#[13,] "202109050630" "879m" 
#[14,] "202109050645" "1060m"
#[15,] "202109050645" "1244m"
#[16,] "202109050645" "1432m"
#[17,] "202109050645" "1623m"
#[18,] "202109050645" "1818m"
#[19,] "202109050645" "186m" 
#[20,] "202109050645" "2017m"
# .....
# It's noticed that the order of the particle release heights are the same
# always start with 1060m and end with 879m then next datetime

# Now do the same for the other instruments
file_names_mc <- readRDS("parallel_mc_particle_files.rds")
file_names_mc_list <- strsplit(file_names_mc, split = "_")
file_names_mc_matrix <- matrix(unlist(file_names_mc_list),nrow=length(file_names_mc_list),ncol=3,byrow=TRUE)
file_dates_mc_matrix<-file_names_mc_matrix[,2:3]

file_names_md <- readRDS("parallel_md_particle_files.rds")
file_names_md_list <- strsplit(file_names_md, split = "_")
file_names_md_matrix <- matrix(unlist(file_names_md_list),nrow=length(file_names_md_list),ncol=3,byrow=TRUE)
file_dates_md_matrix<-file_names_md_matrix[,2:3]

file_names_me <- readRDS("parallel_me_particle_files.rds")
file_names_me_list <- strsplit(file_names_me, split = "_")
file_names_me_matrix <- matrix(unlist(file_names_me_list),nrow=length(file_names_me_list),ncol=3,byrow=TRUE)
file_dates_me_matrix<-file_names_me_matrix[,2:3]

# now we expect that the created matrices using the names from mb,mc,md,me are the same.
# because in the end for a datetime and particle release height pair there should be 4 different particle traj documents. 
# but to make sure check if they are equal
check_1 <-setequal(file_dates_mb_matrix,file_dates_mc_matrix)
check_2<-setequal(file_dates_mb_matrix,file_dates_md_matrix)
check_3<-setequal(file_dates_mb_matrix,file_dates_me_matrix)


if (check_1==TRUE && check_2==TRUE && check_3==TRUE){ # if the name matrices are the same i.e.
  # if there are equal number of trajectory files for each instrument
  # we can easily do matrix operations and then take average

unique_times<-unique(file_names_mb_matrix[,2]) # chr [1:39]
# e.g.
# 202109050630 202109050645 202109050700 202109050715 .... 202109051545 202109051600

for (i in 1:1) { # this is to determine the particle release heights
  mytime<- unique_times[i] # pick a datetime 
  # mytime is "202109050630"
  idx_mytime<-which(file_names_mb_matrix[,2]==mytime)
  # now take a look at the huge matrix file_names_mb_matrix above
  # we want to determine the idx corresponding to this specific datetime
  # e.g. idx_mytime = 1 2 3 4 5 .... 13
  
  meters_of_idx_mytime <- strsplit( file_names_mb_matrix[idx_mytime,3], split = "m")
  # now determine the particle release heights of these idx
  # e.g. list of 13 
  # chr "1060" chr "1244" chr "1432" ....
  meters_of_idx_mytime <- matrix(unlist(meters_of_idx_mytime),nrow=13,ncol=1,byrow=TRUE)
  # now chr [1:13] "1060" "1244" "1432" ...
  meters_of_idx_mytime_num<-  as.numeric(meters_of_idx_mytime)
  # now num [1:13] 1060 1244 1432 1623 1818 186 2017 20 2220 355 527 701 879
  # now we want to get the vertical pressure scaling factors in the same order
  # remember vert_pres is in the increasing order i.e. 
  # first element of vert_pres corresponds to the scaling factor for the first atm. layer (20m) etc.
  sorted_vert_pres<-sort(meters_of_idx_mytime_num,decreasing = FALSE,index.return=TRUE)
  #$x
  #[1]   20  186  355  527  701  879 1060 1244 1432 1623 1818 2017 2220
}

# vert_pres[1] corresponds to 1st atm. layer, vert_pres[2] corresponds to 2nd atm. layer, vert_pres[3] corresponds to 3rd atm. layer etc.
vert_pres_sorted <-vert_pres[match(meters_of_idx_mytime_num,sorted_vert_pres$x)]
# now vert_pres_sorted[1] corresponds to 7th atm. layer i.e. 1060m, 
# vert_pres_sorted[2] corresponds to 8th atm. layer i.e. 1244 etc.

# normalize
vert_pres_sorted <- vert_pres_sorted/ sum(vert_pres_sorted)

# initialize
# these will store the vertical averaged [using the pressure scaling factors] distribution for one timepoint
mb_avg_alts_1 <- c()
mc_avg_alts_1 <- c()
md_avg_alts_1 <- c()
me_avg_alts_1 <- c()

mb_avg_alts_2 <- c()
mc_avg_alts_2 <- c()
md_avg_alts_2 <- c()
me_avg_alts_2 <- c()

mb_avg_alts_3 <- c()
mc_avg_alts_3 <- c()
md_avg_alts_3 <- c()
me_avg_alts_3 <- c()


# for each time point
for (i in 1:(length(file_names_mb_list)/13)) { # length(file_names_mb_list)/13=39, this is the number of timepoints for one day
    
    my_idx <- ((i-1)*13)+seq(1,13) # in each iteration pick the idx corresponding to one time point
    # e.g. mb_hist_ppb_errors_1[,my_idx] is of size [1:20,1:13]
    # 1:13 are the atm. layers and we know that now the particle release heights are matched ie. 
    # vert_pres_sorted[1] and mb_hist_ppb_errors_1[,1] correspond to the same atm. layer etc.
    mb_avg_of_alts_1 <-rowSums(mb_hist_ppb_errors_1[,my_idx]*vert_pres_sorted)
    # mb_avg_of_alts_1 num[1:20] : averaged using the pressure scaling factors, represents the histogram counts i.e. distribution
    # for a time point
    # then concatanate them, in the end we'll store each distribution for each time point
    mb_avg_alts_1<-cbind(mb_avg_alts_1,mb_avg_of_alts_1)
    
    # do the same for others
    # again there are 4 instruments mb mc md me and 3 available inventories for the same area [represented by the subscripts]
    # so this operation needs to be done 12 times. 
    mc_avg_of_alts_1 <-rowSums(mc_hist_ppb_errors_1[,my_idx]*vert_pres_sorted)
    mc_avg_alts_1<-cbind(mc_avg_alts_1,mc_avg_of_alts_1)
    
    md_avg_of_alts_1 <-rowSums(md_hist_ppb_errors_1[,my_idx]*vert_pres_sorted)
    md_avg_alts_1<-cbind(md_avg_alts_1,md_avg_of_alts_1)
    
    me_avg_of_alts_1 <-rowSums(me_hist_ppb_errors_1[,my_idx]*vert_pres_sorted)
    me_avg_alts_1<-cbind(me_avg_alts_1,me_avg_of_alts_1)
    
    
    
    mb_avg_of_alts_2 <-rowSums(mb_hist_ppb_errors_2[,my_idx]*vert_pres_sorted)
    mb_avg_alts_2<-cbind(mb_avg_alts_2,mb_avg_of_alts_2)
    
    mc_avg_of_alts_2 <-rowSums(mc_hist_ppb_errors_2[,my_idx]*vert_pres_sorted)
    mc_avg_alts_2<-cbind(mc_avg_alts_2,mc_avg_of_alts_2)
    
    md_avg_of_alts_2 <-rowSums(md_hist_ppb_errors_2[,my_idx]*vert_pres_sorted)
    md_avg_alts_2<-cbind(md_avg_alts_2,md_avg_of_alts_2)
    
    me_avg_of_alts_2 <-rowSums(me_hist_ppb_errors_2[,my_idx]*vert_pres_sorted)
    me_avg_alts_2<-cbind(me_avg_alts_2,me_avg_of_alts_2)
    
    
    
    
    mb_avg_of_alts_3 <-rowSums(mb_hist_ppb_errors_3[,my_idx]*vert_pres_sorted)
    mb_avg_alts_3<-cbind(mb_avg_alts_3,mb_avg_of_alts_3)
    
    mc_avg_of_alts_3 <-rowSums(mc_hist_ppb_errors_3[,my_idx]*vert_pres_sorted)
    mc_avg_alts_3<-cbind(mc_avg_alts_3,mc_avg_of_alts_3)
    
    md_avg_of_alts_3 <-rowSums(md_hist_ppb_errors_3[,my_idx]*vert_pres_sorted)
    md_avg_alts_3<-cbind(md_avg_alts_3,md_avg_of_alts_3)
    
    me_avg_of_alts_3 <-rowSums(me_hist_ppb_errors_3[,my_idx]*vert_pres_sorted)
    me_avg_alts_3<-cbind(me_avg_alts_3,me_avg_of_alts_3)
  
}

# mb_avg_alts_1 num [1:20,1:39] , 39 are the number of time points
# we need to find the mean of all 39 time points' distributions
# mb_daily_avg_1 num [1:20], now this is the ppb error distribution of all particle trajectory files' error simulations related to instrument mb
mb_daily_avg_1<-rowMeans(mb_avg_alts_1)

# do the same for others
# again there are 4 instruments mb mc md me and 3 available inventories for the same area [represented by the subscripts]
# so this operation needs to be done 12 times. 
mc_daily_avg_1<-rowMeans(mc_avg_alts_1)
md_daily_avg_1<-rowMeans(md_avg_alts_1)
me_daily_avg_1<-rowMeans(me_avg_alts_1)

mb_daily_avg_2<-rowMeans(mb_avg_alts_2)
mc_daily_avg_2<-rowMeans(mc_avg_alts_2)
md_daily_avg_2<-rowMeans(md_avg_alts_2)
me_daily_avg_2<-rowMeans(me_avg_alts_2)

mb_daily_avg_3<-rowMeans(mb_avg_alts_3)
mc_daily_avg_3<-rowMeans(mc_avg_alts_3)
md_daily_avg_3<-rowMeans(md_avg_alts_3)
me_daily_avg_3<-rowMeans(me_avg_alts_3)
  
# now take the average of the 4 distributions related to one day - take average of all instruments' error simulations
# daily_avg_1 num [1:20] ,  ppb error distribution of all particle trajectory files' error simulations when the inventory number 1 is used
# please remember that the subscript represents the inventory used
daily_avg_1 <- (mb_daily_avg_1 + mc_daily_avg_1 + md_daily_avg_1 + me_daily_avg_1)/4
daily_avg_2 <- (mb_daily_avg_2 + mc_daily_avg_2 + md_daily_avg_2 + me_daily_avg_2)/4
daily_avg_3 <- (mb_daily_avg_3 + mc_daily_avg_3 + md_daily_avg_3 + me_daily_avg_3)/4

# save
saveRDS(daily_avg_1,paste(Date,"_daily_city_and_river_smooth.rds",sep=""))
saveRDS(daily_avg_2,paste(Date,"_daily_corrected_elev_bothInstruments_city_and_river_smooth.rds",sep=""))
saveRDS(daily_avg_3,paste(Date,"_daily_corrected_raw_bothInstruments_city_and_river_smooth.rds",sep=""))
# The results for the chosen days can be found on:
# "/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Daily histograms"
}

