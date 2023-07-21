# ---------------------------------------------------------
# TUM - Technichal University of Munich
#
# Authors:  Aydin Uzun
# Date: 2022
# Purpose:  Read the scaling factors from the receptors' .rds documents
# ---------------------------------------------------------
rm(list = ls())
# set directory
setwd("/Volumes/esm/data/Footprint_STILT/Hamburg/ERA5/V2_small") 
directories<-list.dirs(path = "/Volumes/esm/data/Footprint_STILT/Hamburg/ERA5/V2_small", full.names = FALSE, recursive = FALSE)
full_directories<-list.dirs(path = "/Volumes/esm/data/Footprint_STILT/Hamburg/ERA5/V2_small", full.names = TRUE, recursive = FALSE)

# initialize
scaling_matrix <- c()

for (i in 1:length(full_directories)) {
  
  # get into the directories
  setwd(full_directories[i])
  # name of the RDS document
  name_of_RDS_docu <- paste(directories[i],"receptors","ERA5.rds", sep = "_", collapse = NULL)
  receptors_file_RDS <- readRDS(name_of_RDS_docu)
  scaling <- receptors_file_RDS$scaling_factors
  # scaling factor for the corresponding day
  scaling_fac = unique(scaling)
  setwd("/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Scaling_factors")
  # save the scaling factor for the corresponding day
  write.csv(scaling_fac, file = paste(directories[i],"_scaling.csv",sep=""))
  scaling_matrix <- cbind(scaling_matrix,scaling_fac)
}

setwd("/Volumes/esm/11-Thesis/03-Scientific-Internship/2021 FP Aydin Uzun/Data/Scaling_factors")
scaling_matrix<- t(scaling_matrix)
# save the whole scaling factors
write.csv(scaling_matrix, file = "Scaling_factors.csv")


