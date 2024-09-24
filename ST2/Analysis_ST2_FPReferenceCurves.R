
########### Analysis of ST2 (PFS & OS) Reference Curves using FP models ################

# Created by: Hugo Pedder
# Date: 2024-01-29

# Load packages
#devtools::install_github("hugaped/BristolTAG")
library(BristolTAG)
library(survival)
library(ggplot2)
library(tidyverse)
library(readxl)
library(R2jags)
source("NSCLC_functions.R")


#################################################################################
################### REFERENCE FRACTIONAL POLYNOMIAL MODEL - PFS ####################
#################################################################################


###### LUMELung1 - Non-squamous #######

ipd <- read.csv("ST2/ST2 IPD extraction.csv")

# Subset studies/treatments not to be included
ipd <- ipd %>%
  subset(study %in% c("LUME-Lung 1")) %>%
  subset(outcome %in% c("PFS")) %>%
  subset(histology %in% c("Non-Squamous"))

anova <- anova_data(timepoints=seq(0,16, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Docetaxel-Nintedanib"))



# Run 1st order models #

modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST2/JAGS_Results/ST2_PFS_fpoly1_LUMELung1.NS.RDS")

# Save results
saveRDS(modseq1, file="ST2/JAGS_Results/ST2_PFS_fpoly1_LUMELung1.NS.RDS")


# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST2/JAGS_Results/ST2_PFS_fpoly2_LUMELung1.NS.RDS")

# Save results
saveRDS(modseq2, file="ST2/JAGS_Results/ST2_PFS_fpoly2_LUMELung1.NS.RDS")





###### LUMELung1 - Squamous #######

ipd <- read.csv("ST2/ST2 IPD extraction.csv")

# Subset studies/treatments not to be included
ipd <- ipd %>%
  subset(study %in% c("LUME-Lung 1")) %>%
  subset(outcome %in% c("PFS")) %>%
  subset(histology %in% c("Squamous"))

anova <- anova_data(timepoints=seq(0,16, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Docetaxel-Nintedanib"))



# Run 1st order models #

modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST2/JAGS_Results/ST2_PFS_fpoly1_LUMELung1.S.RDS")

# Save results
saveRDS(modseq1, file="ST2/JAGS_Results/ST2_PFS_fpoly1_LUMELung1.S.RDS")


# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST2/JAGS_Results/ST2_PFS_fpoly2_LUMELung1.S.RDS")

# Save results
saveRDS(modseq2, file="ST2/JAGS_Results/ST2_PFS_fpoly2_LUMELung1.S.RDS")





