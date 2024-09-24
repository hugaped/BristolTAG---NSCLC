
########### Analysis of ST1 PFS Reference Curves using FP models ################

# Created by: Hugo Pedder
# Date: 2024-01-08

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
################### REFERENCE FRACTIONAL POLYNOMIAL MODEL ####################
#################################################################################


###### CHECKMATE-057 (Non-squamous, Overall (pooled PDL1 subgroups)) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CHECKMATE057" & PDL1=="Pooled")

anova <- anova_data(timepoints=seq(0,73, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Nivolumab"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CM057.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CM057.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=120000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CM057.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CM057.RDS")







###### CHECKMATE-017 (Aquamous, Overall (pooled PDL1 subgroups)) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CHECKMATE017" & PDL1=="Pooled")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Nivolumab"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CM017.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CM017.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CM017.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CM017.RDS")






###### LIBRETTO-001 (Mixed, Selpercatinib reference) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="LIBRETTO-001")

# Add "imaginary" replicate arm to allow estimation using standard FP models
#ipd <- rbind(ipd, ipd %>% mutate(treatment="dummy"))

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
# jagsdat <- fp_data(anova, trtnames=c("Selpercatinib",
#                                      "dummy"))
jagsdat <- fp_data(anova, trtnames=c("Selpercatinib"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jagsfile = "JAGSmodels/FE_1st_order_model_ref.jags",
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_LIBRETTO-001.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_LIBRETTO-001.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jagsfile = "JAGSmodels/FE_2nd_order_model_ref.jags",
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_LIBRETTO-001.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_LIBRETTO-001.RDS")






###### CODEBREAK200 (Non-squamous, Docetaxel) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CODEBREAK200")

anova <- anova_data(timepoints=seq(0,73, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Sotorasib"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CODEBREAK200.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CODEBREAK200.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=120000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CODEBREAK200.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CODEBREAK200.RDS")










############################################################
############ REFERENCE STUDIES BY PDL1 SUBGROUP ############
############################################################




###### CHECKMATE-057 (Non-squamous, PDL1 >50%) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CHECKMATE057" & PDL1==">10%")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Nivolumab"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CM057_+50pct.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CM057_+50pct.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CM057_+50pct.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CM057_+50pct.RDS")





###### CHECKMATE-057 (Non-squamous, PDL1 1-49%) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CHECKMATE057" & PDL1=="1-9%%")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Nivolumab"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CM057_1-49pct.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CM057_1-49pct.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CM057_1-49pct.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CM057_1-49pct.RDS")






###### CHECKMATE-017 (Aquamous, PDL1 >50%) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CHECKMATE017" & PDL1==">10%")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Nivolumab"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CM017_+50pct.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CM017_+50pct.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CM017_+50pct.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CM017_+50pct.RDS")





###### CHECKMATE-017 (Squamous, PDL1 1-49%) #######

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

ipd <- subset(ipd, study=="CHECKMATE017" & PDL1=="1-9%%")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("Docetaxel",
                                     "Nivolumab"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly1_CM017_1-49pct.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_PFS_fpoly1_CM017_1-49pct.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_PFS_fpoly2_CM017_1-49pct.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_PFS_fpoly2_CM017_1-49pct.RDS")
