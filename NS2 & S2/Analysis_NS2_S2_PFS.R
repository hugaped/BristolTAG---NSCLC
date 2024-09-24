
########### Analysis of NS2 & S2 PFS ################

# Created by: Hugo Pedder
# Date: 2023-10-25

# Load packages
#devtools::install_github("hugaped/BristolTAG")
library(BristolTAG)
library(survival)
library(ggplot2)
library(tidyverse)
library(readxl)
library(R2jags)
source("NSCLC_functions.R")


# Load HR data
hr.df <- read_excel("NS2_S2_HazardRatios.xlsx")

trtnames <- c("Platinum Chemotherapy",
              "Atezolizumab",
              "Pembrolizumab",
              "Pembrolizumab + Chemotherapy")

hr.df <- hr.df %>% mutate(lhr = log(hr),
                          se.lhr = (log(u95_hr) - log(l95_hr)) / (1.96 * 2),
                          histology = factor(histology, labels=c("NS", "S", "Mixed"),
                                             levels=c("non-squamous", "squamous", "mixed")),
                          treatment1 = match(treatment1, trtnames),
                          treatment2 = match(treatment2, trtnames)
                          )


#################################################################################
################### REFERENCE FRACTIONAL POLYNOMIAL MODEL ####################
#################################################################################

###### KEYNOTE 189 #######

ipd <- read.csv("NS2 & S2/NS2 S2 PFS IPD extraction.csv")
ipd <- subset(ipd, study=="KEYNOTE189")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("PPCT",
                                     "Chemotherapy"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_KEYNOTE189.RDS")

# Save results
saveRDS(modseq1, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_KEYNOTE189.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_KEYNOTE189.RDS")

# Save results
saveRDS(modseq2, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_KEYNOTE189.RDS")




###### KEYNOTE 407 #######

ipd <- read.csv("NS2 & S2/NS2 S2 PFS IPD extraction.csv")
ipd <- subset(ipd, study=="KEYNOTE407")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("PPCT",
                                     "Chemotherapy"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_KEYNOTE407.RDS")

# Save results
saveRDS(modseq1, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_KEYNOTE407.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=80000, n.burnin=30000, n.thin=10,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_KEYNOTE407.RDS")

# Save results
saveRDS(modseq2, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_KEYNOTE407.RDS")




#################################################################################
################### FRACTIONAL POLYNOMIAL IGNORING HISTOLOGY ####################
#################################################################################


######## Analysis 1: Fractional Polynomials ignoring histology (incl KEYNOTE-042) ###########

ipd <- read.csv("NS2 & S2/NS2 S2 PFS IPD extraction.csv")

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)
jagsdat <- fp_data(anova, trtnames=c("PPCT",
                                     "Atezolizumab", "Pembrolizumab", "Chemotherapy"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)


# Run 1st order models #

modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=150000, n.burnin=50000, n.thin=25,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_Analysis1.2.RDS")

# Save results
saveRDS(modseq1, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_Analysis1.2.RDS")

# Run 2nd order models #
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=25,
                          jags.seed=890421,
                          savefile="Analysis/NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_Analysis1.2.RDS")

# Save results
saveRDS(modseq2, file="Analysis/NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_Analysis1.2.RDS")




#################################################################################
################### FRACTIONAL POLYNOMIAL WEIGHTED HISTOLOGY ####################
#################################################################################

hist.LHR <- readRDS(file="NS2 & S2/JAGS_Results/NS2S2_PFS_HR_hist.LHR.RDS")
quantile(hist.LHR, probs=c(0.025, 0.5, 0.975))
exp(quantile(hist.LHR, probs=c(0.025, 0.5, 0.975)))


######## Analysis 2: Fractional Polynomial weighted histology (incl KEYNOTE-042) ###########

hist.LHR <- readRDS(file="NS2 & S2/JAGS_Results/NS2S2_PFS_HR_hist.LHR.RDS")
quantile(hist.LHR, probs=c(0.025, 0.5, 0.975))
exp(quantile(hist.LHR, probs=c(0.025, 0.5, 0.975)))

ipd <- read.csv("NS2 & S2/NS2 S2 PFS IPD extraction.csv")
#ipd <- subset(ipd, study!="KEYNOTE042")


anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("PPCT",
                                     "Atezolizumab", "Pembrolizumab", "Chemotherapy"))

# Proportion with squamous histology
jagsdat$phist <- c("IMPOWER110"=0.242,
                   "KEYNOTE024"=0.184,
                   "KEYNOTE042"=0.385,
                   "KEYNOTE189"=0,
                   "KEYNOTE407"=1)

# Set hist.LHR parameters
jagsdat$m.hist.LHR <- mean(hist.LHR)
jagsdat$prec.hist.LHR <- 1/var(hist.LHR)


# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)

# Add initial values for d
for (i in seq_along(inits)) {
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d)))
}


###### Run 1st order models in sequence ######
modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          jagsfile="JAGSmodels/FE_1st_order_MixedHistology.jags",
                          parameters.to.save=c("d", "mu", "hist.LHR", "dev", "totresdev"),
                          n.iter=150000, n.burnin=50000, n.thin=25,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_Analysis2.RDS")

# Save results
saveRDS(modseq1, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_Analysis2.RDS")



###### Run 2nd order models ######
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=2, seed=890421)

# Add initial values for d
for (i in seq_along(inits)) {
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d))) # for 1st coefficient
}


modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2, inits=inits,
                          jagsfile="JAGSmodels/FE_2nd_order_MixedHistology.jags",
                          parameters.to.save=c("d", "mu", "hist.LHR", "dev", "totresdev"),
                          n.iter=150000, n.burnin=50000, n.thin=25,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_Analysis2.RDS")

# Save results
saveRDS(modseq2, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_Analysis2.RDS")




######## Analysis 3: Fractional Polynomial fully weighted histology (incl KEYNOTE-042) ###########

ipd <- read.csv("NS2 & S2/NS2 S2 PFS IPD extraction.csv")
#ipd <- subset(ipd, study!="KEYNOTE042")


anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("PPCT",
                                     "Atezolizumab", "Pembrolizumab", "Chemotherapy"))

# Proportion with squamous histology
jagsdat$phist <- c("IMPOWER110"=0.242,
                   "KEYNOTE024"=0.184,
                   "KEYNOTE042"=0.385,
                   "KEYNOTE189"=0,
                   "KEYNOTE407"=1)

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)

# Add initial values for d (two sets for 2x FP params)
for (i in seq_along(inits)) {
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d)))
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d)))
}


###### Run 1st order models in sequence ######
modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          jagsfile="JAGSmodels/FE_1st_order_MixedHistology_AllFP.jags",
                          parameters.to.save=c("d", "mu", "hist.LHR", "dev", "totresdev"),
                          n.iter=150000, n.burnin=50000, n.thin=25,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_Analysis3.RDS")

# Save results
saveRDS(modseq1, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly1_Analysis3.RDS")



###### Run 2nd order models ######
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=2, seed=890421)

# Add initial values for d
for (i in seq_along(inits)) {
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d))) # for 1st coefficient
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d))) # for 2nd coefficient
  inits[[i]]$d <- cbind(inits[[i]]$d, rep(NA, nrow(inits[[i]]$d))) # for 3rd coefficient
}


modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2, inits=inits,
                          jagsfile="JAGSmodels/FE_2nd_order_MixedHistology_AllFP.jags",
                          parameters.to.save=c("d", "mu", "hist.LHR", "dev", "totresdev"),
                          n.iter=150000, n.burnin=50000, n.thin=25,
                          jags.seed=890421,
                          savefile="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_Analysis3.RDS")

# Save results
saveRDS(modseq2, file="NS2 & S2/JAGS_Results/NS2_S2_PFS_fpoly2_Analysis3.RDS")



