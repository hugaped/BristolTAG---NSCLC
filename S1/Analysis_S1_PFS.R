
########### Analysis of S1 PFS using Bristol TAG ################

# Created by: Hugo Pedder
# Date: 2023-12-07

# Script runs and saves fractional polynomial models. Results can be explored here,
# but use of .Rmd template to generate .docx outputs is easier for the report.

# Load packages
#devtools::install_github("hugaped/BristolTAG")
library(BristolTAG)
library(survival)
library(ggplot2)
library(tidyverse)


#### Load and prepare data ####

ipd <- read.csv("S1/S1 PFS IPD extraction.csv")


anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("PPCT", "PCT"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)





###### Test run for convergence #########
jagsdat$P1 <- -0.5

# Model does not include (or monitor) rankings - can easily be done in R
jagsmod <- jags(data=jagsdat, inits=inits,
                parameters.to.save=c("d", "mu", "dev", "totresdev"),
                model.file=system.file("JAGSmodels", "FE_1st_order_model.jags", package="BristolTAG"),
                jags.seed=890421,
                n.iter=80000, n.burnin=30000, n.thin=10
)




##### Run 1st order models #####

modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="S1/JAGS_Results/S1_PFS_fpoly1.RDS")

# Save results
saveRDS(modseq1, file="S1/JAGS_Results/S1_PFS_fpoly1.RDS")



##### Run 2nd order models #####

modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="S1/JAGS_Results/S1_PFS_fpoly2.RDS")

# Save results
saveRDS(modseq2, file="S1/JAGS_Results/S1_PFS_fpoly2.RDS")




