
########### Analysis of NS1 PFS using Bristol TAG ################

# Created by: Hugo Pedder
# Date: 2023-10-09

# Script runs and saves fractional polynomial models. Results can be explored here,
# but use of .Rmd template to generate .docx outputs is easier for the report.

# Load packages
#devtools::install_github("hugaped/BristolTAG")
library(BristolTAG)
library(survival)
library(ggplot2)
library(tidyverse)



################## ANALYSIS #######################

#### Load and prepare data ####

ipd <- read.csv("NS1/NS1 PFS IPD extraction.csv")
ipd <- subset(ipd, !study %in% c("Egyptian1", "Halmos MAIC"))

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("PPCT", "PCT", "ABCP", "BCP"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)





###### Test run for convergence #########

# jagsdat <- fp_data(anova, trtnames=c("PCT", "PPCT", "ABCP", "BCP"))
# inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=2, seed=890421)
#
# jagsdat$P1 <- -0.5
# jagsdat$P2 <- -0.5
#
# # Model does not include (or monitor) rankings - can easily be done in R
# jagsmod <- jags(data=jagsdat, inits=inits,
#                 parameters.to.save=c("d", "mu", "dev", "totresdev"),
#                 model.file=system.file("JAGSmodels", "FE_2nd_order_model.jags", package="BristolTAG"),
#                 jags.seed=890421,
#                 n.iter=150000, n.burnin=50000, n.thin=10
# )
# # Set attributes to match those of jagsdat
# attr(jagsmod, "trtnames") <- attr(jagsdat, "trtnames")
# attr(jagsmod, "studynames") <- attr(jagsdat, "studynames")
# attr(jagsmod, "ipd") <- attr(jagsdat, "ipd")
#
#
#
# # Change netref
# jagsdat <- fp_data(anova, trtnames=c("BCP", "PPCT", "PCT", "ABCP"))
# inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=2, seed=890421)
#
# jagsdat$P1 <- -0.5
# jagsdat$P2 <- -0.5
#
# # Model does not include (or monitor) rankings - can easily be done in R
# jagsmod2 <- jags(data=jagsdat, inits=inits,
#                 parameters.to.save=c("d", "mu", "dev", "totresdev"),
#                 model.file=system.file("JAGSmodels", "FE_2nd_order_model.jags", package="BristolTAG"),
#                 jags.seed=890421,
#                 n.iter=150000, n.burnin=50000, n.thin=10
# )
# # Set attributes to match those of jagsdat
# attr(jagsmod2, "trtnames") <- attr(jagsdat, "trtnames")
# attr(jagsmod2, "studynames") <- attr(jagsdat, "studynames")
# attr(jagsmod2, "ipd") <- attr(jagsdat, "ipd")
#
#
# # Change netref
# jagsdat <- fp_data(anova, trtnames=c("PPCT", "PCT", "ABCP", "BCP"))
# inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=2, seed=890421)
#
# jagsdat$P1 <- -0.5
# jagsdat$P2 <- -0.5
#
# # Model does not include (or monitor) rankings - can easily be done in R
# jagsmod3 <- jags(data=jagsdat, inits=inits,
#                  parameters.to.save=c("d", "mu", "dev", "totresdev"),
#                  model.file=system.file("JAGSmodels", "FE_2nd_order_model.jags", package="BristolTAG"),
#                  jags.seed=890421,
#                  n.iter=150000, n.burnin=50000, n.thin=10
# )
# # Set attributes to match those of jagsdat
# attr(jagsmod3, "trtnames") <- attr(jagsdat, "trtnames")
# attr(jagsmod3, "studynames") <- attr(jagsdat, "studynames")
# attr(jagsmod3, "ipd") <- attr(jagsdat, "ipd")
#
# surv1 <- survcalc(jagsmod, refstudy="KEYNOTE189")
# surv2 <- survcalc(jagsmod2, refstudy="KEYNOTE189")
# surv3 <- survcalc(jagsmod3, refstudy="KEYNOTE189")
#
# hr1 <- hrcalc(jagsmod)
# hr2 <- hrcalc(jagsmod2)
# hr3 <- hrcalc(jagsmod3)
#
# # AUC are same
# auc1 <- auc(surv1)
# auc2 <- auc(surv2)
# auc3 <- auc(surv3)
#
# # Implied survival are same
# plot(surv1)
# plot(surv2)
# plot(surv3)
#
# # Implied HRs are same
# plot(hr1, reftrt="PCT")
# plot(hr2, reftrt="PCT")
# plot(hr3, reftrt="PCT")
#







##### Run 1st order models #####

path <- "NS1/JAGS_Results/NS1_PFS_fpoly1_EgyptExcl.RDS"
modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile=path,
                          overwrite=TRUE)

# Save results
saveRDS(modseq1, file=path)



##### Run 2nd order models #####

path <- "NS1/JAGS_Results/NS1_PFS_fpoly2_EgyptExcl.RDS"
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile=path)

# Save results
saveRDS(modseq2, file=path)











################## ANALYSIS USING MAIC #######################

#### Load and prepare data ####

ipd <- read.csv("NS1/NS1 PFS IPD extraction.csv")
ipd <- subset(ipd, study %in% c("KEYNOTE189", "Halmos MAIC"))

anova <- anova_data(timepoints=seq(0,72, by=1),
                    df=ipd)

jagsdat <- fp_data(anova, trtnames=c("PPCT", "PCT", "ABCP"))

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=1, seed=890421)





###### Test run for convergence #########
# jagsdat$P1 <- -0.5
#
# # Model does not include (or monitor) rankings - can easily be done in R
# jagsmod <- jags(data=jagsdat, inits=inits,
#                 parameters.to.save=c("d", "mu", "dev", "totresdev"),
#                 model.file=system.file("JAGSmodels", "FE_1st_order_model.jags", package="BristolTAG"),
#                 jags.seed=890421,
#                 n.iter=150000, n.burnin=50000, n.thin=10
# )
#
#


##### Run 1st order models #####

path <- "NS1/JAGS_Results/NS1_PFS_fpoly1_MAIC.RDS"
modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1, inits=inits,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile=path)

# Save results
saveRDS(modseq1, file=path)



##### Run 2nd order models #####

path <- "NS1/JAGS_Results/NS1_PFS_fpoly2_MAIC.RDS"
modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile=path)

# Save results
saveRDS(modseq2, file=path)




