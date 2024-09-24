
########### Analysis of ST1 PFS for all PDL1 thresholds ################

# Created by: Hugo Pedder
# Date: 2023-12-20

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
hr.df <- read_excel("ST1_HazardRatios.xlsx")

hr.df <- hr.df %>% mutate(histology=case_when(histology=="other" ~ "non-squamous",
                                    TRUE ~ histology))

trtnames <- c("Docetaxel",
              "Pemetrexed",
              "Pemetrexed + Carboplatin",
              "Selpercatinib")

hr.df <- hr.df %>% mutate(lhr = log(hr),
                          se.lhr = (log(u95_hr) - log(l95_hr)) / (1.96 * 2),
                          histology = factor(histology, labels=c("NS", "S", "Mixed"),
                                             levels=c("non-squamous", "squamous", "mixed")),
                          treatment1 = match(treatment1, trtnames),
                          treatment2 = match(treatment2, trtnames)
                          )



#############################################################################
################# PROPORTIONAL HAZARDS IGNORING HISTOLOGY ##################
#############################################################################

# Using Selpercatinib MAIC (Genetic matching algorithm)

hr.mix <- subset(hr.df, histology=="NS" |
                   (studyID=="LIBRETTO-001" & Notes=="Genetic-matching")) %>%
  mutate(histology="Mixed")

jagsdat <- getjags_NS2(hr.mix, type="overall")

jagsmod <- jags(data=jagsdat, #inits=inits,
                parameters.to.save=c("d", "resdev", "totresdev"),
                model.file="JAGSmodels/NMA_HR.jags",
                jags.seed=890421,
                n.iter=15000, n.burnin=5000, n.thin=1
)
attr(jagsmod, "trtnames") <- trtnames

d <- jagsmod$BUGSoutput$sims.list$d

plot.df <- data.frame(x=as.vector(d[,]),
                      treatment=factor(rep(rep(trtnames, each=jagsmod$BUGSoutput$n.sims),2))
)

ggplot(plot.df, aes(x=x, y=treatment)) +
  ggdist::stat_halfeye() +
  theme_bw()

saveRDS(jagsmod, file="ST1/JAGS_Results/ST1_PFS_HR_Non-squamous_Genetic.RDS")





# Using Selpercatinib MAIC (TMLE algorithm)

hr.mix <- subset(hr.df, histology=="NS" |
                   (studyID=="LIBRETTO-001" & Notes=="TMLE")) %>%
  mutate(histology="Mixed")

jagsdat <- getjags_NS2(hr.mix, type="overall")

jagsmod <- jags(data=jagsdat, #inits=inits,
                parameters.to.save=c("d", "resdev", "totresdev"),
                model.file="JAGSmodels/NMA_HR.jags",
                jags.seed=890421,
                n.iter=15000, n.burnin=5000, n.thin=1
)
attr(jagsmod, "trtnames") <- trtnames

d <- jagsmod$BUGSoutput$sims.list$d

plot.df <- data.frame(x=as.vector(d[,]),
                      treatment=factor(rep(rep(trtnames, each=jagsmod$BUGSoutput$n.sims),2))
)

ggplot(plot.df, aes(x=x, y=treatment)) +
  ggdist::stat_halfeye() +
  theme_bw()

saveRDS(jagsmod, file="ST1/JAGS_Results/ST1_PFS_HR_Non-squamous_TMLE.RDS")



#################################################################################
################### FRACTIONAL POLYNOMIAL IGNORING HISTOLOGY ####################
#################################################################################


######## Fractional Polynomials ignoring histology (INCLUDING CM078, SPLITTING PDL1 >1% effect) ###########

ipd <- read.csv("ST1/ST1 PFS IPD extraction.csv")

# Subset studies/treatments not to be included
ipd <- ipd %>% subset(!study %in% c("NVALT7", # To be included as histology-specific HR
                                    "GOIRC", # To be included as histology-specific HR
                                    "LIBRETTO-001" # To be included as constant HR
)
) %>%

  # Subset LUMElung1 pooled histology
  subset(!(study %in% "LUME-Lung 1" & histology!="Pooled")) %>%

  # Subset PDL1 for PDL1-targetting studies
  subset((PDL1 %in% c("NR", "Mixed", ">50%", ">10%", "1-49%", "1-9%", "<1%")) |
           (study %in% "CHECKMATE078" & PDL1 %in% ">1%")) %>%

  # Ensure CM078 >1% is at the bottom of the dataset
  mutate(study=case_when(study %in% c("CHECKMATE078") #& PDL1 %in% c(">1%")
                         ~ paste0("ZZ", study, "_>1%"),
                         TRUE ~ study)) %>%

  # Add pdlcat variable
  mutate(pdl1cat=case_when(PDL1 %in% c(">50%", ">10%", ">1%", "Mixed", "NR") ~ 1,
                           PDL1 %in% c("1-49%", "1-9%") ~ 2,
                           PDL1 %in% c("<1%") ~ 3))

# Repeat rows of data for non-PDL1 specific studies
reps <- subset(ipd, PDL1 %in% c("Mixed", "NR"))

# Bind duplicate rows to ensure they are included in all PDL1 analyses
ipd <- rbind(ipd, reps %>% mutate(pdl1cat=2))
ipd <- rbind(ipd, reps %>% mutate(pdl1cat=3))


# FOR NOW, DROP PDL1 <1%
#ipd <- subset(ipd, pdl1cat!=3)


anova <- anova_data_st1(timepoints=seq(0,71, by=1),
                    df=ipd)

# Ensures additional Checkmate rows are ordered last
anova$aggregate <- arrange(anova$aggregate, trialid, pdl1cat)

jagsdat <- fp_data_st1(anova, trtnames=c("Docetaxel",
                                     "Atezolizumab",
                                     "Pembrolizumab",
                                     "Nivolumab",
                                     "Docetaxel-Nintedanib",
                                     "Sotorasib",
                                     "Doctxl-PCT",
                                     "Pacltxl-Carb"))

### Additional bits for split PDL1 ###
# Add proportion of PDL1 thresholds in mixed study (CM078)
jagsdat$p10pct <- c(rep(NA, jagsdat$ns.pdl),
                    0.6349) # Hash out if excluding CM078



# Test run
jagsdat$P1 <- -0.5
inits <- fp_geninits_st1(ns=jagsdat$ns, nt=jagsdat$nt,
                         npdl=3, # Number of PDL1 categories to analyse
                         polyorder=1, seed=890421)
jagsmod <- jags(data=jagsdat, inits=inits,
                parameters.to.save=c("d", "mu", "dev", "totresdev"),
                model.file="JAGSmodels/FE_1st_order_AllPDL1split.jags",
                #model.file=system.file("JAGSmodels", "FE_1st_order_model.jags", package="BristolTAG"),
                jags.seed=890421,
                n.iter=150000, n.burnin=50000, n.thin=10
)




# Run 1st order models #

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits_st1(ns=jagsdat$ns, nt=jagsdat$nt,
                         npdl=3, # Number of PDL1 categories to analyse
                         polyorder=1, seed=890421)


modseq1 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=1,
                          inits=inits,
                          jagsfile="JAGSmodels/FE_1st_order_AllPDL1split.jags",
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_AllPDL_PFS_fpoly1.RDS")

# Save results
saveRDS(modseq1, file="ST1/JAGS_Results/ST1_AllPDL_PFS_fpoly1.RDS")




# Run 2nd order models #

# Initial values (only rerun if you want to change/set these)
inits <- fp_geninits_st1(ns=jagsdat$ns, nt=jagsdat$nt, polyorder=2,
                     npdl=3, # Number of PDL1 categories to analyse
                     seed=890421)

modseq2 <- sequence_fpoly(jagsdat, powers=c(-3,-2,-1,-0.5, 0, 0.5, 1, 2, 3), polyorder=2,
                          n.iter=150000, n.burnin=50000, n.thin=10,
                          inits=inits,
                          jagsfile="JAGSmodels/FE_2nd_order_AllPDL1split.jags",
                          jags.seed=890421,
                          savefile="ST1/JAGS_Results/ST1_AllPDL_PFS_fpoly2.RDS")

# Save results
saveRDS(modseq2, file="ST1/JAGS_Results/ST1_AllPDL_PFS_fpoly2.RDS")




