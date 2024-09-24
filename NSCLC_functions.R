
########### Specific R functions for analysis of NSCLC ################

#' @param df A data frame containing the following variables:
#'  * studyID
#'  * treatment1
#'  * treatment2
#'  * lhr (log-HR for treatment 2 vs treatment 1)
#'  * se.lhr (SE for log-HR for treatment 2 vs treatment 1)
#'  * hist (histology: can take `"squamous"`, `"non-squamous"` or `"mixed"`)
#'  * p.squamous (the proportion of squamous in the study)
#'  @param type can take either `"hist"` for data that allows for histology effect,
#'  or `"overall"` to estimate histology lumped effects using NMA
#'
getjags_NS2 <- function(df, type) {

  df <- arrange(df, histology)

  # Assign studyID
  studynames <- unique(df$studyID)
  df <- df %>% mutate(studyID=match(studyID, studynames))

  # Get jags data
  ns <- dplyr::n_distinct(df$studyID)

  jagsdat <- list()

  if (type=="hist") {
    jagsdat$ns2 <- nrow(subset(df, histology=="Mixed"))
    jagsdat$ns1 <- ns - jagsdat$ns2

    jagsdat$lhr <- matrix(ncol=2, nrow=ns)
    jagsdat$se <- matrix(ncol=2, nrow=ns)

  } else if (type=="overall") {

    jagsdat$lhr <- vector()
    jagsdat$se <- vector()

    jagsdat$ns <- ns
  } else {
    stop("type must be `hist` or `overall`")
  }


  jagsdat$nt <- n_distinct(c(df$treatment1, df$treatment2))

  jagsdat$t <- matrix(ncol=2, nrow=ns)

  if (type=="hist") {
    jagsdat$pS <- vector()
  }

  for (i in 1:ns) {
    study <- subset(df, studyID==i)

    if (!"Mixed" %in% study$histology) {

      if (type=="hist") {
        if ("NS" %in% study$histology) {
          jagsdat$lhr[i,1] <- subset(study, histology=="NS")$lhr
          jagsdat$se[i,1] <- subset(study, histology=="NS")$se.lhr
        }

        if ("S" %in% study$histology) {
          jagsdat$lhr[i,2] <- subset(study, histology=="S")$lhr
          jagsdat$se[i,2] <- subset(study, histology=="S")$se.lhr
        }

        jagsdat$pS <- append(jagsdat$pS, NA)
      }

    } else {

      if (type=="hist") {
        jagsdat$lhr[i,1:2] <- subset(study, histology=="Mixed")$lhr
        jagsdat$se[i,1:2] <- subset(study, histology=="Mixed")$se.lhr

        jagsdat$pS <- append(jagsdat$pS, subset(study, histology=="Mixed")$p.squamous)

      } else if (type=="overall") {
        jagsdat$lhr[i] <- subset(study, histology=="Mixed")$lhr
        jagsdat$se[i] <- subset(study, histology=="Mixed")$se.lhr
      }
    }

    jagsdat$t[i,] <- c(study$treatment1[1], study$treatment2[1])
  }

  jagsdat$studyID <- studynames

  return(jagsdat)
}





#' @param hist Can take either "NS" or "S" to select fractional polynomial
#' parameters corresponding to non-squamous or squamous from shared parameter
#' model
#' @param jagsmod A fractional polynomial model
#'
#' @details
#' Shared parameter FP model must only share info on the 1st FP parameter (the
#' intercept - correspondingish to the HR)
select_histology <- function(jagsmod, hist="NS", analysis=3.1) {
  if (!class(jagsmod) %in% "rjags") {
    stop("jagsmod must be an object of class 'rjags'")
  }

  dlen <- dim(jagsmod$BUGSoutput$sims.list$d)[3]
  if (hist=="NS") {
    if (round(analysis)==3) {
      dropdim <- dim(jagsmod$BUGSoutput$sims.list$d)[3]
      jagsmod$BUGSoutput$sims.list$d <- jagsmod$BUGSoutput$sims.list$d[,,-dropdim]
    } else if (round(analysis)==5) {
      jagsmod$BUGSoutput$sims.list$d <- jagsmod$BUGSoutput$sims.list$d[,,1:(dlen/2)]
    }
  } else if (hist=="S") {
    if (round(analysis)==3) {
      dropdim <- dim(jagsmod$BUGSoutput$sims.list$d)[3]
      jagsmod$BUGSoutput$sims.list$d <-
        jagsmod$BUGSoutput$sims.list$d[,,c(dropdim, 2:(dropdim-1))]
    } else if (round(analysis)==5) {
      #keepdim2 <- dim(jagsmod$BUGSoutput$sims.list$d)[3]/2
      jagsmod$BUGSoutput$sims.list$d <-
        jagsmod$BUGSoutput$sims.list$d[,,((dlen/2)+1):dlen]
    }

  } else {
    stop("hist must take either 'NS' or 'S'")
  }
  return(jagsmod)
}




#' @param pdl1 Can take either `">50%"`, `"1-49%"` or `"<1%"` to select fractional
#' polynomial parameters corresponding to PDL1 categories from shared parameter
#' model
#' @param jagsmod A fractional polynomial model
#'
select_pdl1 <- function(jagsmod, pdl1=">50%") {
  if (!class(jagsmod) %in% "rjags") {
    stop("jagsmod must be an object of class 'rjags'")
  }

  ddim <- dim(jagsmod$BUGSoutput$sims.list$d)[2]

  if (pdl1==">50%") {
    pdl1cat <- 1
  } else if (pdl1=="1-49%") {
    pdl1cat <- 2
  } else if (pdl1=="<1%") {
    pdl1cat <- 3
  } else {
    stop("pdl1 must take either '>50%', '1-49%', or <1%")
  }

  if (pdl1cat>ddim) {
    stop("Dimensions of d do not match those specified in `pdl1`")
  }

  jagsmod$BUGSoutput$sims.list$mu <- jagsmod$BUGSoutput$sims.list$mu[,,pdl1cat,]
  jagsmod$BUGSoutput$sims.list$d <- jagsmod$BUGSoutput$sims.list$d[,,pdl1cat,]

  return(jagsmod)
}





get_relmat <- function(jagsmod, hist="NS", eform=TRUE, pretty=TRUE,
                       digits=2) {

  if (!is.null(hist)) {
    if (hist=="NS") {
      ind <- 1
    } else if (hist=="S") {
      ind <- 2
    }

    d <- jagsmod$BUGSoutput$sims.list$d[,,ind]

  } else {
    d <- jagsmod$BUGSoutput$sims.list$d[,]
  }

  if (pretty==TRUE) {
    mat <- matrix(nrow=dim(d)[2], ncol=dim(d)[2])
  } else {
    mat <- array(dim = c(dim(d)[2], dim(d)[2], jagsmod$BUGSoutput$n.sims))
  }

  for (i in 1:dim(d)[2]) {
    for (k in 1:dim(d)[2]) {
      diff <- d[,i] - d[,k]

      if (eform==TRUE) {
        diff <- exp(diff)
      }

      if (pretty==FALSE) {
        mat[i,k,] <- diff
      } else if (pretty==TRUE) {
        diff <- BristolTAG:::mcmc_sum(diff)
        diff <- round(diff, digits = digits)
        mat[i,k] <- paste0(diff[4], " (", diff[3], ", ", diff[5], ")")
      }
    }
  }
  return(mat)
}




#' Copy from BristolTAG::survcalc() but applies constant HRs to generate survival
#' predictions
#'
#' @param trtnames must be a string of treatment names matching those in d
#' @param d contains MCMC results from a constant HR model (note that if it
#' contains additional dimensions for histology then this should be a subset such
#' that histology has already been selected for)
#' @param jagsmod contains MCMC results from a fractional polynomial analysis
#' of a single reference study
#' @param refstudy is the name of the refernce study on which `jagsmod` was
#' estimated
#'
hr_survival <- function(d, jagsmod, trtnames=NULL, refstudy=NULL,
                        times=seq(1,60, length.out=100)) {
  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }
  if (!all(c("d", "mu") %in% jagsmod$parameters.to.save)) {
    stop("d and mu must be monitored in jagsmod")
  }

  if ("d" %in% names(jagsmod$BUGSoutput$sims.list)) {
    if (dim(jagsmod$BUGSoutput$sims.list$d)[3]==2) {
      polyorder <- 1
    } else if (dim(jagsmod$BUGSoutput$sims.list$d)[3]==3) {
      polyorder <- 2
    }
  } else {
    if (dim(jagsmod$BUGSoutput$sims.list$mu)[3]==2) {
      polyorder <- 1
    } else if (dim(jagsmod$BUGSoutput$sims.list$mu)[3]==3) {
      polyorder <- 2
    }
  }

  # Speeds up computation
  if (jagsmod$BUGSoutput$n.sims > dim(d)[1]) {
    mcmc.index <- 1:dim(d)[1]
  } else if (jagsmod$BUGSoutput$n.sims < dim(d)[1]) {
    mcmc.index <- 1:jagsmod$BUGSoutput$n.sims
  } else {
    mcmc.index <- 1:jagsmod$BUGSoutput$n.sims
  }

  jagsdat <- jagsmod$model$data()
  sims.list <- jagsmod$BUGSoutput$sims.list

  # Get index of reference study
  if (!is.null(attributes(jagsmod)$studynames)) {
    refstudy.ind <- which(attributes(jagsmod)$studynames %in% refstudy)
  } else {
    refstudy.ind <- 1
  }

  #reftrt <- jagsdat$t[refstudy.ind,1]

  haz.df <- data.frame()
  cumhaz.df <- data.frame()
  mort.df <- data.frame()
  S.df <- data.frame()

  mu <- sims.list$mu[, refstudy.ind, ]
  #d  <- sims.list$d

  exponents <- jagsdat$P1
  if (!is.null(jagsdat$P2)) {
    exponents <- c(exponents, jagsdat$P2)
  }

  dt <- diff(c(0,times))

  # Moved outside of the treatments loop
  beta <- mu[mcmc.index,] #+ (d[mcmc.index,k,] - d[mcmc.index,reftrt,])


  loghaz <- get_fp(x = times,
                   params = beta,
                   exponents = exponents)

  #loghaz.d <- array(dim=c(dim(loghaz)[1], dim(d)[2], dim(loghaz)[2]))


  for (k in 1:dim(d)[2]) {

    # Possibly need to ensure same number of samples
    #loghaz.d[,k,] <- loghaz + d[,k]
    loghaz.d <- loghaz + d[mcmc.index,k]

    haz <- exp(loghaz.d)

    dH <- dt * haz # approximate the cumulative hazard (for every MCMC iteration); first calculate the incerments over every interval, then sum up
    H  <- apply(dH, MAR = 2, cumsum)
    #H  <- apply(haz, MAR = 2, cumsum)
    mort <- 1-exp(-H)
    S  <- exp(-H)

    # Create data frames
    S_sum <- t(apply(S, MAR = 1, FUN = BristolTAG:::mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       S_sum)
    S.df <- rbind(S.df, temp)

    mort_sum <- t(apply(mort, MAR = 1, FUN = BristolTAG:::mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       mort_sum)
    mort.df <- rbind(mort.df, temp)

    cumhaz_sum <- t(apply(H, MAR = 1, FUN = BristolTAG:::mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       cumhaz_sum)
    cumhaz.df <- rbind(cumhaz.df, temp)

    haz_sum <- t(apply(haz, MAR = 1, FUN = BristolTAG:::mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       haz_sum)
    haz.df <- rbind(haz.df, temp)
  }

  names(haz.df) <- names(cumhaz.df) <- names(mort.df) <- names(S.df) <-
    c("time", "treatment", "mean", "sd", "2.5%", "50%", "97.5%")

  # Add treatment names
  haz.df$treatment <- trtnames[haz.df$treatment]
  cumhaz.df$treatment <- trtnames[cumhaz.df$treatment]
  mort.df$treatment <- trtnames[mort.df$treatment]
  S.df$treatment <- trtnames[S.df$treatment]

  out <- list(haz=haz.df, cumhaz=cumhaz.df, mort=mort.df, S=S.df,
              P1=jagsdat$P1, P2=jagsdat$P2)

  attr(out, "refstudy") <- refstudy
  attr(out, "trtnames") <- trtnames
  attr(out, "studynames") <- NA
  attr(out, "ipd") <- NA
  class(out) <- "surv.predicts"

  return(out)
}





#' Apply constant HR to baseline hazard from Kaplan-Meier data
#'
#' @param df Is a data frame of IPD containing variables for `time` and `event`
#' @param d A matrix of MCMC HR samples to apply to KM hazard. If there are multiple
#' histology types then it has two columns (and in which case `pS` should also be
#' given a value), otherwise it should have a single column.
#' @param pS A numeric value corresponding to the proportion of squamous patients in
#' each study, used to weight treatment effects. Can be left as `NULL` to ignore
#' histology weighting
#'
consthaz_kmfit <- function(df, d, pS=NULL) {

  if (!is.null(pS)) {
    if (dim(as.matrix(d))[2]!=2) {
      stop("Histology effect assumed: pS has been given a value so dim(d)[2] should be 2")
    }
  }

  KM.est <- survfit(Surv(df$time, df$event) ~ 1, type="kaplan-meier")

  # Sample from KM.est
  n.sims <- 10000
  surv_probs2 <- matrix(nrow=n.sims, ncol=length(KM.est$time)-1)
  km.mat <- matrix(c(log(KM.est$surv), log(KM.est$lower), log(KM.est$upper)), byrow = FALSE, ncol=3)
  for (i in 1:n.sims) {
    kmsamp <- apply(km.mat, MARGIN=1, FUN=function(x) {
      exp(rnorm(1, x[1], (x[3]-x[2])/(2*1.96)))
    }
    )

    surv_probs2[i,] <- addhr_surv(KM.est$time, kmsamp,
                                  d=d[sample(1:nrow(d),1),], pS=pS)
  }

  # Create new KM object
  KM.est2 <- KM.est
  KM.est2$surv <- c(1, apply(surv_probs2, MARGIN=2, median))
  KM.est2$lower <- c(1, apply(surv_probs2, MARGIN=2, FUN=function(x) quantile(x, probs=0.025)))
  KM.est2$upper <- c(1, apply(surv_probs2, MARGIN=2, FUN=function(x) quantile(x, probs=0.975)))



  # # Recalculate survival and 95%CI by applying HR
  # surv_probs2 <- addhr_surv(KM.est$time, KM.est$surv, d=d, pS=pS)
  # lower2 <- addhr_surv(KM.est$time, KM.est$lower, d=d, pS=pS)
  # upper2 <- addhr_surv(KM.est$time, KM.est$upper, d=d, pS=pS)
  #
  # # Create new KM object
  # KM.est2 <- KM.est
  # KM.est2$surv <- c(1,surv_probs2)
  # KM.est2$lower <- c(1,lower2)
  # KM.est2$upper <- c(1,upper2)

  return(list(KM.est, KM.est2))
}




addhr_surv <- function(surv_times, surv_probs, d, pS=NULL) {

  # Calculate the hazard function
  haz1 <- -diff(log(surv_probs)) / diff(surv_times)

  if (!is.null(pS)) {
    d <- (d[1] * (1-pS)) + (d[2] * pS)
  }

  haz2 <- haz1 * exp(d)

  # Cumulative haz and survival for 2nd arm
  cumhaz2 <- cumsum(haz2 * diff(c(surv_times)))
  surv_probs2 <- exp(-cumhaz2)
  return(surv_probs2)
}






comparefits_survplot <- function(jagsmod,
                                 #times=seq(1,60, length.out=100),
                                 d.weights=NULL,
                                 plotinterval=TRUE
) {

  quantity <- "S"
  ipd <- attributes(jagsmod)$ipd

  if (!is.null(d.weights)) {
    if (length(d.weights)!=length(unique(ipd$study))) {
      stop("length(d.weights) must be equal to the number of studies in the dataset")
    }
  }

  # Create data frame of survival predictions for each reference study
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = dplyr::n_distinct(ipd$study), # Maximum value of the progress bar
                              style = 3)   # Character used to create the bar

  plot.df <- data.frame()
  for (s in seq_along(unique(ipd$study))) {
    study <- unique(ipd$study)[s]

    surv <- survcalc(jagsmod, times=seq(1, max(ipd$time[ipd$study==study]), length.out = 100),
                     refstudy=study, d.weight=d.weights[s])

    # Only show survival predictions for treatments within study
    subtrt <- unique(ipd$treatment[ipd$study==study])

    surv <- surv$S %>%
      subset(treatment %in% subtrt) %>%
      dplyr::mutate(refstudy=study)

    plot.df <- rbind(plot.df, surv)

    setTxtProgressBar(pb, s)
  }

  return(plot.df)

  #
  # ##### Create plot #####
  #
  # jagsdat <- jagsmod$model$data()
  # trtnames <- attr(jagsmod, "trtnames")
  # plot.df$treatment <- factor(plot.df$treatment, labels=trtnames, levels=trtnames)
  # plot.df$refstudy <- factor(plot.df$refstudy)
  #
  # # Define colours
  # cols <- RColorBrewer::brewer.pal(dplyr::n_distinct(trtnames), "Set1")
  #
  # # Subset by treatments
  # capt <- paste0("Fractional polynomial; P1 = ", jagsdat$P1, ifelse(!is.null(jagsdat$P2), paste0(", P2 = ", jagsdat$P2), ""))
  #
  # g <- ggplot2::ggplot(plot.df, ggplot2::aes(x=time, ymin=`97.5%`, ymax=`2.5%`, y=`50%`,
  #                                            color=treatment, fill=treatment, linetype=treatment)) +
  #   ggplot2::geom_line() +
  #   ggplot2::facet_wrap(~refstudy) +
  #   ggplot2::xlab("Time") + ggplot2::ylab(quantity) +
  #   ggplot2::theme_bw() +
  #   ggplot2::scale_fill_manual(name = "Treatment", values=cols) +
  #   ggplot2::scale_color_manual(name = "Treatment", values=cols) +
  #   ggplot2::scale_linetype_discrete(name = "Treatment") +
  #   ggplot2::labs(caption=capt)
  #
  # if (plotinterval) {
  #   g <- g + ggplot2::geom_ribbon(alpha=0.3)
  # }
  #
  # return(g)
}





#' ST1-specific version of `BristolTAG::anova_data()`
#'
#' Allows for pdl1cat (multiple PDL1 thresholds) to be aggregated separately
#'
#' @inheritParams BristolTAG::anova_data
#'
anova_data_st1 <- function(timepoints, df){

  if ("pdl1cat" %in% names(df)) {
    pdl1cats <- unique(df$pdl1cat)
  } else {
    pdl1cats <- 1
  }

  dfout <- data.frame()

  for (cat in seq_along(pdl1cats)) {

    if ("pdl1cat" %in% names(df)) {
      dftemp <- subset(df, pdl1cat %in% pdl1cats[cat])
    } else {
      dftemp <- df
    }

    # Split the data at timepoints
    df2 <- survival::survSplit(Surv(time, event) ~., data=dftemp,
                               cut=timepoints[2:length(timepoints)], episode ="timegroup")

    # Calculate offset
    df2$y <- df2$time - df2$tstart

    # Add a variable that equals one for all patients - this is so the number at risk
    # can be calculated when we collapse the data
    df2$n <- 1

    # Collapse data
    df3 <- doBy::summaryBy(y + event + n ~ timegroup + treatment + study, FUN=c(sum, max), data=df2)
    df3 <- subset(df3, select=-c(event.max, n.max))
    names(df3) <- c("spgrp", "treatment", "trialid", "y", "nevents", "natrisk", "y.max")


    # Add in a start time variable
    df3$start <- NA
    for (i in unique(df3$spgrp)) {
      df3$start[df3$spgrp==i] <- timepoints[i]
    }

    # Add in a time variable (i.e. how long since time 0 to max value of y for each row)
    df3$time <- df3$start + df3$y.max

    # Add pdl1cat
    if ("pdl1cat" %in% names(df)) {
      df3$pdl1cat <- pdl1cats[cat]
    }

    dfout <- rbind(dfout, df3)
  }



  # Return the formatted dataset
  out <- list("aggregate"=dfout, "ipd"=df)
  return(out)

}




#' ST1-specific version of `BristolTAG::fp_data()`
#'
#' Allows for data for multiple PDL1 thresholds to be incorporated
#'
#' @inheritParams BristolTAG::fp_data
#'
fp_data_st1 <- function(anova, trtnames) {

  data <- anova$aggregate

  if (any(is.na(match(trtnames, data$treatment))) | any(is.na(match(data$treatment, trtnames)))) {
    stop("trtnames must match the values in data$treatment")
  }

  # Set study and treatment codes
  codes <- trtnames

  data$txCode <- as.numeric(factor(data$treatment, labels = codes,
                                   levels=codes))

  studynames <- unique(data$trialid)
  data$trialid <- as.numeric(factor(data$trialid, labels=studynames, levels=studynames))


  # Order data (using tidyverse)
  if ("pdl1cat" %in% names(data)) {
    data <- dplyr::arrange(data, trialid, pdl1cat, txCode, spgrp)

    #-----------------------------------------------------------------------------
    # Data formatting
    #-----------------------------------------------------------------------------

    # Need to number the treatment arms within each trial
    # data <- data %>%
    #   dplyr::group_by(trialid, pdl1cat, spgrp) %>%
    #   dplyr::mutate(arm=seq(dplyr::n()))
    data <- data %>%
      dplyr::group_by(trialid) %>%
      dplyr::mutate(arm=dense_rank(txCode))


    # Drop uneven spgrp within a trial (test if removing makes a difference by using 1st bit of code)
    # data <- data %>%
    #   dplyr::group_by(trialid, pdl1cat, spgrp) %>%
    #   dplyr::mutate(drop=dplyr::n()) %>%
    #   subset(drop!=1)

    # data <- data %>%
    #   dplyr::group_by(trialid, txCode, arm, pdl1cat) %>%
    #   dplyr::mutate(drop=dplyr::n()) %>%
    #   subset(drop!=1)

    # Maxarm
    # data <- data %>%
    #   dplyr::group_by(trialid, pdl1cat, spgrp) %>%
    #   dplyr::mutate(maxarm=dplyr::n()) %>%
    #   dplyr::ungroup()
    data <- data %>%
      dplyr::group_by(trialid) %>%
      dplyr::mutate(maxarm=max(arm)) %>%
      dplyr::ungroup()

    # Check all arms coded
    #all(data$arm==1 | data$arm==2 | data$arm==3)

    # Length of time intervals
    data$length <- data$time-data$start

    # na
    temp <- data %>% dplyr::select(trialid, maxarm) %>% unique(.)
    na <- temp$maxarm

    # t
    temp <- data %>% dplyr::select(trialid, txCode) %>% unique(.)
    t <- matrix(nrow=length(na), ncol=max(na))
    for (i in seq_along(unique(temp$trialid))) {
      sub <- temp[temp$trialid==unique(temp$trialid)[i],]
      t[i, 1:nrow(sub)] <-
        sub$txCode
    }



  } else {

    stop("Ideally don't use this function for not pdl1cat data")
    # data <- dplyr::arrange(data, trialid, txCode, spgrp)
    #
    # #-----------------------------------------------------------------------------
    # # Data formatting
    # #-----------------------------------------------------------------------------
    #
    # # Need to number the treatment arms within each trial
    # data <- data %>%
    #   dplyr::group_by(trialid, spgrp) %>%
    #   dplyr::mutate(arm=seq(dplyr::n()))
    #
    # # Drop uneven spgrp within a trial (test if removing makes a difference by using 1st bit of code)
    # data <- data %>%
    #   dplyr::group_by(trialid, spgrp) %>%
    #   dplyr::mutate(drop=dplyr::n()) %>%
    #   subset(drop!=1)
    #
    # # data <- data %>%
    # #   group_by(trialid, txCode, arm) %>%
    # #   mutate(drop=n()) %>%
    # #   subset(drop!=1)
    #
    # # Maxarm
    # data <- data %>%
    #   dplyr::group_by(trialid, spgrp) %>%
    #   dplyr::mutate(maxarm=dplyr::n()) %>%
    #   dplyr::ungroup()
    #
    # # Check all arms coded
    # #all(data$arm==1 | data$arm==2 | data$arm==3)
    #
    # # Length of time intervals
    # data$length <- data$time-data$start
    #
    # # na
    # temp <- data %>% dplyr::select(trialid, maxarm) %>% unique(.)
    # na <- temp$maxarm
    #
    # # t
    # temp <- data %>% dplyr::select(trialid, txCode) %>% unique(.)
    # t <- matrix(nrow=length(na), ncol=max(na))
    # for (i in seq_along(unique(temp$trialid))) {
    #   sub <- temp[temp$trialid==unique(temp$trialid)[i],]
    #   t[i, 1:nrow(sub)] <-
    #     sub$txCode
    # }
  }


  # Values for multinomial normal prior (larger than is needed but JAGS code subsets)
  d.mean <- rep(0, 6)
  precarray <- matrix(c(0.0001,0,0,0,0,0,
                           0,0.0001,0,0,0,0,
                           0,0,0.0001,0,0,0,
                           0,0,0,0.0001,0,0,
                           0,0,0,0,0.0001,0,
                           0,0,0,0,0,0.0001
                        ), byrow=TRUE, ncol=6)

  ########### Create JAGS data #########

  jagsdat <- list(s=data$trialid, r=data$nevents, z=data$natrisk, a=data$arm, time=data$time,
                  dt=data$length, N=nrow(data), nt=dplyr::n_distinct(data$treatment), ns=dplyr::n_distinct(data$trialid),
                  mean=d.mean, prec=precarray,
                  npdl=dplyr::n_distinct(data$pdl1cat),
                  pdl=data$pdl1cat,
                  t=t,  na=na)

  # Add ns for distinct and mixed PDL1 thresholds
  if (any(grepl("ZZ", studynames))) {
    len <- sum(grepl("ZZ", studynames))
    jagsdat$ns.pdl <- jagsdat$ns - len
  }

  attributes(jagsdat)$trtnames <- trtnames
  attributes(jagsdat)$studynames <- studynames
  attributes(jagsdat)$ipd <- anova$ipd

  return(jagsdat)
}




#' ST1-specific version of `BristolTAG::fp_geninits()`
#'
#' Allows initial values for multiple PDL1 thresholds to be incorporated
#'
#' @inheritParams BristolTAG::fp_geninits
#'
fp_geninits_st1 <- function(ns, nt, npdl, polyorder=1, seed=890421) {

  if (!is.null(seed)) {
    set.seed <- seed
  }

  init1 <- list()
  init1$d <- array(dim=c(nt, npdl, polyorder+1), 0.1)
  init1$mu <- array(dim=c(ns, npdl, polyorder+1),
                     round(runif(ns*npdl*(polyorder+1), 0.1, 0.6), 1))
  init1$d[1,1:npdl,] <- NA

  init2 <- list()
  init2$d <- array(dim=c(nt, npdl, polyorder+1), 0.2)
  init2$mu <- array(dim=c(ns, npdl, polyorder+1),
                     round(runif(ns*npdl*(polyorder+1), -0.2, 0.7), 1))
  init2$d[1,1:npdl,] <- NA

  init3 <- list()
  init3$d <- array(dim=c(nt, npdl, polyorder+1), round(runif(nt*npdl*(polyorder+1), -0.1,0.1),1))
  init3$mu <- array(dim=c(ns, npdl, polyorder+1),
                     round(runif(ns*npdl*(polyorder+1), -0.1, 0.7), 1))
  init3$d[1,1:npdl,] <- NA

  return(list(init1,
              init2,
              init3))
}


#' Plot multiple survival predictions as different panels
plot_survivalmods <- function(surv.df, all.treats, plot.treats,
                              refstudy) {

  # Define colours
  all.cols <- RColorBrewer::brewer.pal(9, "Set1")
  all.cols <- c(all.cols, RColorBrewer::brewer.pal(4, "Set2"))

  out.df <- surv.df

  out.df$treatment <- factor(out.df$treatment, labels=all.treats, levels=all.treats)

  # Define colours subset by treatments
  cols <- all.cols[unique(as.numeric(out.df$treatment))]

  # Subset by Squamous treatments
  out.df <- subset(out.df, treatment %in% plot.treats)

  out.df$model <- factor(out.df$model)

  capt <-paste("Reference study:", refstudy)

  g <- ggplot2::ggplot(out.df, ggplot2::aes(x=time, ymin=`97.5%`, ymax=`2.5%`, y=`50%`,
                                            color=treatment, fill=treatment, linetype=treatment)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~model) +
    ggplot2::xlab("Time") + ggplot2::ylab("PFS") +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(name = "Treatment", values=cols) +
    ggplot2::scale_color_manual(name = "Treatment", values=cols) +
    ggplot2::scale_linetype_discrete(name = "Treatment") +
    ggplot2::labs(caption=capt)

  return(g)
}




#' Plot implied HRs from multiple models as different panels
plot_hazardmodels <- function(hazardrats, all.treats, reftrt="Docetaxel") {

  # Define colours
  all.cols <- RColorBrewer::brewer.pal(9, "Set1")
  all.cols <- c(all.cols, RColorBrewer::brewer.pal(4, "Set2"))

  out.df <- hazardrats

  out.df$trt1 <- factor(out.df$trt1, labels=all.treats, levels=all.treats)

  out.df <- subset(out.df, trt2==reftrt & trt1 %in% all.treats)
  cols <- all.cols[unique(as.numeric(out.df$trt1))]

  out.df$model <- factor(out.df$model)

  g <- ggplot2::ggplot(out.df, ggplot2::aes(x=time,
                                            #ymin=`97.5%`, ymax=`2.5%`, # Commenting out makes y-axis scales better
                                            y=`50%`,
                                            color=trt1, fill=trt1, linetype=trt1)) +
    ggplot2::geom_line() +
    #ggplot2::geom_ribbon(alpha=0.3) +
    ggplot2::facet_wrap(~model) +
    ggplot2::xlab("Time") +
    ggplot2::ylab("HR") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0("Time-varying relative effects vs treatment ", unique(out.df$trt2))) +
    ggplot2::scale_fill_manual(name = "Treatment", values=cols) +
    ggplot2::scale_color_manual(name = "Treatment", values=cols) +
    ggplot2::scale_linetype_discrete(name = "Treatment")

  return(g)
}





#' Get data frame of proportional hazards model results
get_ph.df <- function(d, times, trt2, trt1) {
  hr.df <- data.frame(trt1=rep(trt1,
                               each=length(times)),
                      trt2=trt2,
                      time=rep(times, length(trt1)),
                      mean=rep(apply(exp(d), MARGIN=2, mean)[2:(length(trt1)+1)],
                               each=length(times)),
                      sd=rep(apply(exp(d), MARGIN=2, sd)[2:(length(trt1)+1)],
                             each=length(times)),
                      l95=rep(apply(exp(d), MARGIN=2, FUN=function(x) quantile(x, probs=0.025))[2:(length(trt1)+1)],
                              each=length(times)),
                      med=rep(apply(exp(d), MARGIN=2, FUN=function(x) quantile(x, probs=0.5))[2:(length(trt1)+1)],
                              each=length(times)),
                      u95=rep(apply(exp(d), MARGIN=2, FUN=function(x) quantile(x, probs=0.975))[2:(length(trt1)+1)],
                              each=length(times))
  ) %>%
    rename("2.5%"=l95,
           "50%"=med,
           "97.5%"=u95)

  return(hr.df)
}




#' ST1-specific version of `BristolTAG::summary.sequence.fpoly()`
#'
#' Gives summary of model fit statistics with residual deviance for each PDL1
#' threshold
#'
#' @inheritParams BristolTAG::summary.sequence.fpoly
#'
sumst1.sequence.fpoly <- function(object, rhat=1.05, pdl.dev=FALSE, ...) {

  if (!class(object) %in% "sequence.fpoly") {
    stop("object must be an object of class 'sequence.fpoly' generated by sequence_fpoly()")
  }

  out.df <- dplyr::tibble("model"=NA, "totresdev"=NA, "pv"=NA, "DIC"=NA, converged=NA, ...)

  if (pdl.dev==TRUE) {
    if (length(dim(object[[1]]$BUGSoutput$sims.list$d))!=4) {
      warning("Dimensions of d do not imply multiple PDL1 effects\nSetting pdl.dev=FALSE")
      pdl.dev <- FALSE
    }
    dev.col <- dim(object[[1]]$BUGSoutput$sims.list$d)[3]
    dev.mat <- matrix(nrow=length(object), ncol=dev.col)
    colnames(dev.mat) <- paste0("resdev", 1:dev.col)
  }

  for (mod in seq_along(object)) {

    # If a model threw an error
    if ("error" %in% names(object[[mod]])) {
      out.df <- out.df %>% dplyr::add_row(model=names(object)[mod],
                                          totresdev=NA,
                                          pv=NA,
                                          DIC=NA,
                                          converged=NA
      )
    } else {

      dat <- object[[mod]]$model$data()

      # Convergence check
      ind <- grep("d\\[", rownames(object[[mod]]$BUGSoutput$summary))
      ind <- append(ind, grep("mu\\[", rownames(object[[mod]]$BUGSoutput$summary)))
      con <- all(object[[mod]]$BUGSoutput$summary[ind, "Rhat"]<rhat)

      out.df <- out.df %>% dplyr::add_row(model=names(object)[mod],
                                          totresdev=round(mean(object[[mod]]$BUGSoutput$mean$totresdev),2),
                                          pv=round(object[[mod]]$BUGSoutput$pD,2),
                                          DIC=round(object[[mod]]$BUGSoutput$DIC,2),
                                          converged=con
      )

      if (pdl.dev==TRUE) {
        for (i in 1:ncol(dev.mat)) {
          temp <- apply(object[[mod]]$BUGSoutput$sims.list$dev[,dat$pdl==i],
                        MARGIN = 1, sum)
          dev.mat[mod, i] <- mean(temp)
        }
      }
    }
  }
  out.df <- out.df[-1,]

  # Add FP powers
  if (class(object) %in% "sequence.fpoly") {
    out.df$P1 <- as.numeric(sapply(out.df$model, FUN=function(x) strsplit(x, "_")[[1]][2]))
    out.df$P2 <- as.numeric(sapply(out.df$model, FUN=function(x) strsplit(x, "_")[[1]][3]))
  }

  if (pdl.dev==TRUE) {
    out.df <- cbind(out.df, dev.mat)

    out.df <- out.df %>%
      select(model, P1, P2,
             totresdev, colnames(dev.mat),
             pv, DIC, converged)

  } else {
    out.df <- out.df %>%
      select(model, P1, P2,
             totresdev,
             pv, DIC, converged)
  }

  return(out.df)
}




plot_st1_survival <- function(modseq.i, times, refstudy, refmod,
                                d, trtnames=c("Docetaxel",
                                              "Pemetrexed",
                                              "Pemetrexed + Carboplatin",
                                              "Selpercatinib"),
                                plot.treats,
                                all.treats=c("Docetaxel", "Atezolizumab", "Pembrolizumab",
                                             "Nivolumab", "Docetaxel-Nintedanib", "Sotorasib",
                                             "Doctxl-PCT", "Pacltxl-Carb", "Pemetrexed",
                                             "Pemetrexed + Carboplatin", "Selpercatinib")) {

  surv <- survcalc(modseq.i,
                   times=times,
                   refstudy=refstudy,
                   refmod=refmod
  )

  surv.hr <- hr_survival(d=d, jagsmod=refmod,
                         times=times,
                         trtnames=trtnames,
                         refstudy=refstudy)

  surv$S <- rbind(surv$S, surv.hr$S)

  # Rearrange so that df is ordered by all.treats
  surv$S <- surv$S %>%
    mutate(temp=factor(treatment, levels=all.treats)) %>%
    arrange(temp, time)

  attributes(surv)$trtnames <- all.treats

  return(surv)

  g <- plot(surv, quantity = "S", plotinterval = FALSE, treats=plot.treats) + geom_line(linewidth=1) +
    ylab("PFS")

  return(g)
}



#' Apply constant HR to a fractional polynomial model to change the network
#' reference
#'
#' Note that this applies a constant HR to the 1st fractional polynomial
#' coefficient (the intercept) to achieve this. The other FP parameters remain
#' relative to the FP model reference.
#'
#' @param d A single column vector containing MCMC samples for logHR of the FP
#' model reference treatment versus the baseline reference treatment
#' @param A character name of the new reference treatment
#' @inheritParams BristolTAG::survcalc
#'
fp_applyhr <- function(d, jagsmod, refmod, reftrt="Selpercatinib") {

  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }
  if (!all(c("d", "mu") %in% jagsmod$parameters.to.save)) {
    stop("d and mu must be monitored in jagsmod")
  }

  adj.fp.d <- jagsmod$BUGSoutput$sims.list$d

  # Ensure row numbers match
  if (nrow(adj.fp.d)>nrow(d)) {
    ind <- c(1:nrow(d), sample(1:nrow(d), nrow(adj.fp.d)-nrow(d)))
    d <- d[ind, 1]
  } else if (nrow(adj.fp.d)<nrow(d)) {
    ind <- sample(1:nrow(d), nrow(adj.fp.d))
    d <- d[ind, 1]
  }

  dimd <- length(dim(jagsmod$BUGSoutput$sims.list$d))

  # For a standard FP model
  if (dimd==3) {

    # Apply constant HR to FP model
    adj.fp.d[,,1] <- apply(adj.fp.d[,,1],
                       MARGIN=2,
                       FUN=function(x) {x - d}) # Assumes consistency

    # Create new array with additional dimension for new reference treatment
    newdim <- dim(adj.fp.d)
    newdim[2] <- newdim[2] +1
    newd <- array(dim=newdim)

    # Add new treatment
    newd[,1,] <- 0
    newd[,2:dim(newd)[2],] <- adj.fp.d

  # For a multiple PDL1 threshold model
  }
  # else if (dimd==4) {
  #   # Apply constant HR to FP model
  #   adj.fp.d[,,,1] <- apply(adj.fp.d[,,,1],
  #                          MARGIN=2,
  #                          FUN=function(x) {x - d}) # Assumes consistency
  #
  #   # Create new array with additional dimension for new reference treatment
  #   newdim <- dim(adj.fp.d)
  #   newdim[2] <- newdim[2] +1
  #   newd <- array(dim=newdim)
  #
  #   # Add new treatment
  #   newd[,1,,] <- 0
  #   newd[,2:dim(newd)[2],,] <- adj.fp.d
  # }
  jagsmod$BUGSoutput$sims.list$d <- newd

  # Apply mu from reference curve and ensure matrix size is the same
  if (nrow(jagsmod$BUGSoutput$sims.list$mu) !=
      nrow(refmod$BUGSoutput$sims.list$mu)) {
    ind <- c(1:dim(refmod$BUGSoutput$sims.list$mu)[1],
             sample(nrow(jagsmod$BUGSoutput$sims.list$mu) -
                      nrow(refmod$BUGSoutput$sims.list$mu),
                    replace=FALSE)
             )
    mu <- refmod$BUGSoutput$sims.list$mu[ind,,]
  } else {
    mu <- refmod$BUGSoutput$sims.list$mu
  }
  jagsmod$BUGSoutput$sims.list$mu <- mu

  attributes(jagsmod)$trtnames <- c(reftrt, attributes(jagsmod)$trtnames)

  return(jagsmod)
}












#' ST1-specific version of `BristolTAG::survcalc()`
#'
#' Allows survcalc to use FP models via a constant HR model to apply to a
#' baseline.
#'
#' @inheritParams BristolTAG::survcalc
#'
survcalc_st1 <- function(jagsmod, refstudy, refmod=jagsmod,
                     times=seq(1,60, length.out=100),
                     n.mcmc=NULL, d.weight=NULL) {

  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }
  if (!"rjags" %in% class(refmod)) {
    stop("refmod must be an object of class rjags")
  }
  if (!all(c("d", "mu") %in% jagsmod$parameters.to.save)) {
    stop("d and mu must be monitored in jagsmod")
  }

  jagsdat <- jagsmod$model$data()
  sims.list <- jagsmod$BUGSoutput$sims.list

  if (is.null(jagsdat$P2)) {
    polyorder <- 1
  } else {
    polyorder <- 2
  }

  # Speeds up computation
  # And ensures the largest MCMC object samples are fully used
  if (!is.null(n.mcmc)) {
    mcmc.index <- sample(1:jagsmod$BUGSoutput$n.sims, size=n.mcmc)
    matsize <- n.mcmc

    mcmc.index.jags <- mcmc.index
    mcmc.index.ref <- mcmc.index

  } else {
    mcmc.index.jags <- 1:jagsmod$BUGSoutput$n.sims
    mcmc.index.ref <- 1:refmod$BUGSoutput$n.sims

    if (jagsmod$BUGSoutput$n.sims > refmod$BUGSoutput$n.sims) {
      mcmc.index.ref <- c(mcmc.index.ref,
                          sample(1:refmod$BUGSoutput$n.sims,
                                 size=jagsmod$BUGSoutput$n.sims - refmod$BUGSoutput$n.sims)
      )
      matsize <- jagsmod$BUGSoutput$n.sims
    } else if (jagsmod$BUGSoutput$n.sims < refmod$BUGSoutput$n.sims) {
      mcmc.index.jags <- c(mcmc.index.jags,
                           sample(1:jagsmod$BUGSoutput$n.sims,
                                  size=refmod$BUGSoutput$n.sims - jagsmod$BUGSoutput$n.sims)
      )
      matsize <- refmod$BUGSoutput$n.sims
    } else {
      matsize <- jagsmod$BUGSoutput$n.sims
    }
    n.sims <- max(jagsmod$BUGSoutput$n.sims, refmod$BUGSoutput$n.sims)
  }

  # Get index of reference study
  if (!identical(jagsmod, refmod)) {
    message("Model for reference curve is different to the treatment effect model")
  }

  refstudy.ind <- which(attr(refmod, "studynames") %in% refstudy)
  if (length(refstudy.ind)==0) {
    if (dim(refmod$BUGSoutput$sims.list$mu)[2]==1) {
      refstudy.ind <- 1
    } else {
      stop("refstudy is not a named study in attr(refmod, 'studynames')")
    }
  }

  reftrt <- jagsdat$t[refstudy.ind,1]

  haz.df <- data.frame()
  cumhaz.df <- data.frame()
  mort.df <- data.frame()
  S.df <- data.frame()

  mu <- refmod$BUGSoutput$sims.list$mu[, refstudy.ind, ]
  d  <- sims.list$d

  # Weight treatment effects
  if (!is.null(d.weight)) {
    d <- weight_shared_d(jagsmod, d.weight=d.weight)
  }

  exponents <- jagsdat$P1
  if (!is.null(jagsdat$P2)) {
    exponents <- c(exponents, jagsdat$P2)
  }

  dt <- diff(c(0,times))

  for (k in 1:jagsdat$nt) {

    beta <- mu[mcmc.index.ref,] +
      (d[mcmc.index.jags,k,] - d[mcmc.index.jags,reftrt,])

    loghaz <- get_fp(x = times,
                     params = beta,
                     exponents = exponents)

    haz <- exp(loghaz)

    dH <- dt * haz # approximate the cumulative hazard (for every MCMC iteration); first calculate the incerments over every interval, then sum up
    H  <- apply(dH, MAR = 2, cumsum)
    #H  <- apply(haz, MAR = 2, cumsum)
    mort <- 1-exp(-H)
    S  <- exp(-H)

    # Create data frames
    S_sum <- t(apply(S, MAR = 1, FUN = mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       S_sum)
    S.df <- rbind(S.df, temp)

    mort_sum <- t(apply(mort, MAR = 1, FUN = mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       mort_sum)
    mort.df <- rbind(mort.df, temp)

    cumhaz_sum <- t(apply(H, MAR = 1, FUN = mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       cumhaz_sum)
    cumhaz.df <- rbind(cumhaz.df, temp)

    haz_sum <- t(apply(haz, MAR = 1, FUN = mcmc_sum))
    temp <- data.frame(time = times,
                       treatment = k,
                       haz_sum)
    haz.df <- rbind(haz.df, temp)
  }

  names(haz.df) <- names(cumhaz.df) <- names(mort.df) <- names(S.df) <-
    c("time", "treatment", "mean", "sd", "2.5%", "50%", "97.5%")

  # Add treatment names
  trtnames <- attr(jagsmod, "trtnames")
  haz.df$treatment <- trtnames[haz.df$treatment]
  cumhaz.df$treatment <- trtnames[cumhaz.df$treatment]
  mort.df$treatment <- trtnames[mort.df$treatment]
  S.df$treatment <- trtnames[S.df$treatment]

  out <- list(haz=haz.df, cumhaz=cumhaz.df, mort=mort.df, S=S.df,
              P1=jagsdat$P1, P2=jagsdat$P2)

  attr(out, "refstudy") <- refstudy
  attr(out, "trtnames") <- trtnames
  attr(out, "studynames") <- attr(jagsmod, "studynames")
  attr(out, "ipd") <- attr(jagsmod, "ipd")
  class(out) <- "surv.predicts"

  return(out)
}










#' ST1-specific version of `BristolTAG::output_coda_fp()`
#'
#' @inheritParams BristolTAG::output_coda_fp
#' @param d.hr A matrix of MCMC logHRs from a constant HR NMA
#' @param d.treatments A character vector of treatment names that match
#' the columns in `d`
#'
#' @export
output_coda_fp_st1 <- function(jagsmod, format="mvn", refstudy=NULL,
                           treatments=attr(jagsmod, "trtnames"),
                           d.hr=NULL, d.treatments=NULL
                           ) {
  if (!"rjags" %in% class(jagsmod)) {
    stop("jagsmod must be an object of class rjags")
  }

  sims.list <- jagsmod$BUGSoutput$sims.list

  trtnames <- attr(jagsmod, "trtnames")
  studynames <- attr(jagsmod, "studynames")

  if (!is.null(refstudy)) {
    mu <- sims.list$mu[,which(studynames %in% refstudy),]

    # Get correct parameter estimates into matrix form
    mu <- apply(mu, MARGIN=2, cbind)
  }

  d <- sims.list$d

  # Loop over FP parameters
  dmat <- matrix(nrow=dim(d)[1])
  cols <- vector()
  for (i in 1:dim(d)[3]) {
    dmat <- cbind(dmat, apply(d[,which(trtnames %in% treatments),i], MARGIN=2, cbind))
    cols <- append(cols, paste(trtnames[which(trtnames %in% treatments)], i, sep="_"))
  }
  dmat <- dmat[,-1]

  colnames(dmat) <- cols

  if (!is.null(refstudy)) {
    colnames(mu) <- paste(refstudy, 1:ncol(mu), sep="_")
    outmat <- cbind(mu, dmat)
  } else {
    outmat <- dmat
  }

  if (format=="mvn") {

    # Means
    param.means <- apply(outmat, MARGIN=2, mean)

    # Covariance - Should reference treatment d's be dropped?
    param.covar <- cov(outmat)

    # Add constant HR NMA parameters (ST1-specific)
    if (!is.null(d.hr)) {
      mean.hr <- apply(d.hr, MARGIN=2, mean, drop=FALSE)
      names(mean.hr) <- d.treatments

      cov.hr <- cov(d.hr)

      # Add constant HR results (covariance with FP parameters is 0)
      ind <- as.numeric(gsub("(.+)(_)(.+)", "\\3", names(param.means)))

      if (d.treatments[1]!=treatments[1]) {
        stop("Reference treatment in jagsmod and d.hr is not the same")
      }
      if (length(d.treatments)!=length(mean.hr)) {
        stop("d.treatments must be the same length as the columns of d.hr")
      }

      d.mean <- mean.hr[2:length(mean.hr)]
      names(d.mean) <- paste(d.treatments[2:length(d.treatments)], 1, sep="_")
      meanpar <- c(param.means[ind==1],
                   d.mean)

      d.cov <- cov.hr[2:nrow(cov.hr), 2:ncol(cov.hr)]
      colnames(d.cov) <- rownames(d.cov) <- d.treatments[2:length(d.treatments)]

      for (k in 2:dim(d)[3]) {
        temp.mean <- rep(0, length(d.mean))
        names(temp.mean) <- paste(d.treatments[2:length(d.treatments)], k, sep="_")
        meanpar <- append(meanpar,
                          c(param.means[ind==k],
                            temp.mean)
        )
      }
      param.means <- meanpar


      # Add covariances
      col <- length(unique(c(treatments, d.treatments))) * dim(d)[3]
      comb.covar <- matrix(ncol=col,nrow=col)
      nams <- unique(c(treatments, d.treatments))
      colnames(comb.covar) <- rownames(comb.covar) <-
        c(paste(nams, 1, sep="_"),
          paste(nams, 2, sep="_"),
          paste(nams, 3, sep="_"))

      for (i in seq_along(rownames(comb.covar))) {
        for (k in seq_along(rownames(comb.covar))) {

          if (rownames(comb.covar)[i] %in% colnames(param.covar)) {
            if (rownames(comb.covar)[k] %in% colnames(param.covar)) {
              comb.covar[i,k] <-
                param.covar[rownames(comb.covar)[i], rownames(comb.covar)[k]]
            } else {
              comb.covar[i,k] <- 0
            }

          } else if (rownames(comb.covar)[k] %in% colnames(param.covar)) {
            if (rownames(comb.covar)[i] %in% colnames(param.covar)) {
              comb.covar[i,k] <-
                param.covar[rownames(comb.covar)[i], rownames(comb.covar)[k]]
            } else {
              comb.covar[i,k] <- 0
            }

          } else {
            if (gsub("(.+)(_)([0-9])", "\\3", rownames(comb.covar)[i]) %in% c("2", "3") |
                gsub("(.+)(_)([0-9])", "\\3", rownames(comb.covar)[k]) %in% c("2", "3")) {
              comb.covar[i,k] <- 0
            } else {
              comb.covar[i,k] <- d.cov[gsub("(.+)(_)([0-9])", "\\1", rownames(comb.covar)[i]),
                                       gsub("(.+)(_)([0-9])", "\\1", rownames(comb.covar)[k])]
            }
          }
        }
      }

      param.covar <- comb.covar
    }

    # Param names
    temp <- paste0("d[", which(trtnames %in% treatments), ",")
    temp <- as.vector(t(sapply(temp, FUN=function(x) {paste0(x, 1:dim(d)[3], "]")})))

    if (!is.null(refstudy)) {
      param.symbol <- c(paste0("mu[", which(studynames %in% refstudy), ",", 1:ncol(mu), "]"),
                        temp)
    } else {
      param.symbol <- temp
    }

    if (is.null(d.hr)) {
      param.symbol <- NULL
    }

    out <- list(mean=param.means, covar=param.covar, params=param.symbol)

  } else if (format=="coda") {

    out <- outmat

  } else {
    stop("format must be 'mvn' or 'coda'")
  }

  return(out)
}





output_coda_ref <- function(fp.mod.1, refstudy) {
  studynames <- attr(fp.mod.1, "studynames")
  sims.list <- fp.mod.1$BUGSoutput$sims.list

  mu <- sims.list$mu[,which(studynames %in% refstudy),]

  # Get correct parameter estimates into matrix form
  mu <- apply(mu, MARGIN=2, cbind)

  colnames(mu) <- paste(refstudy, 1:ncol(mu), sep="_")
  outmat <- mu

  # Means
  param.means <- apply(outmat, MARGIN=2, mean)

  # Covariance - Should reference treatment d's be dropped?
  param.covar <- cov(outmat)

  out <- list(mean=param.means, covar=param.covar)

  return(out)
}
