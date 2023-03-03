#-----------------------------------------------------------------------------#
# fitRR.R                                                                     #
# Estimation functions for Yukon River Chinook run reconstruction             #
#                                                                             #
# Copyright 2019 by Landmark Fisheries Research, Ltd.                         #
#                                                                             #
# This software is provided to DFO in the hope that it will be                #
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                        #
#                                                                             #
# ALL INTELLECTUAL PROPERTY REMAINS WITH LANDMARK FISHERIES RESEARCH, LTD.    #
# THIS SOFTWARE MAY NOT BE REDISTRIBUTED, SUBLICENCED, COPIED, OR SHARED      #
# OUTSIDE OF ESSA TECHNOLOGIES WITHOUT THE EXPRESS WRITTEN CONSENT OF         #
# LANDMARK FISHERIES RESEARCH, LTD.                                           #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE    #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#-----------------------------------------------------------------------------#


# Fit run reconstruction
fitRR <- function( ctlFile="estControlFile.txt",
                   arrivSD=NULL,
                   folder="test",
                   simData=NULL,
                   saveRun=TRUE,
                   fitRR=1 )
{
  folder <- paste0("fits/",folder)

  # Initialize
  source("initChinook.R")

  # Read in control file
  controlTable  <- .readParFile( ctlFile )
  # Create control list
  ctrl <- .createList( controlTable )
  init  <- ctrl$inits

  # Load abundance indices and stock composition data
  load("data/chinookYkData.Rdata")

  # Dimensions
  years  <- ctrl$initYear:ctrl$lastYear
  days   <- chinookYkData$days
  gears  <- chinookYkData$gears
  stocks <- ctrl$stocks
  ages   <- ctrl$ages
  nT     <- length(years)
  nD     <- length(days)
  nH     <- 2
  nG     <- length(gears)
  nS     <- length(stocks)
  nA     <- length(ages)
  nY     <- nT + nA - 1
  nR     <- 3
  nP     <- length(init$arrivMu_p)

  # Create directory for saving plots and report
  suppressWarnings(dir.create(folder))

  # DATA -------------------------------------------------------------------- #

  # Observed numbers by stock, day, year and gear
  n_pdtg <- chinookYkData$n_sdtg[ , ,as.character(years), ]
  # Abundance indices by day, year and gear
  E_dtg  <- chinookYkData$I_dtg[ ,as.character(years), ]
  # Mark-recapture run size indices by year
  I_t    <- chinookYkData$I_t
  seI_t  <- chinookYkData$seI_t

  # Set NAs before first obs and after last obs to 0
  for( t in 1:nT )
  {
    for( g in 1:nG )
    {
      E <- E_dtg[ ,t,g]
      if( sum(!is.na(E)) > 0 )
      {
        E_dtg[is.na(E),t,g] <- 0
      }
    }
  }

  #n_sdtg[ , ,17,2] <- NA

  # Hamazaki 2018
  fishWheelSel_sa <- rbind( c(1,0.295,0.123,0.118), c(0.36,0.129,0.092,0.045) )

  v_asg <- array( data=1, dim=c(nA,2,2) )
  v_asg[ , ,1] <- fishWheelSel_sa

  z_ast <- array( data=1, dim=c(nA,2,nT) )

  catchSD_rh <- array( data=0.1, dim=c(nR,nH) )
  #catchSD_rh[ 2] <- 0.1

  mrCV_t <- rep(0.06,nT)
  mrCV_t[17] <- 0.03

  # Create TMB data object
  data <- list( obsC_rht   = chinookYkData$obsC_rht,
                catchSD_rh = catchSD_rh,
                xB_jgt     = chinookYkData$xB_jgt,
                xC_arht    = chinookYkData$xC_arht,
                n_pdtg     = chinookYkData$n_pdtg,
                E_dtg      = chinookYkData$E_dtg,
                I_t        = chinookYkData$I_t,
                mrCV_t     = mrCV_t,
                day_d      = days,
                aMin       = 4,
                aMax       = 7,
                z_ast      = z_ast,
                RLM_astm   = chinookYkData$RLM_astm,
                mesh_rht   = chinookYkData$mesh_rht,
                v_asg      = v_asg,
                weightI    = 1,
                etaSD      = 1 )


  # Simulated data
  if( !is.null(simData) )
    data <- simData

  # PARAMETERS -------------------------------------------------------------- #

  # Initial conditions
  initMu_p   <- init$arrivMu_p
  sigma_p    <- init$arrivSD_p
  arrivErr_pt <- matrix( data=0, nrow=nS, ncol=nT-1 )
  qE_pg      <- matrix( data=0.05, nrow=nS, ncol=nG )
  qE_pg[ ,2] <- 1
  qI_p       <- rep(0.5,nS)
  errSD_p    <- init$errSD_p
  obsErrSD_g <- rep(0.1,nG)
  cor_pp     <- matrix( data=0, nrow=nS, ncol=nS )
  diag(cor_pp) <- 1 # Correlation matrix must have 1 on diagonal

  lnDisp_tg  <- matrix( data=log(1e-4), nrow=nT, ncol=nG )
  lnDisp_tg[ ,1] <- log(5e-3)
  lnDisp_tg[17,1] <- log(1.)

  initEta_asy <- array( data=0, dim=c(nA-1,2,nY) )
  initEta_asy[1, , ] <- log(4)
  initEta_asy[2, , ] <- log(3)
  initEta_asy[3, , ] <- log(2)

  lnf_rht <- array(-10,dim=c(nR,nH,nT))
  lnf_rht[ ,1, ] <- c(0.1,0.35,0.15)
  lnf_rht[ ,2, ] <- c(0.4,0.05,0.15)
  lnf_rht[ data$obsC_rht==0 ] <- -10

  # Create TMB parameter object
  pars <- list( lnR_py       = array(9,dim=c(nP,nY)),
                lnf_rht      = lnf_rht,
                lnUMSY_p     = rep(log(0.5),nP),
                lnSMSY_p     = rep(log(5e3),nP),
                lnRecSD_p    = rep(log(0.4),nP),
                logitPhi_p   = rep(0,nP),
                pearPars     = c(-0.547,log(0.622),log(0.204),1.920),
                logitpFem_y  = rep(0,nY),
                lnpFemSD     = log(1),
                initEta_asy  = initEta_asy,
                lnArrivMu_p  = log(initMu_p),
                lnArrivSD_p  = log(sigma_p),
                arrivErr_pt  = arrivErr_pt,
                lnErrSD_p    = log(errSD_p),
                logitCor_pp  = logit(cor_pp,lb=-1,ub=1),
                lnqE_pg      = log(qE_pg),
                lnqI_p       = log(qI_p),
                lnDisp_tg    = lnDisp_tg )


  if(!is.null(arrivSD))
    pars$lnArrivSD_p <- rep(log(arrivSD),nS)

  # MAP --------------------------------------------------------------------- #

  qEmap_pg <- matrix( data=1:(nS*nG), nrow=nS, ncol=nG )
  qEmap_pg[ ,1] <- ctrl$map$qFishWheel_p
  qEmap_pg[ ,2] <- NA

  corMap_ps <- matrix( data=1:(nS*nS), nrow=nS, ncol=nS )
  # Always fix diagonal
  diag(corMap_ps) <- NA
  # Optional - define map for lower triangle
  # Let's fix the lower triangle (except diag) at a single value
  for( s in 2:nS )
  {
    if( ctrl$map$corType=="uncor" )
      corMap_ps[s,1:(s-1)] <- NA
    else if( ctrl$map$corType=="single" )
      corMap_ps[s,1:(s-1)] <- nS^2+1
  }
  # Mirror upper & lower triangle
  corMap_ps <- mirrorMatrix(corMap_ps)

  dispMap <- NA*lnDisp_tg
  #dispMap[17,2] <- 1
  #dispMap[22,2] <- 2

  fmap_rht <- array( data=1:(nR*nH*nT), dim=c(nR,nH,nT) )
  fmap_rht[data$obsC_rht==0] <- NA

  map <- list( lnf_rht     = as.factor(fmap_rht),
               pearPars    = as.factor(NA*pars$pearPars),
               logitPhi_p  = as.factor(NA*pars$logitPhi_p),
               #logitPhi_p  = as.factor(rep(1,nP)),
               lnpFemSD    = factor(NA),
               #lnRecSD_p   = as.factor(rep(1,nP)),
               lnArrivMu_p = as.factor(ctrl$map$arrivMu_p),
               lnArrivSD_p = as.factor(ctrl$map$arrivSD_p),
               #arrivErr_pt = as.factor(NA*pars$arrivErr_pt),
               lnErrSD_p   = as.factor(NA*ctrl$map$errSD_p),
               logitCor_pp = as.factor(NA*corMap_ps),
               lnqE_pg     = as.factor(qEmap_pg),
               lnqI_p      = as.factor(ctrl$map$qI_p),
               lnDisp_tg   = as.factor(dispMap) )


  # BUILD AND OPTIMIZE OBJECTIVE FUNCTION ----------------------------------- #

  # Compile and load objective function
  #compile("yukonChinookIntegrated.cpp")
  #dyn.load(dynlib("yukonChinookIntegrated"))

#  load("/Users/rossi/yrcRR/mod1/rpt.Rdata")
#  npars <- rpt[names(pars)]
#browser()
  #pars$lnRunSize_st <- round(npars$lnRunSize_st)
  

  # Build objective function
  obj <- MakeADFun( data       = data,
                    parameters = pars,
                    map        = map,
                    DLL        = "yukonChinookIntegrated",
                    random     = NULL )

  # Set bounds
  low <- obj$par*0-Inf
  upp <- obj$par*0+Inf
  #low[names(low)=="cor_ss"] <- -1
  #upp[names(upp)=="cor_ss"] <- 1

  # Optimization controlsarrivSD
  optCtrl <- list(  eval.max = ctrl$maxFunEval, 
                    iter.max = ctrl$maxIterations )

  # Optimize
  #sink("sink.txt")
  opt <- try( nlminb( start     = obj$par,
                      objective = obj$fn,
                      gradient  = obj$gr,
                      lower     = low,
                      upper     = upp,
                      control   = optCtrl ) )
  #sink()

  rptFE <- obj$report()

  sdrep <- NULL
  errGrad <- NULL

  if( mode(opt)=="character" )
  {
    rpt <- obj$report()
    rpt$opt$convergence <- 1
  }
  else
  {
    # Retrieve optimized parameters and gradients
    par <- data.frame( par  = names(opt$par),
                       val  = opt$par,
                       grad = as.numeric(obj$gr()) )
  
    # Calculate standard errors via delta method
    sink("sink.txt")
    sdobj <- sdreport( obj )
    sdrpt <- summary( sdobj )
    sink()
    if( mode(sdrpt)!="character" )
    {
      colnames(sdrpt) <- c("val","se")
    
      sdrpt <- as.data.frame(sdrpt) %>%
             mutate( par = rownames(sdrpt),
                     lCI = val - qnorm(.95)*se,
                     uCI = val + qnorm(.95)*se ) %>%
             dplyr::select( par, val, se, lCI, uCI )

      errGrad <- filter( par, par=="arrivErr_st" )
      errGrad_st <- matrix( errGrad$grad, nrow=nS, ncol=nT-1 )
    }


    # Build report object
    rpt <- obj$report()
    rpt$opt    <- opt
    rpt$ages  <- ages
    rpt$years  <- years
    rpt$broodYears <- seq( from=years[1]-max(ages), by=1, length.out=nY)
    rpt$gears  <- gears
    rpt$stocks <- stocks
    rpt$regions <- chinookYkData$regions
    rpt$fisheryType <- chinookYkData$fisheryType
    rpt$meshSizes <- chinookYkData$meshSizes
    rpt$par    <- par
    rpt$sdrpt  <- sdrpt
    rpt$errGrad_st <- errGrad_st

#    nObs <- sum(!is.na(data$n_pdtg[1, , , ])) +
#            sum(!is.na(data$E_dtg)) + 
#            length(data$I_t)
#    nPar <- length(opt$par)
#    nll  <- opt$objective
#    aic  <- 2*nPar + 2*nll + 2*nPar*(nPar+1)/(nObs-nPar-1)
#    rpt$aic <- aic
  



#  library(tmbstan)
#  options(mc.cores = 3)
#  mcinit <- list()
#  for( i in 1:3 )
#    mcinit[[i]] <- rnorm(n=length(obj$par),mean=obj$par,sd=1e-5)
#  
#  fit <- tmbstan( obj = obj,
#                  chains = 3,
#                  iter = 1e3,
#                  init = mcinit )
#  rpt$fit <- fit






    if( saveRun )
    {
      save( rpt, file=paste(folder,"/rpt.Rdata",sep="") )
      plotAll(rpt=rpt,folder=folder)
      system( paste("cp ",ctlFile," ",folder,"/estControlFile.txt",sep="") )
    }

  }

  rpt

}