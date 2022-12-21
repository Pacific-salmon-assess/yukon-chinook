######################################################################
# procData.R
#
# Generate R data file for Yukon River run reconstruction for fitRR.R
#
# June 7, 2019
# Authors: Beau Doherty, Steve Rossi, Sean Cox
# Landmark Fisheries Research
######################################################################

library(dplyr)

# Steps
# 1. Read in and process datasets
# 2. Create R data files with arrays for fitRR.r
# 3. Plots

processData <- function()
{
  # --- 1. Read in and process datasets ----
  
  # -- aggregate run-size counts by day, N_gyd --
  counts <- read.csv("data/border-assessment/borderCounts.csv")
  names(counts)[names(counts)=='count_type'] <- 'gear'
  counts$gear <- as.character(counts$gear)
  
  # -- GSI sub-stock campling by day, n_sgyd --
  gsi <- read.csv('data/gsi/gsiSamplesAllProbs.csv')
  stockID <- read.csv('data/gsi/stockIDs.csv') %>% arrange(plotOrder)
  stockID$stockNum <- stockID$plotOrder
  gsi <- dplyr::left_join(gsi, stockID, by='region')
  
  # change gear name to gillnet or fishWheel
  gsi$gear <- as.character(gsi$gear)
  gsi$gear[gsi$gear=='Fish Wheel']    <- 'fishWheel'
  gsi$gear[gsi$gear=='Test Fishery']    <- 'eagle'
  
  # remove rows with no julian_gsi or errors
  gsi <- subset(gsi, !is.na(julian))
  gsi <- subset(gsi, julian <300 & julian >100)
  
  # assign zeros for NA probabilities
  # gsi$prob[is.na(gsi$prob)] <- 0
  
  # Add julian day adjustment for samples/counts from fishwheel site to scale everything relative to Eagle sonar site, which is about 48 km downstream from fishwheel locations (i.e. approx 1-day travel for Chinook)
  fwDayAdj <- 1
  counts$julian[counts$count_type=='fishWheel'] <- counts$julian[counts$count_type=='fishWheel'] +1
  
  gsi$julian_date[gsi$data_label=='YukonRetro'] <- gsi$julian_date[gsi$data_label=='YukonRetro'] +1
  
  # --- 2. Create R data files with arrays for fitRR.r ---
  
  # Generate array with proportions by stock, day, year, gear
  stockNames <- stockID$stock
  stockRegion <- stockID$region
  gsiGear <- c('eagle','fishWheel')
  fDay <- 160
  lDay <- 285
  days <- fDay:lDay
  fYear <- min(gsi$year)
  lYear <- max(gsi$year)
  yrs <- fYear:lYear
  
  # -- mark recapture counts by year, mrI_t --
  mr <- read.csv('data/border-assessment/YkCk_BorderMR_Data.csv')

  # mrI_t index, total return mark-recapture estimates by year for fishwheel
  mr <- mr[mr$Year %in% yrs, ]
  mr$SE[is.na(mr$SE)] <- mr$Estimate[is.na(mr$SE)]*0.06 # Fill in NA SEs with avg SE
  I_t <- mr$Estimate
  seI_t <- mr$SE
  
  # Generate array for counts by stock, day, year, and gear
  n_sdtg <- array(dim=c(length(stockNames),
              length(fDay:lDay),
              length(fYear:lYear),
              length(gsiGear)))
  
  dimnames(n_sdtg) <- list(stockNames=stockNames,
               julianDay=days,
               year = yrs,
               gears=gsiGear)
  
  # Fill n_sdtg array
  # loop over yrs
  for(t in 1:length(yrs))
  {
    # loop over gears
    for (g in 1:length(gsiGear))
    {
      tmp <- subset(gsi, year==yrs[t] & 
                gear == gsiGear[g] &
                !is.na(julian) &
                !is.na(region) &
                prob>0) 
      
      # if no data for given gear, skip
      if(nrow(tmp)==0)
        n_sdtg[,,t,g] <-NA
  
      # Check if any days duplicated across samples
      smpls <- unique(tmp[,c('julian','sample_num')])
  
      if(any(table(smpls$sample_num)>1))
      { 
  
       errors <- smpls$sample_num[duplicated(smpls$sample_num)]
        for(err in errors)
        {
          errDays <- table(tmp$julian[tmp$sample_num==err]) %>%
              sort(decreasing=TRUE)
       
          tmp$julian[tmp$sample_num==err] <- as.integer(names(errDays)[1])
        } 
      }
  
      gsiDat_gt <- tmp
  
      #loop over days
      for(d in 1:length(days))
      {
        tmp <- subset(gsiDat_gt, year== yrs[t] &
                  gear == gsiGear[g] &
                  julian == days[d] &
                  !is.na(julian) &
                  !is.na(region) &
                  prob>0)
  
  
        # if no data for given day, skip
        if(nrow(tmp)==0)
        {
          n_sdtg[,d,t,g] <- NA
        }
        else  
        {
  
          # Renormalize across samples that sum of probs !=1
          nSmpls <- length(unique(tmp$sample_num))
          if( sum(tmp$prob) != nSmpls)
            for(smp in unique(tmp$sample_num))
            {
              nProbs <- tmp[tmp$sample_num==smp,]
              if(sum(nProbs$prob) != 1)
              {
                normProbs <- nProbs$prob/sum(nProbs$prob)
                tmp$prob[tmp$sample_num==smp] <- normProbs
              }
            }
  
          # calculate proportions for each stocks
          n_s <- dplyr::summarize(group_by(tmp,stockNum),
                   expCounts = sum(prob))
          n_s <- dplyr::left_join(stockID,n_s, by='stockNum')
          n_s$expCounts[is.na(n_s$expCounts)] <-0
          
  
          # Check if sum of probs adds to sample nums
          if(round(sum(n_s$expCounts),10) != nSmpls )
            browser(cat('ERROR: sum of normalized GSI probs != sample size'))
  
          n_sdtg[,d,t,g] <- n_s$expCounts
  
        } 
  
      } 
    }
  }
  
  
  # Generate array for index counts
  idxGear <- c('eagle','fishWheel')
  I_dtg <- array(dim=c( length(fDay:lDay),
              length(fYear:lYear),
              length(idxGear)))
  
  dimnames(I_dtg) <- list( julianDay=days,
               year = yrs,
               gears=idxGear)
  
  
  # loop over yrs
  for(t in 1:length(yrs))
  {
    # loop over gears
    for (g in 1:length(idxGear))
    {
      tmp <- subset(counts, year==yrs[t] & 
              gear == idxGear[g] )
  
  
      
      # if no data for given gear, skip
      if(nrow(tmp)==0)
        next()
      
      #loop over days
      for(d in 1:length(days))
      {
        tmp <- subset(counts, year== yrs[t] &
                  gear == idxGear[g] &
                  julian == days[d])
  
  
        # if no data for given day, skip
        if(nrow(tmp)==0)  
          I_dtg [d,t,g] <- NA
        else
        {
          I_dtg [d,t,g] <- sum(tmp$count, na.rm=T)
        }
          
  
      } 
    }
  }
  
  # save list
  chinookYkData <- list(  I_dtg = I_dtg,
              I_t = I_t,
              seI_t = seI_t,
                n_sdtg = n_sdtg,
                stockNames = stockNames,
              stockRegion = stockRegion,
              gears = idxGear,
              fDay = fDay,
              lDay = lDay,
              days = days,
              fYear = fYear,
              lYear = lYear,
              yrs = yrs )
  
  save(chinookYkData, file='data/chinookYkData.Rdata')
  
  # --- 3. Plots ---
  
  # Aggregate GSI samples for all stocks by year
  countsG <- c('eagle','fishWheel')
  gsiG <- c('eagle','fishWheel')
  # gsiG <- c('fishWheel', 'gillnetFW', 'gillnetEagle')
  yrs <- min(gsi$year):max(gsi$year)
  
  gsi_gy <- array(NA,dim=c(length(gsiG),length(yrs)))
  counts_gy <- array(NA,dim=c(length(countsG),length(yrs)))
  
  dimnames(gsi_gy) <- list(gear = gsiG, years=yrs)
  dimnames(counts_gy) <- list(gear = countsG, years=yrs)
  
  for (y in 1:length(yrs))
  { 
    # aggregate gsi samples by year & gear
    for (g in 1:length(gsiG))
    {
      tmp <- subset(gsi, year==yrs[y] & 
                 gear==gsiG[g] &
                 !is.na(julian) &
                 !is.na(region) &
                 prob>0)
      gsi_gy[g,y] <- length(unique(tmp$sample_num))
    }
  
    # aggregate gsi samples by year & gear
    for (g in 1:length(countsG))
    {
      tmp <- subset(counts, year==yrs[y] & 
                  gear==countsG[g])
      counts_gy[g,y] <- sum(tmp$count, na.rm=T)
    }
  
  }
  
  
  # Check gsi_gy and n_sdtg produce same results
  difg1 <- gsi_gy[1,] -apply(n_sdtg[,,,1],c(3),sum,na.rm=T)
  difg2 <- gsi_gy[2,] -apply(n_sdtg[,,,2],c(3),sum,na.rm=T)
  
  # assign NAs instead of zeros for plotting
  gsi_gy[gsi_gy==0] <- NA
  counts_gy[counts_gy==0] <- NA
  
  
  pdf(file='idxAndGSI.pdf')
  
  clrs <- c('#1b9e77', '#d95f02', '#7570b3', 
        '#e7298a','#e6ab02', '#a6761d')
  par(mfrow=c(3,1), mgp=c(2,.6,0), mar=c(2,4,0.5,0.2),
    tck=-0.03)
  
  yMax <- max(gsi_gy[ which(gsiG=='eagle'),],na.rm=T)
  plot(x=yrs, y=gsi_gy[1,], ylab='GSI samples',
     ylim=c(0,yMax), type='n')
  
  nGear <- length(gsiG)+ length(countsG) +1
  legend('topleft',bty='n',
    legend=c(paste('GSI samples',gsiG),
        'Count Index - Fish Wheek',
        'Eagle sonar estimates',
        'Mark-recapture estimates'),
    pch= c(rep(15,length(gsiG)),rep(16,length(countsG)),17),
    col= clrs[1:nGear], cex=1.2)
  
  # GSI samples
  for (g in 1:length(gsiG))
    points(x=yrs,y=gsi_gy[g,], col=clrs[g], pch=15, cex=1.2)
  # points(x=yrs, y=counts_gy[1,],pch=16, col=clrs[4])
  
  # fishwheel counts
  plot(x=yrs, y=counts_gy[which(countsG=='fishWheel'),], ylab='Counts',
     pch=16, col=clrs[length(gsiG)+1], cex=1.2)
  
  # Sonar & Mark-recapture estimates
  yMax <- max(counts_gy[ which(countsG=='eagle'),],na.rm=T)
  plot(x=yrs, y=counts_gy[which(countsG=='eagle'),], 
     ylab='Run Size Estimates',
     ylim=c(0,yMax), pch=16, 
     col=clrs[length(gsiG)+2], cex=1.2)
  points(x=mr$year, y=mr$mark_recapture, 
       col=clrs[nGear], pch=17, cex=1.2)
  
  dev.off()
  
  # Daily Plots of gsi counts:
  par(mfrow=c(5,7),mgp=c(1,0.5,0),
       tck=-0.01, mar=c(2,2,0,0))
  for (y in 1:length(yrs))
  {
    for (g in 1:length(gsiG))
    {
    n_dt  <- apply(n_sdtg[,,y,g],c(2), sum, na.rm=T)
    
    # Check calcs using raw dats
    n_dt2 <- subset(gsi, year==yrs[y] & gear==gsiG[g])
    
    if(nrow(n_dt2) >0)  
    { 
      n_dt2 <- n_dt2[,c('julian','sample_num')] %>% unique()
      n_dt2 <- table(n_dt2$julian)
  
      plot(x=names(n_dt),y=n_dt, type='h',col='black',
         ylim=c(0,50))
        points(x=names(n_dt2),y=n_dt2,col='red',cex=0.1)
    }    
  
    } 
  } 
  
  
  # # Calculate avg. proportions for initial conditions in model
  # avgProp <- dplyr::summarize(group_by(data,stock),
  #               counts=length(gear))
  # avgProp$prop <- avgProp$counts/sum(avgProp$counts)
  
  # N <- dplyr::summarize(group_by(agg,year,count_type),
  #           returns=sum(count,na.rm=T))
  # sonarN <- subset(N, count_type=='eagle')
  
  # avgProp$estAvgReturns <- avgProp$prop*mean(sonarN$returns)
  
  # # Calculate aggregate returns by year
  # I_tg <- apply(I_dtg,c(2,3),sum, na.rm=T)
  # I2 <- read.table('~/Documents/LANDMARK/2018_YukonChinook/subStockModel/data/Yukon_chin_border_passage_indices.txt')
  
  # plot(x=I2$year, y=I2$mark_recapture/1e3, las=1,
  #   ylim=c(0,65),ylab='Numbers (1000s)', xlab='')
  # points(x=yrs[1:26], y=I_tg[1:26,2]/1e3,col='red')
  # legend('topright',bty='n', col=c('black','red'),
  #     legend=c('mark-recapture FW estimates',
  #       'sum of daily FW counts'),
  #     pch=1)
}






