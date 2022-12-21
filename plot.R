# plot.R
# Plots for Yukon River Chinook run reconstruction
# Steven Rossi (spr1@sfu.ca), Beau Doherty, Sean Cox
# June 3, 2019

plotAll <- function( rpt, folder="." )
{
  plotFitI(rpt=rpt,folder=folder)
  plotTotalRunSize(rpt=rpt,folder=folder)
  plotRunSize(rpt=rpt,folder=folder)
  plotArrival(rpt=rpt,folder=folder)
  plotArrivalByYear(rpt=rpt,folder=folder)
  plotCompResid(rpt=rpt,folder=folder)
}

plotCompResid <- function( rpt, d=16:100, folder="." )
{
  pdf( file=paste(folder,"/compResid.pdf",sep=""), height=11, width=8.5 )

  par( mfrow=c(5,1), mar=c(3,12,3,1) )

  obsN_sdtg <- rpt$n_sdtg[ ,d, , ]
  for( g in 1:rpt$nG )
  {
    for( t in 1:rpt$nT )
    {
      P_sd <- apply(obsN_sdtg[ , ,t,g],2,standardize)
      Phat_sd <- rpt$Phat_sdtg[ ,d,t,g]
      if(sum(!is.na(obsN_sdtg[,,t,g]))>0)
      {
        main <- paste( rpt$gears[g], rpt$years[t] )
        res_sd <- P_sd - Phat_sd 
        colnames(res_sd) <- rpt$day_d[d]
        rownames(res_sd) <- rpt$stocks
        bubblePlot( res=res_sd, main=main )
      }
    }
  }

  dev.off()
}

bubblePlot <- function( res, main="" )
{

  # Add legend symbols
  res[2:7,ncol(res)-5] <- c(-1,-0.5,-0.1,0.1,0.5,1)

  stocks <- rownames(res)
  rownames(res) <- 1:nrow(res)
  days <- as.numeric(colnames(res))
  resDF <- melt(res)
  colnames(resDF) <- c("y","x","value")
  resDF <- as.data.frame(resDF) %>%
           filter( value != 0 ) %>%
           mutate( col    = ifelse( value>0, "black", rgb(191,191,191,maxColorValue=255) ),
                   bgcol  = ifelse( value>0, NA, rgb(191,191,191,50,maxColorValue=255) ),
                   radius = abs(value) )
  
  plot( x=range(days), y=c(0,nrow(res)+1), type="n", xlab="",
        ylab="", axes=FALSE, main=main )
  symbols( x=resDF$x, y=resDF$y, circles=resDF$radius, add=TRUE,
           inches=FALSE, fg=resDF$col, bg=resDF$bgcol )
  axis( side=1 )
  axis( side=2, las=2, labels=stocks, at=1:nrow(res) )
  box()

  # Legend
  xpos <- rep(rev(days)[5],6)
  ypos <- 7:2
  rads <- c(1,0.5,0.1,0.1,0.5,1)
  labs <- rads*c(1,1,1,-1,-1,-1)
  cols <- rep("black",6)
  cols[4:6] <- "grey"
  bgcol <- c(rep(NA,3),
             rep(rgb(191,191,191,50,maxColorValue=255),3))
  text( x=xpos, y=ypos, labels=labs, pos=4 )
  rect( xleft=xpos[1]-4, xright=xpos[1]+5, ybottom=1, ytop=8 )

}

plotFitI <- function( rpt, folder="." )
{
  dims <- list( c(3,4), c(4,6) )
  hei <- c(6,8)
  wid <- c(8,11)
  ylab <- c("Abundance (1000s)","Relative abundance")

  lims <- list( c(18,76), c(15,111) )

  z <- 0
  for( g in 1:rpt$nG )
  {
    pdf( file=paste(folder,"/indexFitsg",g,".pdf",sep=""),
         height=hei[g], width=wid[g] )
    par( mfrow=dims[[g]], mar=c(2,2,1,1), oma=c(0,3,0,0) )
    for( t in 1:rpt$nT )
    {
      I_d <- 1e-3*rpt$E_dtg[ ,t,g]
      E_d <- 1e-3*rpt$Ihat_dtg[ ,t,g]
      if( sum(!is.na(I_d))>0 )
      {
        plot( x=rpt$day_d, y=I_d, xlim=rpt$day_d[lims[[g]]],
              ylim=c(0,1.1*max(I_d,E_d,na.rm=1)),
              las=1, xlab="", ylab="" )
        grid()
        box()
        lines( x=rpt$day_d, y=E_d, lwd=2 )
        legend( x="topright", bty="n", legend=rpt$years[t] )
      }
    } # next t

    mtext( side=2, outer=TRUE, line=1,
           text=ylab[g] )
    dev.off()
  } # next g
}

plotFitIMulti <- function( rptFiles=c("mod1test4/rpt.Rdata",
                                      "mod2test4/rpt.Rdata",
                                      "mod3test4/rpt.Rdata") )
{
  load(rptFiles[1])
  cols <- c("blue","red","green")
  dims <- list( c(4,3), c(6,4) )
  hei <- c(8,9)
  wid <- c(8,8)
  ylab <- c("Abundance (1000s)","Relative abundance")
  ltys <- c(1,2,5)
  lims <- list( c(18,76), c(15,111) )

  z <- 0
  for( g in 1:rpt$nG )
  {
    pdf( file=paste("indexFitsg",g,".pdf",sep=""),
         height=hei[g], width=wid[g] )
    par( mfrow=dims[[g]], mar=c(2,2,1,1), oma=c(2,3,0,0) )
    for( t in 1:rpt$nT )
    {
      E_id <- matrix( data=NA, nrow=3, ncol=rpt$nD )
      for( i in 1:3 )
      {
        load(rptFiles[i])
        E_id[i, ] <- 1e-3*rpt$Ihat_dtg[ ,t,g]
      }
      I_d <- 1e-3*rpt$E_dtg[ ,t,g]
      if( sum(!is.na(I_d))>0 )
      {
        ymax <- 1.1*max(I_d,E_id,na.rm=1)
        plot( x=rpt$day_d, y=I_d, xlim=rpt$day_d[lims[[g]]],
              ylim=c(0,ymax),
              las=1, xlab="", ylab="" )
        plotbg()
        box()
        points( x=rpt$day_d, y=I_d )
        for( i in 1:3 )      
          lines( x=rpt$day_d, y=E_id[i, ], lwd=1.5, col=cols[i], lty=ltys[i] )
        legend( x="topright", bty="n", legend=rpt$years[t] )

        if( g==1 & rpt$years[t]==2005 )
        {
          legend( x=210, y=0.9*ymax, col=cols, lty=ltys, lwd=1.5, bty="n",
                  legend=c("RR_base","RR_oneCor","RR_fullCor"), cex=0.9 )
        }
        else if( g==2 & rpt$years[t]==2007 )
        {
          plot( x=c(0,1), y=c(0,1), type="n", axes=FALSE, xlab="", ylab="" )
          legend( x="topleft", col=cols, lty=ltys, lwd=1., bty="n",
                  legend=c("RR_base","RR_oneCor","RR_fullCor"), cex=1 )
        }

      }
    } # next t

    mtext( side=1, outer=TRUE, line=0.5, text="Julian day" )
    mtext( side=2, outer=TRUE, line=1, text=ylab[g] )
    dev.off()
  } # next g
}

plotTotalRunSize <- function( rpt, folder="." )
{
  pdf( file=paste(folder,"/totalRunSize.pdf",sep=""), height=5, width=7 )

  par( mar=c(3,5,1,1) )

  #t <- !is.na(rpt$I_t)
  t <- rep(TRUE,rpt$nT)
  yr  <- rpt$years[t]
  I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_s)
  E_t <- colSums(exp(rpt$lnRunSize_st))[t]*1e-3
  sonarN_t <- colSums(rpt$E_dtg[ ,t,1])*1e-3
  ymax <- max(I_t,E_t,sonarN_t,na.rm=TRUE)
  
  plot( x=yr, y=I_t, type="n", las=1, yaxs="i", xlab="",
        ylab="Total border passage (1000s)", ylim=c(0,1.1*ymax) )
  grid()
  box()
  
  if( is.finite(rpt$sdrpt[1,5]) )
  {
    Ese <- filter(rpt$sdrpt,par=="runSize_t")[t, ]
    segments( x0=yr, y0=Ese$lCI*1e-3, y1=Ese$uCI*1e-3, col="grey70", lwd=4 )
  }
  
  points( x=yr, y=E_t, pch=16, col="grey40" )
  points( x=yr, y=I_t, pch=0, lwd=1.5 )
  points( x=yr, y=sonarN_t, pch=1, lwd=1.5, col="red" )

  legend( x="bottomleft", bty="n",
          legend=c("Run reconstruction estimates","Mark-recapture estimates","Sonar counts"),
          pch=c(NA,NA), lwd=c(0,4), col=c("grey70","black","red"), lty=c(1,0,0) )
  legend( x="bottomleft", bty="n",
          legend=c("Run reconstruction estimates","Mark-recapture estimates","Sonar counts"),
          pch=c(16,0,1), lwd=c(1.5,1.5), col=c("grey40","black","red"), lty=c(0,0,0) )

  dev.off()
}

plotTotalRunSizeMulti <- function( rptFiles=c("mod1test2/rpt.Rdata",
                                              "mod2test2/rpt.Rdata",
                                              "mod3test2/rpt.Rdata") )
{
  load(rptFiles[1])

  pdf( file="totalRunSize.pdf", height=5, width=7 )

  cols <- brewer.pal(9,"Blues")[c(4,5,7)]

  par( mar=c(3,5,1,1) )

  #t <- !is.na(rpt$I_t)
  t <- rep(TRUE,rpt$nT)
  yr  <- rpt$years[t]
  I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_s)
  sonarN_t <- colSums(rpt$E_dtg[ ,t,1])*1e-3
  ymax <- max(I_t,sonarN_t,na.rm=TRUE)

  plot( x=yr, y=I_t, type="n", las=1, yaxs="i", xlab="",
        ylab="Total run size (1000s)", ylim=c(0,120) )
  grid()
  box()
  
  jtr <- c(-0.2,0,0.2)

  for( i in 1:3 )
  {
    load(rptFiles[i])
    E_t <- colSums(exp(rpt$lnRunSize_st))[t]*1e-3
    Ese <- filter(rpt$sdrpt,par=="runSize_t")[t, ]
    segments( x0=yr+jtr[i], y0=Ese$lCI*1e-3, y1=Ese$uCI*1e-3, col=cols[i], lwd=2 )
    points( x=yr+jtr[i], y=E_t, col="black", pch=14+i, cex=0.5 )
  }


  points( x=yr, y=I_t, pch=0, lwd=1.5 )
  points( x=yr, y=sonarN_t, pch=1, lwd=1.5, col="red" )

  legend( x="bottomleft", bty="n",
          legend=c("Run reconstruction estimates","Mark-recapture estimates","Sonar counts"),
          pch=c(NA,NA), lwd=c(0,4), col=c("grey70","black","red"), lty=c(1,0,0) )
  legend( x="bottomleft", bty="n",
          legend=c("Run reconstruction estimates","Mark-recapture estimates","Sonar counts"),
          pch=c(16,0,1), lwd=c(1.5,1.5), col=c("grey40","black","red"), lty=c(0,0,0) )

  dev.off()
}

plotFig6 <- function( rptFiles=c("mod1test2/rpt.Rdata",
                                 "mod2test2/rpt.Rdata",
                                 "mod3test2/rpt.Rdata") )
{
  load(rptFiles[1])

  pdf( file="figure6.pdf", height=6, width=5 )

  cols <- brewer.pal(9,"Blues")[c(3,5,7)]

  par( mfrow=c(3,1), mar=c(2,3,1,1), oma=c(0,2,0,0) )

  t <- !is.na(rpt$I_t)
  nT <- sum(t)
  yr  <- rpt$years[t]
  I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_s)

  Elow_it <- matrix( data=NA, nrow=3, ncol=nT )
  Eupp_it <- matrix( data=NA, nrow=3, ncol=nT )
  Emle_it <- matrix( data=NA, nrow=3, ncol=nT )

  legs <- c( expression(alpha~"= 10"),
             expression(alpha~"= 150"),
             expression(alpha~"= 500") )

  for( i in 1:3 )
  {
    load(rptFiles[i])
    E_t <- colSums(exp(rpt$lnRunSize_st))[t]*1e-3
    Ese <- filter(rpt$sdrpt,par=="runSize_t")[t, ]
    Elow_it[i, ] <- Ese$lCI*1e-3
    Eupp_it[i, ] <- Ese$uCI*1e-3
    Emle_it[i, ] <- E_t

    plot( x=yr, y=I_t, type="n", las=1, yaxs="i", xlab="",
          ylab="", ylim=c(0,150) )
    grid()
    box()
    segments( x0=yr, y0=Elow_it[i, ], y1=Eupp_it[i, ], col="grey60", lwd=2.5 )
    points( x=yr, y=Emle_it[i, ], col="black", pch=16 )

    if( i==1 )
      points( x=2001, y=142, pch=2, lwd=1.5, col="red" )

    points( x=yr, y=I_t, pch=1, lwd=1.5 )

    legend( x="topleft", legend=legs[i], bty="n", cex=1.4 )

#    legend( x="bottomleft", bty="n",
#            legend=c("Run reconstruction","Mark-recapture"),
#            pch=c(NA,NA), lwd=c(2.5,0), col=c("grey60","black"), lty=c(1,0) )
#    legend( x="bottomleft", bty="n",
#            legend=c("Run reconstruction","Mark-recapture"),
#            pch=c(16,1), lwd=c(1,1.5), col=c("black","black"), lty=c(0,0) )
    if( i>1 )
      axis( side=3, labels=NA )
  }

  mtext( side=2, outer=TRUE, line=0, text="Total border passage (1000s)" )

  dev.off()
}


plotRunSize <- function( rpt, folder="." )
{
  x <- rpt$years
  runSize_st <- exp(rpt$lnRunSize_st)*1e-3

  if( is.finite(rpt$sdrpt[1,5]) )
  {
    Rse <- filter(rpt$sdrpt,par=="runSize_st")
    Rlow_st <- matrix( data=Rse$lCI*1e-3, nrow=rpt$nS, ncol=rpt$nT )
    Rupp_st <- matrix( data=Rse$uCI*1e-3, nrow=rpt$nS, ncol=rpt$nT )

    Rse_t <- filter(rpt$sdrpt,par=="runSize_t")
    Rlow_t <- Rse_t$lCI*1e-3
    Rupp_t <- Rse_t$uCI*1e-3
  }
  else
  {
    Rlow_st <- matrix( data=NA, nrow=rpt$nS, ncol=rpt$nT )
    Rupp_st <- matrix( data=NA, nrow=rpt$nS, ncol=rpt$nT )
    Rlow_t <- rep(NA,rpt$nT)
    Rupp_t <- rep(NA,rpt$nT)
  }

  pdf( file=paste(folder,"/runSize.pdf",sep=""), height=6.5, width=8 )
  par( mfrow=c(3,3), mar=c(3,2,2,1), oma=c(0,2,0,0) )

  for( s in 1:rpt$nS )
  {
    plot( x=range(x), y=c(0,1.1*max(runSize_st[s,],Rupp_st[s,],na.rm=TRUE)),
          type="n", las=1, main=rpt$stocks[s], yaxs="i" )
    grid()
    box()
    segments( x0=x, y0=Rlow_st[s, ], y1=Rupp_st[s, ], col="grey70", lwd=4 )
    points( x=x, y=runSize_st[s, ], pch=16 )
  }

  runSize_t <- colSums(runSize_st)

  plot( x=range(x), y=c(0,1.1*max(runSize_t,Rupp_t,na.rm=1)),
        type="n", las=1, main="Total", yaxs="i" )
  grid()
  box()
  segments( x0=x, y0=Rlow_t, y1=Rupp_t, col="grey70", lwd=4 )
  points( x=x, y=runSize_t, pch=16 )

  mtext( side=2, text="Run size ('000s)", outer=TRUE, line=0.5 )

  dev.off()

}

plotRunSizeMulti <- function( rptFiles=c("mod1/rpt.Rdata",
                                         "mod2/rpt.Rdata",
                                         "mod3/rpt.Rdata")
                            )
{
  stks <- c("L.Mstem","W.Donjek","Pelly","Stewart","Carmacks","Teslin","M.Mstem","U.Mstem")
  cols <- brewer.pal(3,"Set1")[c(1,3,2)]
  #col2 <- c(rgb(225,31,39,150,maxColorValue=255),
  #          rgb(59,127,162,150,maxColorValue=255),
  #          rgb(81,174,79,150,maxColorValue=255))
  #cols <- brewer.pal(3,"Dark2")
  #col2 <- c(rgb(38,157,120,150,maxColorValue=255),
  #          rgb(215,95,28,150,maxColorValue=255),
  #          rgb(117,113,177,150,maxColorValue=255))

  load(rptFiles[1])

  x <- rpt$years

  Rlow_ist <- array( data=NA, dim=c(3,rpt$nS,rpt$nT) )
  Rmle_ist <- array( data=NA, dim=c(3,rpt$nS,rpt$nT) )
  Rupp_ist <- array( data=NA, dim=c(3,rpt$nS,rpt$nT) )
  for( i in 1:3 )
  {
    load(rptFiles[i])
    Rmle_ist[i, , ] <- exp(rpt$lnRunSize_st)*1e-3
    Rse <- filter(rpt$sdrpt,par=="runSize_st")
    Rlow_ist[i, , ] <- matrix( data=Rse$lCI*1e-3, nrow=rpt$nS, ncol=rpt$nT )
    Rupp_ist[i, , ] <- matrix( data=Rse$uCI*1e-3, nrow=rpt$nS, ncol=rpt$nT )
  }

  jtr <- c(-0.2,0,0.2)

  pdf( file="runSizeMult.pdf", height=8, width=7 )
  par( mfrow=c(8,1), mar=c(0,2,0,1), oma=c(2,3,2,2) )

  for( s in 1:rpt$nS )
  {
    ymax <- 1.2*max(Rupp_ist[ ,s,],na.rm=TRUE)
    plot( x=range(x), y=c(0,ymax),
          type="n", axes=FALSE, yaxs="i" )
    grid()
    box()
    legend( x=1983, y=1.05*ymax, bty="n", legend=stks[s], cex=1.5 )
    if( s==1 )
    {
      legend( x="topright", bty="n", col=cols, pch=NA, lwd=2.5,
              legend=c("RR_base","RR_oneCor","RR_fullCor") )
      legend( x="topright", bty="n", col="black", pch=15:17, lwd=0,
              legend=c("RR_base","RR_oneCor","RR_fullCor") )
      axis( side=3 )
    }
    if( s == rpt$nS )
      axis( side=1 )
    if( (s %% 2)==1 )
      axis( side=2, las=1 )
    else
      axis( side=4, las=1 )

    for( t in 1:rpt$nT )
      segments( x0=x[t]+jtr[-3], x1=x[t]+jtr[-1],
                y0=Rlow_ist[-3,s,t], y1=Rupp_ist[-1,s,t],
                col="grey" )

    for( i in 1:3 )
    {
      segments( x0=x+jtr[i],
                y0=Rlow_ist[i,s, ], y1=Rupp_ist[i,s, ],
                col=cols[i], lwd=2.5 )
      points( x=x+jtr[i], y=Rmle_ist[i,s, ], col="black", pch=14+i, cex=0.5 )
    }
  }

  mtext( side=2, text="Border passage (1000s)", outer=TRUE, line=1 )

  dev.off()

}

plotArrivalMulti <- function( rptFiles=c("mod1test2/rpt.Rdata",
                                         "mod2test2/rpt.Rdata",
                                         "mod3test2/rpt.Rdata")
                            )
{
  cols <- brewer.pal(3,"Dark2")
  col2 <- c(rgb(38,157,120,150,maxColorValue=255),
            rgb(215,95,28,150,maxColorValue=255),
            rgb(117,113,177,150,maxColorValue=255))

  load(rptFiles[1])

  x <- rpt$years

  Rlow_ist <- array( data=NA, dim=c(3,rpt$nS,rpt$nT) )
  Rmle_ist <- array( data=NA, dim=c(3,rpt$nS,rpt$nT) )
  Rupp_ist <- array( data=NA, dim=c(3,rpt$nS,rpt$nT) )
  for( i in 1:3 )
  {
    load(rptFiles[i])
    Rmle_ist[i, , ] <- rpt$mu_st
    Rse <- filter(rpt$sdrpt,par=="mu_st")
    Rlow_ist[i, , ] <- matrix( data=Rse$lCI, nrow=rpt$nS, ncol=rpt$nT )
    Rupp_ist[i, , ] <- matrix( data=Rse$uCI, nrow=rpt$nS, ncol=rpt$nT )
  }

  jtr <- c(-0.2,0,0.2)

  pdf( file="arrivalMult.pdf", height=6, width=6.5 )
  par( mfrow=c(4,2), mar=c(0,2,0,1), oma=c(4,2,1,0) )

  for( s in 1:rpt$nS )
  {
    plot( x=range(x), y=c(0.98*min(Rlow_ist[ ,s, ]),1.02*max(Rupp_ist[ ,s,],na.rm=TRUE)),
          type="n", yaxs="i", axes=FALSE )
    legend( x="topright", legend=rpt$stocks[s], bty="n" )
    if( s > 6 )
      axis( side=1 )
    axis( side=2, las=1 )
    grid()
    box()
    for( i in 1:3 )
    {
      segments( x0=x+jtr[i],
                y0=Rlow_ist[i,s, ], y1=Rupp_ist[i,s, ],
                col=col2[i], lwd=1.5 )
      points( x=x+jtr[i], y=Rmle_ist[i,s, ], col=cols[i], pch=16, cex=0.5 )
    }
  }

  dev.off()

}

plotArrival <- function( rpt, folder="." )
{
  cols <- matlab.like2(rpt$nT)
  d <- rpt$day_d

  pdf( file=paste(folder,"/arrivalTiming.pdf",sep=""), height=9, width=9 )
  par( mfrow=c(ceiling(rpt$nS/2),2), mar=c(2,4,1,1), oma=c(3,2,0,0) )

  for( s in 1:rpt$nS )
  {
    rho_dt <- rpt$rho_dst[,s,]
    maxy <- max(rho_dt)
    plot( x=range(d), y=c(0,1.15*maxy), type="n", xlab="",
          ylab="", las=1 )
    plotbg()
    box()
    for( t in 1:rpt$nT )
      lines( x=d, y=rho_dt[ ,t], col=cols[t], lwd=1 )

    #if( s==1 )
    {
      y0 <- 0.15*maxy
      y1 <- maxy
      tseq <- seq(1,rpt$nT,by=10)
      yseq <- seq(y0,y1,length=rpt$nT)[tseq]
      legend_image <- as.raster(matrix(cols, ncol=1))
      rasterImage( image=legend_image,
                   xleft=rev(d)[16],
                   xright=rev(d)[12],
                   ybottom=y0,
                   ytop=y1 )
      text( x      = rev(d)[5],
            y      = yseq,
            labels = rpt$years[tseq] )
      legend( x="topleft", legend=rpt$stocks[s], bty="n" )

                   
    }

  }

  mtext( side=1, text="Julian day", outer=TRUE, line=1, cex=1.3 )
  mtext( side=2, text="Daily border passage proportions", outer=TRUE, line=0, cex=1.3 )

  dev.off()
}

plotArrivalByYear <- function( rpt, folder="." )
{
  controlTable  <- .readParFile( "estControlFile.txt" )
  ctrl <- .createList( controlTable )
  stocks <- ctrl$stks

  cols <- brewer.pal(rpt$nS,"Set1")
  d <- rpt$day_d
  i <- 1

  for( t in 1:rpt$nT )
  {
    if( t %in% c(1,17) )
    {
      pdf( file=paste(folder,"/arrivalTimingByYear",i,".pdf",sep=""), height=7, width=9 )
      par( mfrow=c(4,4), mar=c(2,2,0.5,0.5), oma=c(3,3,0,0) )
    }

    N_ds <- rpt$N_dst[ , ,t]*1e-3
    maxy <- max(N_ds)
    plot( x=range(d), y=c(0,1.15*maxy), type="n", xlab="",
          ylab="", las=1 )
    plotbg()
    box()
    for( s in 1:rpt$nS )
      lines( x=d, y=N_ds[ ,s], col=cols[s], lwd=1.5 )

    legend( x="topleft", bty="n", legend=rpt$years[t] )
    legend( x="topright", bty="n", col=cols, lwd=1.5, legend=stocks, cex=0.6 )

    if( t %in% c(rpt$nT,rpt$nT/2) )
    {
      mtext( side=2, text="Daily border passage (1000s)", outer=TRUE, line=1, cex=1.3 )
      mtext( side=1, text="Julian day", outer=TRUE, line=1, cex=1.3 )
    }
  
    if( t %in% c(16,rpt$nT) )
    {
      dev.off()
      i <- i + 1
    }

  }
}

plotCompN <- function( rpt, folder="." )
{
  n_dtg <- apply(rpt$n_sdtg,2:4,sum)
  d <- rpt$day_d

  for( g in 1:rpt$nG )
  {
    pdf( file=paste(folder,"/sampleSize-",rpt$gears[g],".pdf",sep=""),
         height=6, width=8 )

    ymax <- max(n_dtg[ , ,g],na.rm=TRUE)

    if( g==1 )
      par( mfrow=c(3,3), mar=c(0,0,0,0), oma=c(2,5,1,1) )
    else
      par( mfrow=c(4,6), mar=c(0,0,0,0), oma=c(2,5,1,1) )

    for( t in 1:rpt$nT )
    {
      if( sum(!is.na(n_dtg[ ,t,g]))>0 )
      {
        plot( x=range(d), y=c(0,1.05*ymax), type="n", axes=FALSE, yaxs="i" )
        axis( side=2, las=1 )
        grid()
        box()
        rect( xleft=d-0.3, xright=d+0.3, ybottom=0, ytop=n_dtg[ ,t,g] )
        legend( x="topleft", bty="n", legend=rpt$years[t], cex=1.3 )
      }
    }
  
    mtext(side=2,line=2.5,cex=1.5,outer=TRUE,
          text=paste(rpt$gears[g],"GSI sample size"))

    dev.off()
  }

}

plotQuants <- function(y_qs,ylim=NULL)
{
  if(is.null(ylim))
    ylim <- range(y_qs)
  clrs <- brewer.pal(ncol(y_qs),'Dark2')
  plot( x=c(0.5,ncol(y_qs)+0.5), y=ylim, axes=FALSE,
        las=1, type='n',xlab='', ylab='MRE' )
  grid()
  box()
  abline( h=0, lty=3 )
  segments( x0=1:ncol(y_qs), y0=y_qs[1, ], y1=y_qs[3, ],
            lwd=1.5, col=clrs )
  points( x=1:ncol(y_qs), y=y_qs[2, ],
          pch=16, col=clrs, cex=1, lwd=1.5 )
}

plotStats <- function()
{
  ests <- c("mod1","mod2","mod3")
  sims <- c("sim_OM_base","sim_OM_incFW","sim_OM_incS","sim_OM_incFWS")
  simName <- c("OM_base","OM_incFW","OM_incS","OM_incFWS")

  controlTable  <- .readParFile( "estControlFile.txt" )
  ctrl <- .createList( controlTable )
  stocks <- ctrl$stks

  nM <- length(sims)
  nE <- length(ests)
  nS <- length(stocks)

  # meqs - simulator, estimator, quantile, stock
  MRE_meqs <- array( data=NA, dim=c(nM,nE,3,nS) )
  cvg_me <- array( data=NA, dim=c(nM,nE) )
  med_me <- array( data=NA, dim=c(nM,nE) )

  for( m in 1:nM )
    for( e in 1:nE )
    {
      load( paste("simOutputs",sims[m],ests[e],"perf.Rdata",sep="/") )
      MRE_meqs[m,e, , ] <- apply( perf$stats$MRE_is, 2, quantile, c(0.05,0.5,0.95) )
      cvg_me[m,e] <- length(perf$reps)
      med_me[m,e] <- mean(perf$stats$MRE_is)
    }

  par( mfcol=c(nE,nM), mar=c(0,0,0,0), oma=c(6,5,2,8) )

  for( m in 1:nM )
    for( e in 1:nE )
    {
      plotQuants(MRE_meqs[m,e, , ],ylim=range(MRE_meqs))
      #legend("bottomright",legend=sims[m],bty="n",cex=1.5)
      if( m==1 )
        axis( side=2, las=1 )
      if( m==nM )
      {
        #axis( side=4, las=1 )
        mtext( side=4, text=c("RR_base","RR_oneCor","RR_fullCor")[e], line=.5, cex=1., las=1 )
      }
      if( e==1 )
        mtext(side=3,text=simName[m],line=0.2)
      if( e==nE )
        axis( side=1, at=1:nS, labels=stocks, las=2 )

      #legend("topright",bty="n",legend=cvg_me[m,e])
      #legend("topleft",bty="n",legend=med_me[m,e])

    }

    mtext( side=2, text="Median relative error", line=3, cex=1.5, outer=TRUE )
}

plotStatsVar <- function()
{
  ests <- c("mod1-uncor","mod2-singleCor","mod3-fullCor")
  sims <- c("sim_OM_base","sim_OM_incFW","sim_OM_incS","sim_OM_incFWS")

  controlTable  <- .readParFile( "estControlFile.txt" )
  ctrl <- .createList( controlTable )
  stocks <- ctrl$stks

  nM <- length(sims)
  nE <- length(ests)
  nS <- length(stocks)

  # meqs - simulator, estimator, quantile, stock
  var_meqs <- array( data=NA, dim=c(nM,nE,3,nS) )
  cvg_me <- array( data=NA, dim=c(nM,nE) )

  for( m in 1:nM )
    for( e in 1:nE )
    {
      load( paste("simOutputs",sims[m],ests[e],"perf.Rdata",sep="/") )
      var_meqs[m,e, , ] <- apply( perf$stats$CV_is, 2, quantile, c(0.05,0.5,0.95) )
      cvg_me[m,e] <- length(perf$reps)
    }

  par( mfcol=c(nE,nM), mar=c(0,0,0,0), oma=c(6,6,2,6) )

  for( m in 1:nM )
    for( e in 1:nE )
    {
      plotQuants(var_meqs[m,e, , ],ylim=range(var_meqs))
      #legend("bottomright",legend=sims[m],bty="n",cex=1.5)
      if( m==1 )
        axis( side=2, las=1 )
      if( m==nM )
      {
        #axis( side=4, las=1 )
        mtext( side=4, text=c("base","oneCor","fullCor")[e], line=.5, cex=1.2, las=1 )
      }
      if( e==1 )
        mtext(side=3,text=sims[m],line=0.2)
      if( e==nE )
        axis( side=1, at=1:nS, labels=stocks, las=2 )

      #legend("topright",bty="n",legend=cvg_me[m,e])

    }

    mtext( side=2, text="Residual variance", line=3.5, cex=1.5, outer=TRUE )
}

plotStatsMu <- function()
{
  ests <- c("mod1","mod2","mod3")
  sims <- c("sim_OM_base","sim_OM_incFW","sim_OM_incS","sim_OM_incFWS","sim_OM_big")

  controlTable  <- .readParFile( "estControlFile.txt" )
  ctrl <- .createList( controlTable )
  stocks <- ctrl$stks

  nM <- length(sims)
  nE <- length(ests)
  nS <- length(stocks)

  # meqs - simulator, estimator, quantile, stock
  MRE_meqs <- array( data=NA, dim=c(nM,nE,3,nS) )
  cvg_me <- array( data=NA, dim=c(nM,nE) )
  med_me <- array( data=NA, dim=c(nM,nE) )

  for( m in 1:nM )
    for( e in 1:nE )
    {
      load( paste("simOutputs",sims[m],ests[e],"perf.Rdata",sep="/") )
      MRE_meqs[m,e, , ] <- apply( perf$stats$muMRE_is, 2, quantile, c(0.05,0.5,0.95) )
      cvg_me[m,e] <- length(perf$reps)
      med_me[m,e] <- mean(perf$stats$muMRE_is)
    }

  par( mfcol=c(nE,nM), mar=c(0,0,0,0), oma=c(6,5,2,7) )

  for( m in 1:nM )
    for( e in 1:nE )
    {
      plotQuants(MRE_meqs[m,e, , ],ylim=range(MRE_meqs))
      #legend("bottomright",legend=sims[m],bty="n",cex=1.5)
      if( m==1 )
        axis( side=2, las=1 )
      if( m==nM )
      {
        axis( side=4, las=1 )
        mtext( side=4, text=c("base","cor","RW")[e], line=3.5, cex=1.2, las=1 )
      }
      if( e==1 )
        mtext(side=3,text=sims[m],line=0.2)
      if( e==nE )
        axis( side=1, at=1:nS, labels=stocks, las=2 )

      #legend("topright",bty="n",legend=cvg_me[m,e])
      #legend("topleft",bty="n",legend=med_me[m,e])

    }

    mtext( side=2, text="Relative error", line=3, cex=1.5, outer=TRUE )
}

plotDevSD <- function()
{
  controlTable  <- .readParFile( "estControlFile.txt" )
  ctrl <- .createList( controlTable )
  stocks <- ctrl$stks
  cols <- brewer.pal(3,"Set1")
  ests <- c("mod1test4","mod2test4","mod3test4")
  low_is <- matrix( data=0, nrow=length(ests), ncol=8 )
  mle_is <- matrix( data=0, nrow=length(ests), ncol=8 )
  upp_is <- matrix( data=0, nrow=length(ests), ncol=8 )
  for( i in 1:length(ests) )
  {
    load(paste(ests[i],"/rpt.Rdata",sep=""))
    sdrpt <- filter(rpt$sdrpt,par=="errSD_s")
    low_is[i, ] <- sdrpt$lCI
    mle_is[i, ] <- sdrpt$val
    upp_is[i, ] <- sdrpt$uCI
  }
  x <- 1:8
  par( mar=c(5,6,1,1) )
  plot( x=c(0,9), y=c(min(low_is)*0.95,max(upp_is)*1.05), xaxs="i",
        type="n", yaxs="i", axes=FALSE,
        xlab="", ylab="Process error SD\n" )
  grid()
  box()
#  rect( xleft=x-0.1, xright=x+0.1, ybottom=0, ytop=sd_is[1, ], col=cols[1] )
#  rect( xleft=x-0.3, xright=x-0.1, ybottom=0, ytop=sd_is[2, ], col=cols[2] )
#  rect( xleft=x+0.1, xright=x+0.3, ybottom=0, ytop=sd_is[3, ], col=cols[3] )
  jtr <- c(-0.15,0,0.15)
  for( i in 1:3 )
  {
    segments( x0=1:8+jtr[i], y0=low_is[i, ], y1=upp_is[i, ], col=cols[i], lwd=4 )
    points( x=1:8+jtr[i], y=mle_is[i, ], pch=14+i, cex=0.8 )
  }
  axis( side=1, las=2, at=1:8, labels=stocks, cex.axis=0.8 )
  axis( side=2, las=1 )
  legend( x="topleft", legend=c("RR_base","RR_oneCor","RR_fullCor"),
          bty="n", col=cols, lwd=4, pch=c(22,21,23), pt.bg="black", pt.lwd=0,
          pt.cex=0.8 )

}

plotDailyAvg <- function( rptFile="mod1/rpt.Rdata" )
{
  load( rptFile )
  stks <- c("L.Mstem","W.Donjek","Pelly","Stewart","Carmacks","Teslin","M.Mstem","U.Mstem")
  par( mfrow=c(4,2), mar=c(2,2,1,1), oma=c(4,4,0,0) )

  for( s in 1:rpt$nS )
  {
    x1_d <- numeric( rpt$nD )
    x2_d <- numeric( rpt$nD )
    for( d in 1:rpt$nD )
    {
      x1_d[d] <- mean(rpt$rho_dst[d,s, ])
      x2_d[d] <- median(rpt$rho_dst[d,s, ])
    }
    x2_d <- x2_d/sum(x2_d)
    plot( x=c(170,270), y=c(0,1.05*max(x1_d,x2_d)), type="n", las=1 )
    grid()
    box()
    lines( x=rpt$day_d, y=x1_d, lwd=1.5 )
    lines( x=rpt$day_d, y=x2_d, lwd=1.5, col="red" )
    
    legend( x="topright", legend=stks[s], bty="n", cex=1.2 )

    if( s==1 )
      legend( x="bottomright", legend=c("Mean","Median"), lwd=1.5,
              col=c("black","red"), bty="n" )
  }

  mtext( side=1, text="Julian day", outer=TRUE, line=2, cex=1.2 )
  mtext( side=2, text="Proportion passing border", outer=TRUE, line=2, cex=1.2 )

}









