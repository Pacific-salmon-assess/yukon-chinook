#-----------------------------------------------------------------------------#
# plot.R                                                                      #
# Plotting script for integrated Yukon River Chinook run reconstruction       #
#                                                                             #
# Copyright 2023 by Landmark Fisheries Research, Ltd.                         #
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


plotAll <- function( rpt, folder="." )
{
  plotFitI(rpt=rpt,folder=folder)
  #plotTotalRunSize(rpt=rpt,folder=folder)
  plotRunSize(rpt=rpt,folder=folder)
  plotArrival(rpt=rpt,folder=folder)
  #plotArrivalByYear(rpt=rpt,folder=folder)
  #plotCompResid(rpt=rpt,folder=folder)
  plotCatchFit(rpt=rpt,folder=folder)
  plotCatchAgeFit(rpt=rpt,folder=folder)
  plotCatchAgeFitBubble(rpt=rpt,folder=folder)
  plotBorderAgeFit(rpt=rpt,folder=folder)
  plotBorderAgeFitBubble(rpt=rpt,folder=folder)
  plotBorderPassage(rpt=rpt,folder=folder)
  plotRecruitment(rpt=rpt,folder=folder)
  plotStockRec(rpt=rpt,folder=folder)
  plotSelectivity(rpt=rpt,folder=folder)
  plotSelectivityAvg(rpt=rpt,folder=folder)
  plotSelectivityByLength(rpt=rpt,folder=folder)
  plotpFem(rpt=rpt,folder=folder)
  plotProbReturnAtAge(rpt=rpt,folder=folder)
  plotAvgRunTiming(rpt=rpt,folder=folder)
}

plotAvgRunTiming <- function( rpt, folder="." )
{
  pdf( file=paste(folder,"/fig22-avgRunTiming.pdf",sep=""), height=5, width=7 )
  cols <- brewer.pal(rpt$nP,"Set1")

  x <- rpt$day_d

  y_dpt <- rpt$rho_dpt
  y_dp  <- apply( y_dpt, 1:2, mean )

  par(mar=c(5,7,1,1),oma=c(0,0,0,0))

  plot( x=range(x), c(1,rpt$nP+1.5), type="n", axes=FALSE,
        xlab="Ordinal date", ylab="" )
  plotbg()
  axis( side=1 ) 
  axis( side=2, at=rpt$nP:1, labels=rpt$stocks, las=1, cex.axis=0.6 ) 
  for( p in 1:rpt$nP )
  {
    y <- 1.5*(y_dp[ ,p]/max(y_dp[ ,p]))
    polygon( x=c(x,rev(x)), y=rpt$nP-p+1+c(y,rep(0,length(x))),
             border=NA, col=cols[p] )
  }

  dev.off()
}

plotAnnualRunTiming <- function( rpt, folder="." )
{
  #pdf( file=paste(folder,"/annualRunTiming.pdf",sep=""), height=5, width=7 )
  cols <- brewer.pal(rpt$nP,"Set1")

  mu_pt <- rpt$mu_pt
  #for( p in 1:rpt$nP )
  #  mu_pt[p, ] <- mu_pt[p, ]/mu_pt[p,1]

  #mu_pt <- rpt$arrivErr_pt

  par(mar=c(5,7,1,1),oma=c(0,0,0,0))

  plot( x=range(rpt$years), y=range(mu_pt), type="n", las=1,
        xlab="Year", ylab="Mean date of arrival (ordinal)" )
  plotbg()
  #axis( side=1 ) 
  #axis( side=2, at=rpt$nP:1, labels=rpt$stocks, las=1, cex.axis=0.6 ) 
  for( p in 1:rpt$nP )
  {
    lines( x=rpt$years, y=mu_pt[p, ], lwd=1.5, col=cols[p] )
  }

  #dev.off()
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

bubblePlot2 <- function( res, xlab="", ylab="", xax=TRUE, yax=TRUE,
                         x1=1985, xby=10, yby=1, main="",
                         yaxlab=NULL )
{
  colnames(res) <- 1:ncol(res)
  resDF <- melt(res)
  colnames(resDF) <- c("x","y","value")
  resDF <- as.data.frame(resDF) %>%
           filter( value != 0 ) %>%
           mutate( col    = ifelse( value>0, "black", "grey" ),
                   bg     = ifelse( value>0, NA, scales::alpha("grey",0.6) ),
                   radius = abs(value) )
  plot( x=c(1,nrow(res)), y=c(0,ncol(res)+1), type="n", xlab=xlab, ylab=ylab, axes=FALSE, main=main )
  if(xax)
    axis( side=1, at=seq(1,nrow(res),by=xby), labels=seq(x1,x1-1+nrow(res),by=xby) )
  axis( side=1, at=seq(1,nrow(res),by=xby), labels=FALSE )
  
  if( yax )
  {
    if( is.null(yaxlab) )
      axis( side=2, las=2 )
    else
      axis( side=2, las=2, at=1:ncol(res), labels=yaxlab )
  }
  else
    axis( side=2, labels=NA )
  box()

  symbols( x=resDF$x, y=resDF$y, circles=resDF$radius, add=TRUE, inches=0.07, bg=resDF$bg, fg=resDF$col )

}

plotFitI <- function( rpt, folder="." )
{
  dims <- list( c(4,6), c(3,6) )
  hei <- c(6,5)
  wid <- c(8,8)
  #ylab <- c("Abundance (1000s)","Relative abundance")

  lims <- list( c(18,76), c(15,111) )

  z <- 0
  for( g in 1:rpt$nG )
  {
    pdf( file=paste(folder,"/fig",g+12,"indexFitsg",g,".pdf",sep=""),
         height=hei[g], width=wid[g] )
    par( mfrow=dims[[g]], mar=c(2,2,1,1), oma=c(3,3,0,0) )
    for( t in 1:rpt$nT )
    {
      I_d <- 1e-3*rpt$E_dtg[ ,t,g]/exp(rpt$lnqE_pg[1,g])
      E_d <- 1e-3*rpt$Ihat_dtg[ ,t,g]/exp(rpt$lnqE_pg[1,g])
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

    mtext( side=1, outer=TRUE, line=1.5, text="Ordinal date" )
    mtext( side=2, outer=TRUE, line=1.5, text="Abundance (1000s)" )
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

plotMrFit <- function( rpt )
{
  plot( x=rpt$years, y=rpt$I_t )
  lines( x=rpt$years, y=rpt$mrIhat_t )
}

plotBorderPassage <- function( rpt, folder="." )
{
  pdf( file=paste(folder,"/fig12-borderPassage.pdf",sep=""), height=5, width=7 )

  par( mar=c(3,5,1,1) )

  #t <- !is.na(rpt$I_t)
  t <- rep(TRUE,rpt$nT)
  yr  <- rpt$years[t]
  I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_p[1])
  E_t <- apply(rpt$N_asrpt[ , ,3, , ],4,sum)*1e-3
  sonarN_t <- colSums(rpt$E_dtg[ ,t,2],na.rm=1)*1e-3
  ymax <- max(I_t,E_t,sonarN_t,na.rm=TRUE)

  if( is.finite(rpt$sdrpt[1,5]) )
  {
    Rse <- filter(rpt$sdrpt,par=="borderPass_t")
    Rlow_t <- Rse$lCI*1e-3
    Rupp_t <- Rse$uCI*1e-3
  }
  else
  {
    Rlow_t <- NA*I_t
    Rupp_t <- NA*I_t
  }

  sonarN_t[sonarN_t==0] <- NA

  plot( x=yr, y=I_t, type="n", las=1, yaxs="i", xlab="",
        ylab="Total border passage (1000s)", ylim=c(0,1.1*ymax) )
  grid()
  box()
  
  segments( x0=yr,
            y0=Rlow_t, y1=Rupp_t,
                col="grey80", lwd=2.5 )

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

plotTotalRunSize <- function( rpt, folder="." )
{
  pdf( file=paste(folder,"/totalRunSize.pdf",sep=""), height=5, width=7 )

  par( mar=c(3,5,1,1) )

  t <- !is.na(rpt$I_t)
  t <- rep(TRUE,rpt$nT)
  yr  <- rpt$years[t]
  I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_p)
  E_t <- apply(rpt$N_asrpt[,,3,,],4,sum)*1e-3
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
  runSize_pt <- rpt$runSize_pt*1e-3

  if( is.finite(rpt$sdrpt[1,5]) )
  {
    Rse <- filter(rpt$sdrpt,par=="runSize_pt")
    Rlow_pt <- matrix( data=Rse$lCI*1e-3, nrow=rpt$nP, ncol=rpt$nT )
    Rupp_pt <- matrix( data=Rse$uCI*1e-3, nrow=rpt$nP, ncol=rpt$nT )

    Rse_t <- filter(rpt$sdrpt,par=="runSize_t")
    Rlow_t <- Rse_t$lCI*1e-3
    Rupp_t <- Rse_t$uCI*1e-3
  }
  else
  {
    Rlow_pt <- matrix( data=NA, nrow=rpt$nP, ncol=rpt$nT )
    Rupp_pt <- matrix( data=NA, nrow=rpt$nP, ncol=rpt$nT )
    Rlow_t <- rep(NA,rpt$nT)
    Rupp_t <- rep(NA,rpt$nT)
  }
  

  pdf( file=paste(folder,"/fig15-runSize.pdf",sep=""), height=6.5, width=8 )
  par( mfrow=c(3,3), mar=c(0.5,3,0,0), oma=c(4,2,1,1) )

  for( p in 1:rpt$nP )
  {
    plot( x=range(x), y=c(0,1.1*max(runSize_pt[p,],Rupp_pt[p,],na.rm=TRUE)),
          type="n", axes=FALSE, ylab="" )
    if( p < 7 )
      axis( side=1, labels=NA )
    else
      axis( side=1 )
    axis( side=2, las=1 )
    grid()
    box()
    segments( x0=x, y0=Rlow_pt[p, ], y1=Rupp_pt[p, ], col="grey70", lwd=4 )
    points( x=x, y=runSize_pt[p, ], pch=16 )
    par(font=2)
    legend(x="topright",legend=rpt$stocks[p],bty="n")
    par(font=1)
  }

  runSize_t <- colSums(runSize_pt)

  plot( x=range(x), y=c(0,1.1*max(runSize_t,Rupp_t,na.rm=1)),
        type="n", las=1, yaxs="i" )
  grid()
  box()
  segments( x0=x, y0=Rlow_t, y1=Rupp_t, col="grey70", lwd=4 )
  points( x=x, y=runSize_t, pch=16 )
  par(font=2)
  legend(x="topright",legend="Total",bty="n")
  par(font=1)

  mtext( side=1, text="Year", outer=TRUE, line=2.5 )
  mtext( side=2, text="Run size ('000s)", outer=TRUE, line=0.5 )

  dev.off()

}

plotRunSizeMulti <- function( rptFiles=c("fits/mod1/rpt.Rdata",
                                         "fits/mod2/rpt.Rdata",
                                         "fits/mod3/rpt.Rdata")
                            )
{
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

  Rlow_ist <- array( data=NA, dim=c(3,rpt$nP+1,rpt$nT) )
  Rmle_ist <- array( data=NA, dim=c(3,rpt$nP+1,rpt$nT) )
  Rupp_ist <- array( data=NA, dim=c(3,rpt$nP+1,rpt$nT) )
  for( i in 1:3 )
  {
    load(rptFiles[i])
    Rmle_ist[i,1:rpt$nP, ] <- rpt$runSize_pt*1e-3
    Rse <- filter(rpt$sdrpt,par=="runSize_pt")
    Rlow_ist[i,1:rpt$nP, ] <- matrix( data=Rse$lCI*1e-3, nrow=rpt$nP, ncol=rpt$nT )
    Rupp_ist[i,1:rpt$nP, ] <- matrix( data=Rse$uCI*1e-3, nrow=rpt$nP, ncol=rpt$nT )
  
    Rmle_ist[i,rpt$nP+1, ] <- colSums(rpt$runSize_pt)*1e-3
    Rse <- filter(rpt$sdrpt,par=="runSize_t")
    Rlow_ist[i,rpt$nP+1, ] <- Rse$lCI*1e-3
    Rupp_ist[i,rpt$nP+1, ] <- Rse$uCI*1e-3

  }
  jtr <- c(-0.2,0,0.2)

  stks <- c(rpt$stocks,"Total")

  pdf( file="runSizeMult.pdf", height=8, width=7 )
  par( mfrow=c(rpt$nP+1,1), mar=c(0,2,0,1), oma=c(2,3,2,2) )

  for( s in 1:(rpt$nP+1) )
  {
    ymax <- 1.2*max(Rupp_ist[ ,s,],na.rm=TRUE)
    plot( x=range(x), y=c(0,ymax),
          type="n", axes=FALSE, yaxs="i" )
    grid()
    box()
    legend( x="topright", bty="n", legend=stks[s], cex=1.5 )
    if( s==1 )
    {
      legend( x="topleft", bty="n", col=cols, pch=NA, lwd=2.5,
              legend=c("S mod","E mod","EM mod") )
      legend( x="topleft", bty="n", col="black", pch=15:17, lwd=0,
              legend=c("S mod","E mod","EM mod") )
      axis( side=3 )
    }
      
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

  axis( side=1 )

  mtext( side=2, text="Run size (1000s)", outer=TRUE, line=1 )

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
  par( mfrow=c(ceiling(rpt$nP/2),2), mar=c(2,4,1,1), oma=c(3,2,0,0) )

  for( s in 1:rpt$nP )
  {
    rho_dt <- rpt$rho_dpt[,s,]
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

  cols <- brewer.pal(rpt$nP,"Set1")
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

  for( s in 1:rpt$nP )
  {
    x1_d <- numeric( rpt$nD )
    x2_d <- numeric( rpt$nD )
    for( d in 1:rpt$nD )
    {
      x1_d[d] <- mean(rpt$rho_dpt[d,s, ])
      x2_d[d] <- median(rpt$rho_dpt[d,s, ])
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

plotbg <- function(col=rgb(235,235,235,maxColorValue=255),lty=1,gridcol="white")
{
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=col)
  grid( col=gridcol, lty=lty )
  box()
}

plotHarvestData <- function()
{
  harv <- read.csv(file="data/harvest/harv-processed.csv")
  harv$catch <- harv$catch*1e-3
  yr <- sort(unique(harv$Year))
  nT <- length(yr)
  cols <- brewer.pal(4,"Set1")
  maxy <- max(harv$catch)

  par( mfrow=c(2,2), mar=c(0,3,0,0), oma=c(4.5,4,1,1) )

  harvt <- list()
  harvt[[1]] <- filter( harv, Fishery=="US_Comm" )
  harvt[[2]] <- filter( harv, Fishery=="US_Subsistence" )
  harvt[[3]] <- filter( harv, Fishery=="CDN_Comm" )
  harvt[[4]] <- filter( harv, Fishery=="CDN_Subsistence" )

  labs <- c("US Commercial","US Subsistence","CA Commercial","CA Subsistence")

  for( i in 1:4 )
  {
    maxy <- max(harvt[[i]]$catch)
    plot( x=range(yr), y=c(0,maxy), axes=FALSE, type="n" )
    axis( side=2, las=1 )
    grid()
    box()

    for( a in 4:7 )
    {
      z <- filter(harvt[[i]],variable==a)
      lines( z$Year, z$catch, lwd=2, col=cols[a-3] )
    }

    par(font=2)
    legend( x="topright", bty="n", legend=labs[i], cex=1.2 )
    par(font=1)
    legend( x="topright", bty="n", legend=c("","",paste("Age",4:7)),
            col=c(NA,NA,cols), lwd=c(0,0,2,2,2,2) )

    if( i>2 )
      axis( side=1 )

  }

  mtext( side=1, text="Year", outer=TRUE, line=2.5, cex=1.2 )
  mtext( side=2, text="Catch (1000s)", outer=TRUE, line=1, cex=1.2 )

}

plotHarvestDatav1 <- function()
{
  harv <- read.csv(file="data/harvest/harv-processed.csv")
  yr <- sort(unique(harv$Year))
  nT <- length(yr)

  maxy <- max(harv$catch)

  par( mfrow=c(5,ceiling(nT/5)), mar=c(0,0,0,0), oma=c(4.5,5.5,1,1) )

  for( t in 1:nT )
  {
    harvt <- filter( harv, Year==yr[t] )
    harvUC <- filter( harvt, Fishery=="US_Comm" )$catch
    harvUS <- filter( harvt, Fishery=="US_Subsistence" )$catch
    harvCC <- filter( harvt, Fishery=="CDN_Comm" )$catch
    harvCS <- filter( harvt, Fishery=="CDN_Subsistence" )$catch

    plot( x=c(3.5,7.5), y=c(0,maxy), axes=FALSE, type="n" )
    grid()
    box()
    if(length(harvUC))
      rect( xleft=4:7-0.3, xright=4:7-0.15, ybottom=0, ytop=harvUC, border="black", col=NA )
    if(length(harvUS))
      rect( xleft=4:7-0.15, xright=4:7, ybottom=0, ytop=harvUS, border="grey70", col=NA )
    if(length(harvCC))
      rect( xleft=4:7, xright=4:7+0.15, ybottom=0, ytop=harvCC, border=NA, col="black" )
    if(length(harvCS))
      rect( xleft=4:7+0.15, xright=4:7+0.3, ybottom=0, ytop=harvCS, border=NA, col="grey70" )

    if( t>(nT-ceiling(nT/5)) )
      axis( side=1 )
    if( t %% ceiling(nT/5) == 1 )
      axis( side=2, las=1 )

    par(font=2)
    legend( x="topright", bty="n", legend=yr[t] )
    par(font=1)

  }

  mtext( side=1, text="Age (yr)", outer=TRUE, line=2.5, cex=1.2 )
  mtext( side=2, text="Catch", outer=TRUE, line=3.5, cex=1.2 )

}

plotCatchFit <- function( rpt, folder="." )
{
  obsC_rht <- rpt$obsC_rht*1e-3
  C_rht <- apply( rpt$C_asrht, 3:5, sum )*1e-3
  x <- rpt$years

  pdf( file=paste0(folder,"/fig3-catchFit.pdf"), height=6.5, width=7 )
  par( mfrow=c(rpt$nR,rpt$nH), mar=c(0.5,3,0,1), oma=c(4,3,1,1) )
  for( r in 1:rpt$nR )
  {
    for( h in 1:rpt$nH )
    {
      plot( x, C_rht[r,h, ], ylim=c(0,1.1*max(obsC_rht[r,h, ],C_rht[r,h, ])), axes=FALSE, type="n" )
      grid()
      box()
      rect( xleft=x-0.3, xright=x+0.3, ybottom=0, ytop=obsC_rht[r,h, ], col="grey80", border=NA )
      points( x, C_rht[r,h, ] )

      if( r < rpt$nR )
        axis( side=1, labels=NA )
      else
        axis( side=1 )

      axis( side=2, las=1 )

      par(font=2)
      legend(x="topright",legend=paste(rpt$regions[r],rpt$fisheryType[h]),bty="n")
      par(font=1)
      legend( x="topright", legend=c("","Obs.","Est."), pch=c(NA,15,1),
              col=c("black","grey80","black"), pt.cex=c(1,1.3,1), bty="n" )
    }
  }
  mtext( side=1, text="Year", outer=TRUE, line=2.5 )
  mtext( side=2, text="Landings ('000s)", outer=TRUE, line=1.5 )
  dev.off()
}

plotCatchAgeFit <- function( rpt, folder="." )
{
  f <- 6
  x <- 1:rpt$nA
  for( r in 1:rpt$nR )
  {
    for( h in 1:rpt$nH )
    {
      tseq <- which(!is.na(rpt$xC_arht[1,r,h, ]))
      n <- length(tseq)
      if( n < rpt$nT )
        ncl <- 9
      else
        ncl <- 10
      pdf( file=paste0(folder,"/fig",f,"-catchAgeFit-r",r,"-h",h,".pdf"), height=4, width=7 )
      par( mfrow=c(ceiling(n/ncl),ncl), mar=c(0,0,0,0), oma=c(4,4,2,1) )
      obs_at <- rpt$xC_arht[ ,r,h, ]
      est_at <- rpt$xhatC_arht[ ,r,h, ]
      for( t in 1:n )
      {
        plot( x=c(0.5,max(x)+0.5), y=c(0,1.05), type="n", axes=FALSE )
        grid()
        box()
        rect( xleft=x-0.3, xright=x+0.3, ybottom=0, ytop=obs_at[ ,tseq[t]]/sum(obs_at[ ,tseq[t]]), border=NA, col="grey80" )
        points( x=x, y=est_at[ ,tseq[t]] )
        lines( x=x, y=est_at[ ,tseq[t]] )
        if( t>(n-ncl) )
        {
          axis( side=1, at=1:rpt$nA, labels=NA )
          axis( side=1, at=1:rpt$nA, labels=rpt$ages, cex.axis=0.9, line=-0.2, lty=0 )
        }
  
        if( t %% ncl == 1 )
          axis( side=2, las=1 )
        box()
        legend(x="topright",legend=rpt$years[t],bty="n")
      }
      plot( x=0, y=0, type="n", axes=FALSE, xlab="", ylab="" )
      legend( x="bottomright", bty="n", legend=c("Obs.","Est."),
              col=c("grey80","black"), pch=c(15,1) )
      mtext( side=1, outer=TRUE, line=2.5, text="Age class" )
      mtext( side=2, outer=TRUE, line=2.5, text="Proportion" )  
      mtext( side=3, outer=TRUE, text=paste(rpt$regions[r],rpt$fisheryType[h]) )
      dev.off()
      f <- f+1
    }
  } 
}

plotCatchAgeFitBubble <- function( rpt, folder="." )
{
  x <- 1:(rpt$nA*2)
  yaxlab <- rpt$ages
  pdf( file=paste0(folder,"/catchAgeFitBubble.pdf"), height=7, width=7 )
  par( mfrow=c(rpt$nR,rpt$nH), mar=c(0.5,0.5,0,0), oma=c(5,4,1,1) )
  for( r in 1:rpt$nR )
  {
    for( h in 1:rpt$nH )
    {
      obs_at <- apply( rpt$xC_arht[ ,r,h, ], 2, standardize )
      est_at <- rpt$xhatC_arht[ ,r,h, ]      
      res_at <- log(est_at) - log(obs_at)
      bubblePlot2(t(res_at),xax=r==rpt$nR,yax=h==1,yaxlab=yaxlab)
    }
  }
  mtext( side=1, text="Year", outer=TRUE, line=2.5 )
  mtext( side=2, text="Age class", outer=TRUE, line=2 )
  dev.off() 
}

plotBorderAgeFit <- function( rpt, folder="." )
{
  x <- 1:(rpt$nA*2)
 
  labs <- c(paste0("M",rpt$ages),paste0("F",rpt$ages))

  for( g in 1:rpt$nG )
  {
    tseq <- which(!is.na(rpt$xB_jgt[1,g, ]))
    n <- length(tseq)
    pdf( file=paste0(folder,"/fig",g+3,"-borderAgeFit",g,".pdf"), height=6, width=7 )
    par( mfrow=c(ceiling(n/4),4), mar=c(0,0,0,0), oma=c(4,4,2,1) )
    obs_jt <- rpt$xB_jgt[ ,g,tseq]
    est_jt <- rbind(rpt$xB_asgt[ ,1,g,tseq],rpt$xB_asgt[ ,2,g,tseq])
    for( t in 1:n )
    {
      plot( x=c(0.5,max(x)+0.5), y=c(0,0.85), type="n", axes=FALSE )
      grid()
      box()
      rect( xleft=x-0.3, xright=x+0.3, ybottom=0, ytop=obs_jt[ ,t]/sum(obs_jt[ ,t]), border=NA, col="grey80" )
      points( x=x, y=est_jt[ ,t] )
      lines( x=x, y=est_jt[ ,t] )
      if( t>(n-4) )
      {
        axis( side=1, at=1:(2*rpt$nA), labels=NA)
        axis( side=1, at=1:(2*rpt$nA), labels=labs, cex.axis=0.65, line=-0.25, lty=0 )
      }

      if( t %% 4 == 1 )
        axis( side=2, las=1 )
      box()
      legend(x="topright",legend=rpt$years[tseq[t]],bty="n")
    }
    plot( x=0, y=0, type="n", axes=FALSE, xlab="", ylab="" )
    legend( x="bottomright", bty="n", legend=c("Obs.","Est."),
              col=c("grey80","black"), pch=c(15,1) )
    mtext( side=3, outer=TRUE, text=rpt$gears[g], line=0.2 )
    mtext( side=1, outer=TRUE, line=2.5, text="Age/sex class" )
    mtext( side=2, outer=TRUE, line=2.5, text="Proportion" )  
  dev.off() 
  }
}


plotBorderAgeFitBubble <- function( rpt, folder="." )
{
  x <- 1:(rpt$nA*2)
  yaxlab <- c(paste0("M",rpt$ages),paste0("F",rpt$ages))
  labs <- c(paste0("M",1:rpt$nA),paste0("F",1:rpt$nA))
  pdf( file=paste0(folder,"/borderAgeFitBubble.pdf"), height=5, width=6 )
  par( mfrow=c(2,1), mar=c(0.5,2,0,0), oma=c(3,3,1,1) )
  for( g in 1:rpt$nG )
  {
    obs_jt <- rpt$xB_jgt[ ,g, ]
    obs_jt <- apply( obs_jt, 2, standardize )
    est_jt <- rbind(rpt$xB_asgt[ ,1,g, ],rpt$xB_asgt[ ,2,g, ])
    res_jt <- log(est_jt) - log(obs_jt)
    bubblePlot2(t(res_jt),xax=g==rpt$nG,yaxlab=yaxlab)
  }
  mtext( side=1, text="Year", outer=TRUE, line=2 )
  mtext( side=2, text="Age/sex class", outer=TRUE, line=2 )
  dev.off() 
}

plotRecruitment <- function( rpt, folder="." )
{
  pdf( file=paste0(folder,"/recruitment.pdf"), height=7, width=6 )
  par( mfrow=c(rpt$nP/2,2), mar=c(0,0,0,0), oma=c(4,6,2,1) )
  for( p in 1:rpt$nP )
  {
    exp_y <- rpt$R_py[p, ]
    rea_y <- exp(rpt$lnR_py[p, ])

    plot( x=c(1,rpt$nY), y=c(0,max(exp_y,rea_y)), type="n", axes=FALSE )
    grid()
    box()
    points( exp_y )
    lines( rea_y )
    if( p>(rpt$nP-2) )
      axis( side=1 )

    if( p %% 2 == 1 )
      axis( side=2, las=1 )
    box()

    legend( x="topright", legend=c("Expected","Realized"), bty="n",
            pch=c(1,NA), lty=c(0,1) )

  }
  mtext( side=1, outer=TRUE, line=2.5, text="Year" )
  mtext( side=2, outer=TRUE, line=4.5, text="Recruitment" )  
  dev.off() 
}


plotSelectivity <- function( rpt, folder="." )
{
  v_astm <- rpt$v_astm
  cols <- brewer.pal(rpt$nM,"Set1")
  pdf( file=paste0(folder,"/selectivity.pdf"), height=7, width=6 )
  par( mfrow=c(ceiling(rpt$nT/5),5), mar=c(0,0,0,0), oma=c(4,6,2,1) )
  for( t in 1:rpt$nT )
  {
    plot( x=range(rpt$ages)+c(-0.25,0.25), y=c(0,1.08), type="n", axes=FALSE )
    grid(ny=rpt$nA)
    box()
    for( m in 1:rpt$nM )
    {
      lines( x=rpt$ages, y=v_astm[ ,1,t,m], col=cols[m] )
      lines( x=rpt$ages, y=v_astm[ ,2,t,m], col=cols[m], lty=2 )
    }
    legend(x="bottomleft",legend=rpt$years[t],bty="n")

    if( t>(rpt$nT-5) )
      axis(side=1,at=rpt$ages)
    if( t %% 5 == 1 )
      axis( side=2, las=1 )
  }
  plot( x=0,y=0,axes=FALSE,type="n" )
  labs <- c(paste(rpt$meshSizes, "in, Male"))
  legend( x="bottomright", legend=labs, bty="n", col=c(cols),
          lty=c(rep(1,3)) )
  plot( x=0,y=0,axes=FALSE,type="n" )
  labs <- c(paste(rpt$meshSizes, "in, Female"))
  legend( x="bottomright", legend=labs, bty="n", col=c(cols),
          lty=c(rep(2,3)) )

  mtext( side=1, outer=TRUE, line=2.5, text="Age (yr)" )
  mtext( side=2, outer=TRUE, line=4.5, text="Selectivity" )  
  dev.off() 
}

plotSelectivityAvg <- function( rpt, folder="." )
{
  v_astm <- rpt$v_astm

  vAll_asm <- apply( rpt$v_astm, c(1,2,4), mean )
  vf10_asm <- apply( rpt$v_astm[ , ,1:10, ], c(1,2,4), mean )
  vl10_asm <- apply( rpt$v_astm[ , ,(rpt$nT-9):rpt$nT, ], c(1,2,4), mean )

  cols <- brewer.pal(3,"Set1")
  cols <- c(cols[1],"black",cols[2])
  pdf( file=paste0(folder,"/fig17-selectivityAvg.pdf"), height=5.5, width=3.5 )
  par( mfrow=c(3,1), mar=c(0.5,0,0,0), oma=c(4,6,2,1) )
  
  for( m in 1:rpt$nM )
  {
    plot( x=range(rpt$ages), y=c(0,1.08), type="n", axes=FALSE )
    axis( side=1, labels=NA, at=rpt$ages )
    axis( side=2, las=1 )
    plotbg()
    points( x=rpt$ages, y=vf10_asm[ ,1,m]/max(vf10_asm[ ,1,m]), col=cols[1], pch=16, cex=1.1 )
    lines( x=rpt$ages, y=vf10_asm[ ,1,m]/max(vf10_asm[ ,1,m]), col=cols[1] )
    points( x=rpt$ages, y=vf10_asm[ ,2,m]/max(vf10_asm[ ,2,m]), col=cols[1], pch=18, cex=1.1 )
    lines( x=rpt$ages, y=vf10_asm[ ,2,m]/max(vf10_asm[ ,2,m]), col=cols[1], lty=2 )

    points( x=rpt$ages, y=vAll_asm[ ,1,m]/max(vAll_asm[ ,1,m]), col=cols[2], pch=16, cex=1.1 )
    lines( x=rpt$ages, y=vAll_asm[ ,1,m]/max(vAll_asm[ ,1,m]), col=cols[2] )
    points( x=rpt$ages, y=vAll_asm[ ,2,m]/max(vAll_asm[ ,2,m]), col=cols[2], pch=18, cex=1.1 )
    lines( x=rpt$ages, y=vAll_asm[ ,2,m]/max(vAll_asm[ ,2,m]), col=cols[2], lty=2 )

    points( x=rpt$ages, y=vl10_asm[ ,1,m]/max(vl10_asm[ ,1,m]), col=cols[3], pch=16, cex=1.1 )
    lines( x=rpt$ages, y=vl10_asm[ ,1,m]/max(vl10_asm[ ,1,m]), col=cols[3] )
    points( x=rpt$ages, y=vl10_asm[ ,2,m]/max(vl10_asm[ ,2,m]), col=cols[3], pch=18, cex=1.1 )
    lines( x=rpt$ages, y=vl10_asm[ ,2,m]/max(vl10_asm[ ,2,m]), col=cols[3], lty=2 )

    par(font=2)
    legend(x="bottomleft",legend=paste(rpt$meshSizes[m],"in. mesh"),bty="n")
    par(font=1)

    legend( x=5.2, y=0.3, bty="n", pt.cex=1.1,
            pch=c(16,18), col=c("black","black"), lty=c(1,2),
            legend=c("Male","Female") )
    legend( x=6.1, y=0.4, bty="n",
            pch=c(15,15,15), col=cols,
            legend=c("First 10 yrs","All yrs", "Last 10 yrs") )

  }
  axis( side=1, at=rpt$ages )

  mtext( side=1, outer=TRUE, line=2.5, text="Age (yr)" )
  mtext( side=2, outer=TRUE, line=4.5, text="Selectivity" )  
  dev.off() 
}


plotStockRec <- function( rpt, folder=".", sType=1 )
{
  S_pt <- rpt$Z_pt
  R_py <- exp(rpt$lnR_py)

  x <- as.numeric(substr(rpt$years,3,3))
  x <- match(x,unique(x))

  cols <- c( rgb(225,31,39,maxColorValue=255),
             rgb(251,173,104,maxColorValue=255),
             rgb(255,254,76,maxColorValue=255),
             rgb(172,217,232,maxColorValue=255),
             rgb(71,118,178,maxColorValue=255) )
  
  pdf( file=paste0(folder,"/fig20-stockRec.pdf"), height=7, width=6 )
  par( mfrow=c(ceiling(rpt$nP/2),2), mar=c(3,4,0,0), oma=c(2,2,1,1) )
  for( p in 1:rpt$nP )
  {
    S <- rep(NA,rpt$nY)
    R <- rep(NA,rpt$nY)
    for( y in 8:rpt$nY )
    {
      R[y] <- R_py[p,y]
      S[y] <- S_pt[p,y-7]
    }

    plot( x=c(0,1.1*max(S*1e-3,na.rm=1)), y=c(0,1.15*max(R*1e-3,na.rm=1)),
          las=1, xlab="", ylab="", type="n" )
    plotbg()
    points( x=S*1e-3, y=R*1e-3, col=cols[x], pch=16 )
    points( x=S*1e-3, y=R*1e-3 )
    S <- seq(from=0,to=1.2*max(S,na.rm=1),length.out=1e4)
    R <- rpt$alpha_p[p]*S*exp(-rpt$beta_p[p]*S)
    lines(S*1e-3,R*1e-3)
    par(font=2)
    legend(x="topright",legend=rpt$stocks[p],bty="n")
    par(font=1)
    if(p==1)
    {
      legend( x="topright", legend=c("","1980s","1990s","2000s","2010s","2020s"),
          col=c(NA,cols), pch=16, cex=0.95, bty="n", lty=0 ) 
      legend( x="topright", legend=c("","1980s","1990s","2000s","2010s","2020s"),
          bty="n", pch=c(NA,rep(1,5)), cex=0.95 ) 
    }
  }

  xtext <- "Spawner abundance ('000s)"
  if( sType==2 )
    xtext <- "Total egg production (millions)"
  else if( sType==3 )
    xtext <- "Total egg mass (mt)"

  mtext( side=1, outer=TRUE, line=0.5, text=xtext )
  mtext( side=2, outer=TRUE, line=0.5, text="Recruits ('000s)" )  
  dev.off() 
}


plotpFem <- function( rpt, folder="." )
{
  pdf( file=paste0(folder,"/fig18-pFem.pdf"), height=4.5, width=5 )
  x <- rpt$broodYears[1:(length(rpt$broodYears)-2)]
  y <- rpt$pFem_y[1:(length(rpt$broodYears)-2)]

  if( is.finite(rpt$sdrpt[1,5]) )
  {
    Rse <- filter(rpt$sdrpt,par=="pFem_y")
    Rlow_t <- Rse$lCI[1:(length(rpt$broodYears)-2)]
    Rupp_t <- Rse$uCI[1:(length(rpt$broodYears)-2)]
  }
  else
  {
    Rlow_t <- rep(NA,rpt$nT)
    Rupp_t <- rep(NA,rpt$nT)
  }

  plot( x=x, y=y, type="n", las=1, xlab="Brood Year", ylab="Proportion female", ylim=c(0,1) )
  grid()
  box()
  segments( x0=x, y0=Rlow_t, y1=Rupp_t, lwd=6, col="grey80"  )
  lines( x=x, y=y )
  points( x=x, y=y, pch=16, cex=0.85 )
  dev.off()
}


plotRelativeReproOutput <- function()
{
  folders <- c("mod2","mod3")

  ylabs <- c("Eggs per spawner ('000s)",
             "Egg mass per spawner (kg)")

  pdf( file="figs/fig21-relReproOutput.pdf", height=6, width=5 )
  par( mfrow=c(2,1), mar=c(0.5,4,0,0), oma=c(4,0,1,1) )

  for( i in 1:2 )
  {
    load(paste0("fits/",folders[i],"/rpt.Rdata"))
    x <- rpt$years
    y <- rpt$relReproOutput_t
  
    if( is.finite(rpt$sdrpt[1,5]) )
    {
      Rse <- filter(rpt$sdrpt,par=="relReproOutput_t")
      Rlow_t <- Rse$lCI
      Rupp_t <- Rse$uCI
    }
    else
    {
      Rlow_t <- rep(NA,rpt$nT)
      Rupp_t <- rep(NA,rpt$nT)
    }

    plot( x=range(x), y=range(y,Rlow_t,Rupp_t,na.rm=TRUE),
          type="n", axes=FALSE, xlab="Year",
          ylab=ylabs[i] )
    axis( side=1, labels=NA )
    axis( side=2, las=1 )
    grid()
    box()
    segments( x0=x, y0=Rlow_t, y1=Rupp_t, lwd=6, col="grey80"  )
    lines( x=x, y=y )
    points( x=x, y=y, pch=16, cex=0.8 )

  }
  axis( side=1 )

  mtext(side=1,text="             Year",outer=TRUE,line=2.5)

  dev.off()
}


plotMETF <- function( rpt, folder="." )
{
  S_pt <- rpt$Z_pt
  R_py <- exp(rpt$lnR_py)

  cols <- brewer.pal(rpt$nM,"Set1")
  pdf( file=paste0(folder,"/stockRec.pdf"), height=7, width=6 )
  par( mfrow=c(ceiling(rpt$nP/2),2), mar=c(3,4,0,0), oma=c(2,2,1,1) )
  for( p in 1:rpt$nP )
  {
    S <- rep(NA,rpt$nY)
    R <- rep(NA,rpt$nY)
    for( y in 8:rpt$nY )
    {
      R[y] <- R_py[p,y]
      S[y] <- S_pt[p,y-7]
    }

    plot( x=S*1e-3, y=R*1e-3, las=1 )
    grid()
    box()
    S <- seq(from=0,to=1.2*max(S,na.rm=1),length.out=1e4)
    R <- rpt$alpha_p[p]*S*exp(-rpt$beta_p[p]*S)
    lines(S*1e-3,R*1e-3)
    legend(x="topright",legend=rpt$stocks[p],bty="n")
  }

  mtext( side=1, outer=TRUE, line=0.5, text="Spawners ('000s)" )
  mtext( side=2, outer=TRUE, line=0.5, text="Recruits ('000s)" )  
  dev.off() 
}

plotSelectivityByLength <- function( rpt, folder="." )
{
  cols <- brewer.pal(3,"Set1")[1:2]
  load("data/chinookYkData.Rdata")
  L_as <- apply(chinookYkData$rawMETF_ast,1:2,mean,na.rm=TRUE)
  pars <- list(lambda=rpt$pearLambda,theta=rpt$pearTheta,sigma=rpt$pearSigma,tau=rpt$pearTau)
  
  lengthSeq <- 300:1200
  mm_per_inch <- 25.4 
  
  # Pearson function
  Pearson = function(fish_length, mesh_stretch_length, lambda, theta, sigma, tau, standardize = FALSE) {
    
    # enforce parameter constraints
    if (any(sigma <= 0)) stop ("sigma must be > 0")
    if (any(theta <= 0)) stop ("theta must be > 0")
    
    # perimeter = length of one side * 4
    mesh_perimeter = mesh_stretch_length * 2
    
    # calculate RLM: ratio of length to mesh perimeter
    # units don't matter so long as they are the same
    rlm = fish_length/mesh_perimeter
    
    # separate calculation into 4 "terms"
    # "t2" is used in two places exactly, so this makes it a bit easier and avoids potential typos
    t1 = (1 + lambda^2/(4 * theta^2))^theta
    t2 = rlm - (sigma * lambda)/(2 * theta) - tau
    t3 = (1 + t2^2/sigma^2)^-theta
    t4 = exp(-lambda * (atan(t2/sigma) + atan(lambda/(2 * theta))))
    
    # calculate the entire Pearson function
    v = t1 * t3 * t4
    
    # standardize the output so only one age/sex is fully vuln for a gear
    if (standardize) out = v/max(v) else out = v
    
    # return the output
    return(out)
  }
  
  # calculate the selectivity for each length in the sequence for two mesh sizes
  #v_8 = do.call(Pearson, append(pars, list(fish_length = lengthSeq, mesh_stretch_length = 8 * mm_per_inch)))
  #v_6 = do.call(Pearson, append(pars, list(fish_length = lengthSeq, mesh_stretch_length = 6 * mm_per_inch)))
  
  # plot the curves
  par( mar=c(4,4,1,1), oma=c(0,0,0,0) )
  pdf( file=paste0(folder,"/fig16-selectivityByLength.pdf"), height=6.5, width=6.5 )
  plot( x=range(lengthSeq), y=c(0,1), las=1, type="n", xlab="Fish Length (METF; mm)", ylab="Selectivity")
  plotbg()
  for( s in 1:2 )
    for( a in 1:rpt$nA )
        abline( v=L_as[a,s], col=cols[s], lty=a )
  maxS <- numeric(rpt$nM)
  for( m in 1:rpt$nM )
  {
    y <- do.call(Pearson, append(pars, list(fish_length = lengthSeq, mesh_stretch_length = rpt$meshSizes[m] * mm_per_inch)))
    lines( x=lengthSeq, y=y, lwd=1.5, lty=m )
    maxS[m] <- lengthSeq[which.max(y)]
  }

  leg <- c(paste0("Male age-",rpt$ages),paste0("Female age-",rpt$ages))
  legend("topleft", legend=paste(rpt$meshSizes,"in"), lty=1:rpt$nM, bty="n", lwd=1.5)
  legend("topright", legend=leg, lty=c(1:rpt$nA,1:rpt$nA), bty="n", col=c(rep(cols[1],4),rep(cols[2],4)))
  dev.off()
  maxS
#  plot(v_8 ~ lengthSeq, type = "l", ylim = c(0,1), xlab = "Fish Length (mm)", ylab = "Relative Selectivity")
#  lines(v_6 ~ lengthSeq, lty = 2)
#  legend("topleft", legend = c("8", "6"), title = "Mesh (in)", lty = c(1, 2), bty = "n")
  
#  # get the Kusko mean length at age/sex data, calculate average post-2010, and plot female MLA
#  Kusko_MLA = read.csv("https://raw.githubusercontent.com/bstaton1/esc-qual-ms-analysis/master/2-model-fit/inputs/esc-mean-length.csv")
#  means = colMeans(Kusko_MLA[Kusko_MLA$year > 2010,-1])
#  means = means[1:4] # only females
#  
#  # draw the mean length-at-age
#  abline(v = means)
#  points(x = means, y = rep(0.05, length(means)), pch = 15, col = "white", cex = 2)
#  text(x = means, y = rep(0.05, length(means)), labels = toupper(names(means)))
  
  # when applied to the mean lengths of discrete ages in an estimation model
  # it is helpful to standardize such that one age/sex class has selectivity = 1
  # this allows F to be identifiable and interpretted as mortality of fully selected fish
  # we estimated time-constant selectivity parameters, but used time-varying mean length data
  # and time- and gear-varying F parameters
  # here is an example
  
#  # get the means again
#  means = colMeans(Kusko_MLA[Kusko_MLA$year > 2010,-1])
#  
#  # apply the selectivity function to only the mean lengths-at-age/sex
#  # this time use standardize = TRUE
#  v_8 = do.call(Pearson, append(use_params, list(fish_length = means, mesh_stretch_length = 8 * mm_per_inch, standardize = TRUE)))
#  v_6 = do.call(Pearson, append(use_params, list(fish_length = means, mesh_stretch_length = 6 * mm_per_inch, standardize = TRUE)))
#  
#  # select the mesh size you'd like to see
#  v_use = v_6
}


plotProbReturnAtAge <- function( rpt, folder="." )
{
  eta_asy <- rpt$eta_aspy[ , ,1, ]
  
  cols <- brewer.pal( rpt$nA, "Set1" )

  x <- rpt$broodYears

  bump_sa <- matrix( data=0, nrow=2, ncol=rpt$nA )
  bump_sa[1,3] <- 0.02
  bump_sa[2,2] <- 0.08
  bump_sa[2,rpt$nA] <- 0.04

  pdf( file=paste0(folder,"/fig19-probReturnAtAge.pdf"), height=6.5, width=7 )
  par( mfrow=c(2,1), mar=c(0.5,3,0,1), oma=c(4,2,1,1) )
  for( s in 1:2 )
  {
    plot( x=range(x)+c(0,2.6), y=c(0,1), axes=FALSE, type="n" )
    grid()
    box()
    for( a in 1:rpt$nA ) 
      lines( x, eta_asy[a,s, ], col=cols[a], lwd=1.5 )
    axis( side=1, labels=NA )
    axis( side=2, las=1 )

    text( x=x[rpt$nY]+2.1, y=eta_asy[ ,s,rpt$nY]+bump_sa[s, ],
      col=cols, labels=paste0("Age-",rpt$ages), cex=0.9 )
    par(font=2)
    legend( x="topleft", legend=c("Male","Female")[s], bty="n" )
    par(font=1)
  }
  axis( side=1 )
  mtext( side=1, text="Brood year", outer=TRUE, line=2.5 )
  mtext( side=2, text="Proportion returning-at-age", outer=TRUE, line=0.5 )
  dev.off()
}







