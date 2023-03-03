processASL <- function( fillNAs=TRUE )
{
	ages <- 4:7
	mmPerInch <- 25.4

# Hamazaki 2018
  fishWheelSel_sa <- rbind( c(1,0.295,0.123,0.118),
  													c(0.36,0.129,0.092,0.045) )

	# SAMPLE SIZES
	# mesh     N
	#	2.75    13
	#	4       36
	#	5.25  2426
	#	5.5     32
	#	5.75     2
	#	6.5   2380
	#	7.25     1
	#	7.5   3142
	#	8.5   1620

	df <- read.csv("~/yukon-chinook/data/border-assessment/aslEagle.csv") %>%
				mutate( age=ageFresh+ageSalt+1,
								meshMM=mesh*mmPerInch*2,
								RLM=length/meshMM) %>%
				filter( ASLProjectType=="Test Fishing",
								!is.na(age),
								age>3,
								sexID!=-1,
								mesh%in%c(5.25,6.5,7.5,8.5) ) %>%
				select( sampleYear, mesh, sexID, age, length, RLM )
browser()
	dfLen <- df %>%
					 group_by( sampleYear, sexID, age )

	yr <- min(df$sampleYear):max(df$sampleYear)
	
	meshIn <- sort(unique(df$mesh))
	nMesh <- length(meshIn)

	# get means across years
	dfmean <- df %>%
						group_by(age,sexID) %>%
						summarize( length=mean(length,na.rm=TRUE) )

	RLM_astm <- array( data=NA, dim=c(length(ages),2,length(yr),nMesh),
					  				 dimnames=list(NULL,c("Male","Female"),yr,meshIn) )
	L_ast  <- array( data=NA, dim=c(length(ages),2,length(yr)),
					  				 dimnames=list(NULL,c("Male","Female"),yr) )
	nL_ast  <- array( data=0, dim=c(length(ages),2,length(yr)),
					  				  dimnames=list(NULL,c("Male","Female"),yr) )	
	Lmean_as  <- array( data=NA, dim=c(length(ages),2),
					  				 dimnames=list(NULL,c("Male","Female")) )	

	par( mfrow=c(2,2), mar=c(0.5,0.5,0,0), oma=c(4,5,1,1) )
	for( m in 1:nMesh )
	{
		for( i in 1:length(ages) )
		{
			dfi <- filter( df, age==ages[i] ) %>%
						 group_by(sampleYear,sexID,mesh) %>%
						 summarize(length=mean(length),RLM=mean(RLM))
			#dfm <- filter( dfi, sexID==1 )$length
			#dfmx <- filter( dfi, sexID==1 )$sampleYear
			#dfmr <- filter( dfi, sexID==1 )$RLM
			#dff <- filter( dfi, sexID==2 )$length
			#dffx <- filter( dfi, sexID==2 )$sampleYear
			#dffr <- filter( dfi, sexID==2 )$RLM
		
			for( s in 1:2 )
			{
				dfim <- filter( dfi, sexID==s, mesh==meshIn[m] )
				RLM_astm[i,s,dfim$sampleYear-min(yr)+1,m] <- dfim$RLM

				dfimLen <- filter( dfLen, sexID==s, age==ages[i] ) %>% group_by( sampleYear ) %>% summarize( n=length(length), length=mean(length) )
				L_ast[i,s,dfimLen$sampleYear-min(yr)+1] <- dfimLen$length

				dfimn <- dfimLen
				nL_ast[i,s,dfimn$sampleYear-min(yr)+1] <- dfimn$n



#				RLM_astm[i,s,is.na(RLM_astm[i,s, ,m]),m] <- mean(RLM_astm[i,s, ,m],na.rm=TRUE)


				if( m==1 )
				{
					Lmean_as[i,s] <- filter(dfmean,sexID==s,age==ages[i])$length
				}

			}
		}
	}


	save( Lmean_as, file="data/modelObjects/Lmean_as.Rdata" )
	save( L_ast, file="data/modelObjects/L_ast.Rdata" )
	save( nL_ast, file="data/modelObjects/nL_ast.Rdata" )
	save( RLM_astm, file="data/modelObjects/RLM_astm.Rdata" )

}
# units for mesh and length ?
# mesh presumably inches
# length = mm?
# mesh = "stretch length"


plotRLM <- function()
{
	RLM_astm <- processASL()

	for( i in 1:length(ages) )
	{
		dfi <- filter( df, age==ages[i] ) %>%
					 group_by(sampleYear,sexID,meshIn) %>%
					 summarize(length=mean(length),RLM=mean(RLM))
		dfm <- filter( dfi, sexID==1 )$length
		dfmx <- filter( dfi, sexID==1 )$sampleYear
		dfmr <- filter( dfi, sexID==1 )$RLM
		dff <- filter( dfi, sexID==2 )$length
		dffx <- filter( dfi, sexID==2 )$sampleYear
		dffr <- filter( dfi, sexID==2 )$RLM
	
		#yrng <- range( dfmr, dffr, na.rm=1 )
		yrng <- c(1.5,3.4)

		plot( range( yr ), yrng, type="n", axes=FALSE )
		plotbg()
		abline(h=filter(dfmean,sexID==1,age==ages[i])$RLM,lty=2)
		abline(h=filter(dfmean,sexID==2,age==ages[i])$RLM,lty=2,col="red")
		points( dfmx, dfmr, lwd=2 ,pch=0 )
		points( dffx, dffr, col="red", lwd=2 ,pch=1 )

		RLM_ast[i,1,dfmx-min(yr)+1] <- dfmr
		RLM_ast[i,2,dffx-min(yr)+1] <- dffr

		legx <- "topright"
		if( i > 2 )
		{
			axis( side=1 )
			legx <- "bottomright"
		}
		else
			axis( side=1, labels=FALSE )
		if( i %in% c(1,3) )
			axis( side=2, las=1 )
		else
			axis( side=2, labels=FALSE )

		par(font=2)
		legend(x="topleft",legend=paste("Age",ages[i]),bty="n")
		par(font=1)

		legend(x=legx,bty="n",legend=c("Male","Female"),pch=0:1,
					col=c("black","red"),lwd=2,lty=0)

	}	
		mtext( side=1, text="Year", line=2, outer=TRUE )
	mtext( side=2, text="Mean ratio of fish length to mesh perimeter", line=3.2, outer=TRUE )
}

plotSel <- function()
{
	load("data/chinookYkData.Rdata")
	RLM_astm<-chinookYkData$RLM_astm
	colm <- brewer.pal(9,"YlGn")[8:5]
	colf <- brewer.pal(9,"YlOrRd")[8:5]
	nT <- dim(RLM_astm)[3]
	nM <- dim(RLM_astm)[4]
	yr <- dimnames(RLM_astm)[[3]]
	meshIn <- dimnames(RLM_astm)[[4]]
	par( mfrow=c(4,4), mar=c(0.5,0.5,0,0), oma=c(4,5,1,1) )
	for( t in 1:nT )
	{
		plot( c(4,7), c(0,1), type="n", axes=FALSE )
		plotbg()
		for( m in 1:nM )
		{
			selM_a <- getPearsonSel(RLM_astm[ ,1,t,m])
			selF_a <- getPearsonSel(RLM_astm[ ,2,t,m])
			lines( 4:7, selM_a, lwd=2, col=colm[m], lty=m )
			lines( 4:7, selF_a, lwd=2, col=colf[m], lty=m )
		}

		if( t>10 )
			axis( side=1, labels=4:7, at=4:7 )
		else
			axis( side=1, labels=NA, at=4:7 )

		if( t%in%seq(1,nT,by=4) )
			axis( side=2, las=1 )
		else
			axis( side=2, labels=NA )			

		par(font=2)
		legend(x="topright",legend=yr[t],bty="n",cex=1.2)
		par(font=1)

	}

	maleLeg <- paste0("Male, mesh=",meshIn," in")
	femaleLeg <- paste0("Female, mesh=",meshIn," in")
	plot( x=0, y=0, axes=FALSE, type="n" )
	legend( x="bottomleft", bty="n", col=c(colm,colf), lwd=2, lty=c(1:4,1:4),
					legend=c(maleLeg,femaleLeg) )

	mtext( side=1, text="Age (yr)", line=2, outer=TRUE )
	mtext( side=2, text="Selectivity", line=3.2, outer=TRUE )

}


plotLAA <- function()
{
	load( "data/modelObjects/Lmean_as.Rdata" )
	load( "data/modelObjects/L_ast.Rdata" )
	load( "data/modelObjects/nL_ast.Rdata" )

	par( mfcol=c(4,2), mar=c(0.5,0.5,0,0), oma=c(4,5,1,5) )

	nA <- dim(L_ast)[1]
	nS <- dim(L_ast)[2]
	yr <- as.numeric(dimnames(L_ast)[[3]])

	labs <- c("Male","Female")
	laba <- 4:7

	gap <- c(8,10,6,20)

	for( s in 1:2 )
	{
		for( a in 1:4 )
		{

			# Sample size
			plot( x=range(yr), y=c(1,1.09)*range(nL_ast[a, , ]), axes=FALSE, type="n" )
			rect( xleft=yr-0.3, xright=yr+0.3, ybottom=0, ytop=nL_ast[a,s, ],
						border=NA, col="grey85" )
			if( s==2 )
				axis( side=4, las=1 )

			# Length
			par(new=TRUE)
			yrng <- c(0.98,1.02)*range(L_ast[a, , ],na.rm=TRUE)
			plot( x=yr, y=L_ast[a,s, ], axes=FALSE, type="n",
						ylim=yrng )
			grid()
			box()

			abline( h=Lmean_as[a,s], lty=2 )

			points( x=yr, y=L_ast[a,s, ], pch=16 )

			if( a<4 )
				axis( side=1, labels=NA )
			else
				axis( side=1 )
			if( s==1 )
				axis( side=2, las=1 )
			else
				axis( side=2, labels=NA )
		
			par(font=2)
			text( x=mean(yr), y=yrng[2]-gap[a], labels=paste0(labs[s],", age ",laba[a]), cex=1.2 )
			par(font=1)

		}

	}

	mtext( side=1, text="Year", line=2, outer=TRUE )
	mtext( side=2, text="Mean length (mm; points)", line=3.2, outer=TRUE )
	mtext( side=4, text="Sample size (bars)", line=3.2, outer=TRUE )

}


plotLAA2 <- function()
{
	load( "data/modelObjects/METF_ast.Rdata" )

	par( mfcol=c(4,2), mar=c(0.5,0.5,0,0), oma=c(4,5,1,1) )

	nA <- dim(METF_ast)[1]
	nS <- dim(METF_ast)[2]
	yr <- as.numeric(dimnames(METF_ast)[[3]])

	labs <- c("Male","Female")
	laba <- 4:7

	gap <- c(8,10,6,20)

	for( s in 1:2 )
	{
		for( a in 1:4 )
		{

			yrng <- c(0.95,1.05)*range(METF_ast[a, , ],na.rm=TRUE)
			plot( x=yr, y=METF_ast[a,s, ], axes=FALSE, type="n",
						ylim=yrng )
			grid()
			box()

			points( x=yr, y=METF_ast[a,s, ], pch=16 )

			if( a<4 )
				axis( side=1, labels=NA )
			else
				axis( side=1 )
			if( s==1 )
				axis( side=2, las=1 )
			else
				axis( side=2, labels=NA )
		
			par(font=2)
			text( x=mean(yr), y=yrng[2]-gap[a], labels=paste0(labs[s],", age ",laba[a]), cex=1.2 )
			par(font=1)

		}

	}

	mtext( side=1, text="Year", line=2, outer=TRUE )
	mtext( side=2, text="METF (mm; points)", line=3.2, outer=TRUE )

}















