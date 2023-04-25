#-----------------------------------------------------------------------------#
# ypr.R                                                                       #
# YPR script for integrated Yukon River Chinook run reconstruction            #
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


calcYPR <- function()
{
	#mods <- c("mod1-spawners","mod2-eggNum","mod3-eggMass")
	mods <- c("mod1","mod2","mod3")
	nI <- length(mods)

	load(paste0("fits/",mods[1],"/rpt.Rdata"))

	C_ipm <- array( data=NA, dim=c(nI,rpt$nP,rpt$nM) )
	S_ipm <- array( data=NA, dim=c(nI,rpt$nP,rpt$nM) )
	F_ipm <- array( data=NA, dim=c(nI,rpt$nP,rpt$nM) )

	for( i in 1:length(mods) )
	{
		load(paste0("fits/",mods[i],"/rpt.Rdata"))
		z_as <- apply( rpt$z_ast, 1:2, mean )

		for( m in 1:rpt$nM )
		{	
			v_as   <- apply( rpt$v_astm[ , , ,m], 1:2, mean )
			v_as   <- v_as/max(v_as)
	
			for( p in 1:rpt$nP )
			{
				eta_as <- apply( rpt$eta_aspy[ , ,p, ], 1:2, mean )
				eta_as <- eta_as/sum(eta_as)
		
				calcCeq <- function(Fmax)
				{
					U_as <- 1-exp(-Fmax*v_as)
					reproOutputPerSpawner <- sum((1-U_as)*z_as*eta_as)
					alpha <- exp( log(rpt$alpha_p[p]) + (rpt$recSD_p[p]^2)/(2*(1-rpt$phi_p[p]^2)) )
					Req <- log( alpha*reproOutputPerSpawner ) / (rpt$beta_p[p]*reproOutputPerSpawner)
					Neq_as <- Req * eta_as
					Ceq_as <- Neq_as * U_as
					Seq_as <<- Neq_as * (1-U_as)
					sum(Ceq_as)
				}
		
  			opt <- optimize( calcCeq, c(0,5), maximum=TRUE )
				
  			C_ipm[i,p,m] <- opt$objective
  			F_ipm[i,p,m] <- opt$maximum
  			S_ipm[i,p,m] <- sum(Seq_as)
	
  			#x <- seq(0,5,by=0.01)
  			#y <- mapply(calcCeq,x)
  			#plot(x,y)
  		}
		}
	}

	ypr <- list( C_ipm=C_ipm, F_ipm=F_ipm, S_ipm=S_ipm, stocks=rpt$stocks )
	save( ypr, file="ypr.Rdata" )
}

plotYPR <- function()
{
	calcYPR()
	load("ypr.Rdata")	

	nI <- dim(ypr[[1]])[1]
	nP <- dim(ypr[[1]])[2]
	nM <- dim(ypr[[1]])[3]

	cols <- brewer.pal(nI,"Set1")
	agrey <- alpha("grey85",alpha=0.3)

	p <- 1:nP
	pev <- -1:(nP+2)
	pev <- pev[pev%%2==0]

	jtr <- seq(-0.15,0.15,length.out=nI)

	par( mfrow=c(2,nM-1), mar=c(0.5,0,0,0), oma=c(11,5,2,1) )

	Smax <- max(ypr$S_ipm)*1e-3
	Cmax <- max(ypr$C_ipm)*1e-3

	for( m in c(1,nM) )
	{
		S_ip <- ypr$S_ipm[ , ,m]*1e-3

		plot( x=c(0.5,nP+0.5), y=c(0,1.1*Smax), type="n", axes=FALSE )
		axis( side=1, labels=NA, at=1:nP )
		if( m==1 )
			axis( side=2, las=1 )
		abline( h=seq(0,Smax,by=10), col="grey90" )

		rect( xleft=pev-0.5, xright=pev+0.5, ybottom=-2*min(S_ip), ytop=2*max(S_ip), border=NA, col=agrey )
		box()
		for( i in 1:nI )
			points( x=1:nP+jtr[i], y=S_ip[i, ], pch=14+i, col=cols[i] )
	
		if( m==1 )
			legend( x="topleft", bty="n", col=cols, pch=seq(15,by=1,length=nI),
							legend=c("S model","E model","EM model") )

	}

	for( m in c(1,nM) )
	{
		C_ip <- ypr$C_ipm[ , ,m]*1e-3

		plot( x=c(0.5,nP+0.5), y=c(0,Cmax), type="n", axes=FALSE )
		axis( side=1, labels=NA )
		if( m==1 )
			axis( side=2, las=1 )
		abline( h=seq(0,Cmax,by=10), col="grey90" )
		rect( xleft=pev-0.5, xright=pev+0.5, ybottom=-2*min(C_ip), ytop=2*max(C_ip), border=NA, col=agrey )
		box()
		for( i in 1:nI )
			points( x=1:nP+jtr[i], y=C_ip[i, ], pch=14+i, col=cols[i] )
		axis( side=1, labels=ypr$stocks, at=1:nP, las=2 )

	}

	mtext( side=2, text="Catch (1000s)                        Escapement (1000s)", outer=TRUE, line=3.5 )
	mtext( side=3, text="6 inch mesh                                        8 inch mesh", outer=TRUE, line=0.5 )

}











