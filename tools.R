#-----------------------------------------------------------------------------#
# tools.R                                                                     #
# Helper functions for integrated Yukon River Chinook run reconstruction      #
#                                                                             #
# Copyright 2023 by Landmark Fisheries Research, Ltd.                         #
#                                                                             #
# This software is provided to Essa Technologies in the hope that it will be  #
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

getParSummary <- function( folder="fits/mod1" )
{
  load(paste0(folder,"/rpt.Rdata"))
  z <- data.frame( stock = rpt$stocks,
                   UMSY = format(round(filter(rpt$sdrpt,par=="UMSY_p")$val,2),nsmall=2),
                   UMSYse = format(round(filter(rpt$sdrpt,par=="UMSY_p")$se,2),nsmall=2),
                   SMSY = round(filter(rpt$sdrpt,par=="SMSY_p")$val),
                   SMSYse = round(filter(rpt$sdrpt,par=="SMSY_p")$se),
                   recSD = format(round(filter(rpt$sdrpt,par=="recSD_p")$val,2),nsmall=2),
                   recSDse = format(round(filter(rpt$sdrpt,par=="recSD_p")$se,2),nsmall=2),
                   phi = format(round(filter(rpt$sdrpt,par=="phi_p")$val,2),nsmall=2),
                   phise = format(round(filter(rpt$sdrpt,par=="phi_p")$se,2),nsmall=2),
                   errSD = format(round(filter(rpt$sdrpt,par=="errSD_p")$val,3),nsmall=3),
                   errSDse = format(round(filter(rpt$sdrpt,par=="errSD_p")$se,3),nsmall=3),
                   mu = format(round(filter(rpt$sdrpt,par=="meanMu_p")$val,3),nsmall=3),
                   muse = format(round(filter(rpt$sdrpt,par=="meanMu_p")$se,3),nsmall=3))
  write.csv(z,file=paste0(folder,"/parSummary.csv"),row.names=FALSE)
}

getRunTimingSummary <- function( rpt, folder="." )
{
  mu_pt <- rbind(rpt$mu_pt,colMeans(rpt$mu_pt))
  x <- as.numeric(substr(rpt$years,3,3))
  x <- match(x,unique(x))
  z <- data.frame( stock = c(rpt$stocks,"Average"),
                   y1980s = format(round(rowMeans(mu_pt[ ,x==1]),1),nsmall=1),
                   y1990s = format(round(rowMeans(mu_pt[ ,x==2]),1),nsmall=1),
                   y2000s = format(round(rowMeans(mu_pt[ ,x==3]),1),nsmall=1),
                   y2010s = format(round(rowMeans(mu_pt[ ,x>3]),1),nsmall=1) )
  write.csv(z,file=paste0(folder,"/runTimingSummary.csv"),row.names=FALSE)
}

getNLLSummary <- function( rpt )
{
  c( C=sum(rpt$nllC_rht), xC=sum(rpt$nllxC_rht), xB=rpt$nllxB, R=sum(rpt$nllR_py),
     MR=rpt$weightI*rpt$nllMR, nlp=rpt$nlp, I=sum(rpt$nllI_tg), P=sum(rpt$nllP_g),
     total=rpt$objFun )
}

loadModel <- function( name )
{
  if(name %in% names(getLoadedDLLs()))
    dyn.unload(dynlib(name))          # unlink the C++ code if already linked
  compile(paste0(name,".cpp"),flags = "")
  dyn.load(dynlib(name))          # Dynamically link the C++ code
}

rnbinom2 <- function( n, mu, sd )
{
  prob <- mu/(sd*sd)
  size <- mu*prob/(1-prob)
  rnbinom( n=n, size=size, prob=prob )
}

dnbinom2 <- function( x, mu, sd )
{
  prob <- mu/(sd*sd)
  size <- mu*prob/(1-prob)
  dnbinom( x=x, size=size, prob=prob )
}

dinvGamma <- function( x, a1=2, a2=1 )
  a1*log(a2) - lgamma(a1) - (a1+1)*log(x) - a2/x;

standardize <- function( x )
  x/sum(x,na.rm=1)

# Standardize columns of matrix
standCol <- function(x)
  apply( x, 2, standardize )

# Bounded logit transformation
logit<-function( x, lb=0, ub=1 )
{
  xB <- (x-lb)/(ub-lb)    # x bounded
  log(xB/(1.-xB))
}

# Bounded inverse logit transformation
invlogit<-function( x, lb=0, ub=1 )
  lb + (ub-lb)*exp(x)/(1.+exp(x))

plotbg <- function(col=rgb(235,235,235,maxColorValue=255))
{
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=col)
  grid( col="white", lty=1 )
}

# Pearson selectivity - modified from B. Staton's script
# sigma > 0, theta > 0
getPearsonSel <- function( RLM,
                           lambda      = -0.547,
                           theta       = 0.622,
                           sigma       = 0.204,
                           tau         = 1.920,
                           standardize = TRUE )
{
 
  # enforce parameter constraints
  if(any(sigma <= 0))
    stop("sigma must be > 0")
  if(any(theta <= 0))
  stop("theta must be > 0")
  
  t1 <- (1 + lambda^2/(4 * theta^2))^theta
  t2 <- RLM - (sigma * lambda)/(2 * theta) - tau
  t3 <- (1 + t2^2/sigma^2)^-theta
  t4 <- exp(-lambda * (atan(t2/sigma) + atan(lambda/(2 * theta))))
  v  <- t1 * t3 * t4
  
  # standardize the output so only one age/sex is fully vuln for a gear
  if(standardize)
    v <- v/max(v)
  
  return(v)
}


# Map lower triangle to upper triangle
mirrorMatrix <- function(x)
{
  X <- nrow(x)
  for( s in 1:(X-1) )
    x[s,(s+1):X] <- x[(s+1):X,s]
  x
}

# lisread()
# lisread: Function to read a list of data objects from a file.
# The initial characters "##" denote a comment line (ignored).
# The initial characters "# " denote a variable name.
# All other lines must contain scalars or vectors of numbers.
# Furthermore, all rows in a matrix must contain the same number of
# columns. Row and column vectors are not converted to matrices.
#
# fname  : File name.
# quiet  : If true, shut up about reporting progress.
# result : List object with components named in the file.

# Original functions courtesy of Jon Schnute.
# Modifications by A.R. Kronlund.
lisread <- function( fname,quiet=TRUE )
{
  lis2var <- function( x )
  {
    # lis2var: Makes global variables from the components of a list
    # x      : list object with named components.
    # result : global variables with names and contents extracted from x.

    namx <- names( x )
    nx <- length( namx )
    if (nx > 0) for (i in 1:nx)
    {
      if (namx[i] != "")
      {
        cmd <- paste( namx[i],"<<- x[[i]]" )
        eval( parse(text=cmd) )
      }
    }
    namx[namx != ""]
  }

  # numvecX functions:
  #
  # Function to convert a single string with white spaces into a numeric
  # vector. The string is parsed into separate components and converted
  # to char and numeric. A direct conversion to numeric fails.

  numvec <- function( x )
  {
    # Deprecated.
    xp <- parse( text=x,white=T )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec2 <- function( x )
  {
    # Patch required for S6.0 where white=T option is defunct in parse.
    # Deprecated:  xp <- parse( text=x,white=T )
    # ARK 30-Oct-03: R text connections get used up, must open/close.
    tc <- textConnection( x )
    xp <- scan( tc )
    close( tc )
    xc <- as.character( xp )
    as.numeric( xc )
  }

  numvec3 <- function( x,quiet )
  {
    # ARK 12-Jan-06: Patch to read characters because Rashmi asked nicely.
    # This is a largely untested hack, no expressed or implied warrantee.
    tc <- textConnection( x )
    xp <- scan( tc, what="character",quiet=quiet )
    close( tc )
    xc <- as.character( xp )
    if ( !all(is.na(as.numeric(xc))) )
      xc <- as.numeric( xc )
    xc
  }

  #------------------------------------------------------------------#

  file <- scan( fname, what=character(), sep="\n", quiet=quiet )

  f2 <- file[ regexpr("##",file)==-1 ]           # remove comments
  nf2 <- length( f2 )                            # number of lines
  llab <- regexpr( "#",f2 )==1                   # identifies label lines
  vlab <- substring( f2[llab],3 )                # variable labels

  # ARK 30-Oct-03 R does not coerce logical to character for grep.
  ilab <- grep( "TRUE",as.character(llab) )      # label indices

  nvar <- length( vlab )                         # number of variables

  if( nvar==1 )
    nrow <- c( nf2 + 1 ) - ilab - 1
  else
    nrow <- c( ilab[2:nvar],nf2+1 ) - ilab - 1     # number of rows in var i
  zout <- list( NULL )

  for ( i in 1:nvar )
  {
    i1 <- ilab[i] + 1
    i2 <- i1 + nrow[i] - 1                       # range of lines for var i
    zstr <- paste(f2[i1:i2],collapse=" ")
#    zvec <- numvec2(zstr)                        # numeric vector
    zvec <- numvec3(zstr,quiet)                  # numeric or character vector

    nz <- length(zvec)
    zrow <- nrow[i]
    zcol <- nz / zrow                            # dimensions
    if ( (zrow>1) & (zcol>1) )                   # a true matrix
      zvec <- matrix( zvec,nrow=zrow,ncol=zcol,byrow=T )
    zout[[i]] <- zvec
    if ( !quiet )
      cat( "vlab = ", vlab[i], "\n" )
  }
  names(zout) <- vlab
  zout
}

# .readParFile   (reads an ASCII file with 1 comment line, header, data frame)
# Purpose:      Reads an ASCII file: 1 comment, 1 header, space-delimited
#               data frame usually containing columns "parameter" and "value".
# Parameters:   parFile is a character string indicating the input file.
# Returns:      result, a data frame.
# Source:       A.R. Kronlund
.readParFile <- function( parFile="inputFile.par" )
{
  # Read the file and store as a dataframe.
  result <- read.table( file=parFile, as.is=TRUE, header=TRUE, skip=1,
                        quote="",sep=" " )
  result
}

.createList <- function( obj )
{
  # Input  a data frame with columns "parameter" and "value".
  # Output a list with elements named as parameter and values in "value".

  result <- list()

  # Shut off whining, then coerce to numeric to let NA indicate non-numerics.
  options( warn=-1 )
  numericVal <- as.numeric( obj[,"value"] )
  options( warn=0 )

  for ( i in 1:nrow(obj) )
  {
    # Value is numeric, build the parse string.
    if ( !is.na(numericVal[i]) )
      listText <- paste( "result$",obj[i,"parameter"],"=",
                    obj[i,"value"],sep="" )
    # Value is character, build the parse string.
    else
      listText <- paste( "result$",obj[i,"parameter"],"=",
                  obj[i,"value"], sep="" )

    # ARK: At one point I had this code, presumably to handle a different
    #      input format convention, perhaps assuming "value" was all character.
    #                   sQuote(obj[i,"value"]),sep="" )
    
    # Evaluate the parsed string.
    eval( parse( text=listText ) )
  }
  result
}