library(parallel)
library(coda)
library(rjags)
library(runjags)
library(compute.es)

#header functions to help with calculations

nCores = detectCores() 
if ( !is.finite(nCores) ) { nCores = 1 } 
if ( nCores > 4 ) { 
  nChainsDefault = 4  # because JAGS has only 4 rng's.
  runjagsMethodDefault = "parallel"
}
if ( nCores == 4 ) { 
  nChainsDefault = 3  # save 1 core for other processes.
  runjagsMethodDefault = "parallel"
}
if ( nCores < 4 ) { 
  nChainsDefault = 3 
  runjagsMethodDefault = "rjags" # NOT parallel
}

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ,
                     saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("skyblue","black","royalblue","steelblue")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
}


#data for meta-analysis of beta blockers (Gelman 7.7.3 using data from section 5.6)
meta = matrix( c(
  1,     3, 39,        3, 38,
  2,     14, 116,      7, 114, 
  3,     11, 93,         5, 69, 
  4,     127, 1520,  102, 1533, 
  5,     27, 365,      28, 355, 
  6,     6, 52,          4, 59, 
  7,     152, 939,     98, 945, 
  8,     48, 471,      60, 632, 
  9,     37, 282,      25, 278, 
  10,    188, 1921,  138, 1916, 
  11,    52, 583,      64, 873, 
  12,    47, 266,      45, 263, 
  13,    16, 293,       9, 291, 
  14,    45, 883,      57, 858, 
  15,    31, 147,      25, 154, 
  16,    38, 213,      33, 207, 
  17,    12, 122,      28, 251, 
  18,    6, 154,        8, 151, 
  19,    3, 134,        6, 174, 
  20,    40, 218,      32, 209, 
  21,    43, 364,      27, 391, 
  22,    39, 674,      22, 680 
), ncol=5, byrow=TRUE )

colnames(meta) = c( "StudyID" , "ControlDeaths" , "ControlTotal" , "TreatedDeaths" , "TreatedTotal" )
Meta = data.frame(meta)
# Package the data for JAGS:
dataList = list(
  nS = length(Meta$StudyID) ,
  zC = Meta$ControlDeaths ,
  nC = Meta$ControlTotal , 
  zT = Meta$TreatedDeaths ,
  nT = Meta$TreatedTotal) 


# Define the JAGS model:
#
# thetaC[s] is success probability in control group, study s.
# thetaT[s] is success probability in treatment group, study s.
# rho[s] is the difference of log-odds between groups:
#   rho[s] = logit(thetaT[s]) - logit(thetaC[s])
#   Re-arranged, the equation expresses the relation of thetaT[s] to thetaC[s]:
#   thetaT[s] = logistic( rho[s] + logit(thetaC[s]) )
#   This relation is a "natural" way to represent the dependency of the 
#   probabilities between groups because the relation is (i) symmetric with
#   respect to what outcome is defined as success because logit(.5+p) =
#   -logit(.5-p), and (ii) symetric with respect to which group is defined as
#   the treatment by reversing the sign of rho. 
#   Note that rho[s] is the so-called log odds ratio across groups: 
#   rho[s] = log( (thetaT[s]/(1-thetaT[s])) / (thetaC[s]/(1-thetaC[s])) )
# thetaComega is modal thetaC[s] across control studies.
# thetaCkappa is concentration (consistency) of thetaC[s] across studies.
# rhoMu is the mean rho[s] across studies.
# rhoSD is the standard deviation rho[s] across studies.
#

modelString = "
model {
for ( s in 1:nS ) {
zC[s] ~ dbin( thetaC[s] , nC[s] )
zT[s] ~ dbin( thetaT[s] , nT[s] )
thetaC[s] ~ dbeta( thetaComega*(thetaCkappa-2)+1 ,
(1-thetaComega)*(thetaCkappa-2)+1 )
thetaT[s] <- ilogit( rho[s] + logit( thetaC[s] ) ) # ilogit is logistic
rho[s] ~ dnorm( rhoMu , 1/rhoSD^2 )
}
# Prior rhoMeta, rhoSD:
rhoMu ~ dnorm( 0 , 1/10^2 )
rhoSD ~ dgamma(1.64,0.64) # mode=1,sd=2
# Prior on thetaComega and thetaCkappa:
thetaComega ~ dbeta(1.01,1.01) 
thetaCkappa <- thetaCkappaMinusTwo + 2
thetaCkappaMinusTwo ~ dgamma(2.618,0.162) # mode=10 , sd=10
# Derived variables:
for ( s in 1:nS ) { 
mu[s] <- thetaT[s] / thetaC[s]  # risk ratio
}
thetaTomega <- ilogit( rhoMu + logit( thetaComega ) )
muMeta <- thetaTomega / thetaComega  # risk ratio
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Variables to monitor:
parameters = c( "thetaC" , "thetaComega" , "thetaCkappa" , 
                "thetaT" , "thetaTomega" ,
                "rho" , "rhoMu" , "rhoSD" , 
                "mu" , "muMeta" )
# Initial values for MCMC chains:
initsList = list( 
  thetaC=rep( sum(dataList$zC)/sum(dataList$nC) , dataList$nS ) ,
  thetaComega=sum(dataList$zC)/sum(dataList$nC) , 
  rho = rep( 0 , dataList$nS ) ,
  rhoMu = 0 ,
  rhoSD = 1 ,
  thetaCkappaMinusTwo = 5
)

# Run the chains:
adaptSteps = 1000 
burnInSteps = 1000
numSavedSteps = 30000
thinSteps=10
nChains = 3

#Run model and save posterior simulations
runJagsOut <- run.jags( method=c("rjags","parallel")[2] ,
                        model="TEMPmodel.txt" , 
                        monitor=parameters , 
                        data=dataList ,  
                        inits=initsList , 
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps , 
                        sample=ceiling(numSavedSteps/nChains) ,
                        thin=thinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )

codaSamples = as.mcmc.list( runJagsOut )

for ( parName in c( "thetaC[1]" , "thetaComega" , "thetaCkappa" , 
                    "thetaT[1]" , "thetaTomega" ,
                    "rho[1]" , "rhoMu" , "rhoSD" ,
                    "mu[1]" , "muMeta"  ) ) {
  diagMCMC( codaObject=codaSamples , parName=parName )
}
save( codaSamples , file=paste0("BDA hw","-Mcmc.Rdata") )
mcmcMat = as.matrix(codaSamples)

# Graphical view of selected variables:
for ( parName in c( "thetaC[1]" , "thetaComega" , "thetaCkappa" , 
                    "thetaT[1]" , "thetaTomega" ,
                    "rho[1]" , "rhoMu" , "rhoSD" ,
                    "mu[1]" , "muMeta"  ) ) {
  diagMCMC( codaObject=codaSamples , parName=parName )
}

( summaryMCMC = summary( runJagsOut )) 

summaryMCMC2 = data.frame(summaryMCMC)

logit = function(a){
  return(log((a/(1-a))))
}


#### MLE frequentist model:
Meta$pc = Meta$ControlDeaths/Meta$ControlTotal
Meta$pt = Meta$TreatedDeaths/Meta$TreatedTotal

library(lme4)
freq_m <- glmer(cbind(TreatedDeaths, ControlDeaths) ~  (1 | StudyID),data = Meta, family = binomial(link = "logit"))

#compute diagnostic statistics: AIC, DIC, WAIC, and LOO-CV:
AIC = summary(freq_m)$AIC[1]
print(AIC)
#AIC = 134

DIC = extract(runJagsOut, "DIC")
#DIC stat = 285

library(loo)
mcall <- rbind(runJagsOut$mcmc[[1]], runJagsOut$mcmc[[2]], runJagsOut$mcmc[[3]])

WAIC = waic(mcall)
WAIC_stat = -2*(WAIC$estimates[1] - WAIC$estimates[2])
print(WAIC_stat)
#WAIC = 156.06

LOO = loo(mcall)
LOO_stat = -2*LOO$estimates[3]
print(LOO_stat)
#LOO = 143.67
