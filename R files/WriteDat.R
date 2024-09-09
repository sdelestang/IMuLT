WriteDataFile <- function(Report,TheData,DatFile,Nsim)
{
 # Extract stuff needed
 Narea <- TheData$Narea
 Year1 <- TheData$Year1
 Year2 <- TheData$Year2
 Nyear <- Year2-Year1+1
 Nstep <- TheData$Nstep
 Nage <- TheData$Nage
 Nsex <- TheData$Nsex
 Nfleet <- TheData$Nfleet
 Nlen <- TheData$Nlen
 MaxLen <- max(Nlen)
 
 if (Nsim == -1) write("# Original Data File\n",DatFile)
 if (Nsim == 0) write("# Expected values\n",DatFile)
 if (Nsim > 0) write(paste("# Psuedo data set",Nsim,"\n"),DatFile)
 
 write(paste("# First year of the assessment\n",Year1,sep=""),DatFile,append=T)
 write(paste("# Last year of the assessment\n",Year2,sep=""),DatFile,append=T)
 write(paste("# Burn-in\n",TheData$BurnIn,sep=""),DatFile,append=T)
 write(paste("# Time steps per year\n",Nstep,sep=""),DatFile,append=T)
 write(paste("# Number of areas of data included in the file\n",Narea,sep=""),DatFile,append=T)
 write(paste("# Number of sexes\n",Nsex,sep=""),DatFile,append=T)
 write(paste("# Number of ages\n",Nage-1,sep=""),DatFile,append=T)
 write(paste("# Number of fleets\n",Nfleet,sep=""),DatFile,append=T)
 
 write(paste("# Number of size-classes (males then females)"),DatFile,append=T)
 write(Nlen,DatFile,append=T)
 write("# The Time steps",DatFile,append=T)
 write(paste(c("#",seq(from=1,to=Nstep))),DatFile,append=T,ncol=Nstep+1,)
 write(Data$TimeStepLen[1,],DatFile,ncol=Nstep,append=T)
 write(paste("# Loop counter for initial conditions\n",TheData$Num_Iteration,sep=""),DatFile,append=T)
 write(paste("# Years over which to tune\n",TheData$Tune_Years,sep=""),DatFile,append=T)
 
 write(paste("\n# Lower Length Bins (one more than number of size-classes"),DatFile,append=T)
 for (Isex in 1:Nsex)
  write(TheData$LowLenBin[Isex,1:(Nlen[Isex]+1)],DatFile,append=T,ncol=Nlen[Isex]+1)
 
 # Catches
 write(paste("\n# Catch data\n",Nfleet*Nyear*Nstep,"\n# Year step fleet catch",sep=""),DatFile,append=T)
 CatOut <- matrix(0,nrow=Nfleet*Nyear*Nstep,ncol=4)
 Ipnt <- 0
 for (Iyear in 1:Nyear)
  for (Istep in 1:Nstep)
   for (Ifleet in 1:Nfleet)
    { Ipnt <- Ipnt + 1; CatOut[Ipnt,] <- c(Iyear+Year1-1,Istep,Ifleet,Data$Catch[Iyear,Istep,Ifleet])  }
  write(t(CatOut),DatFile,ncol=4,append=T)
  
 # Cpue
 NcpueDataSeries <- TheData$NcpueDataSeries
 write("\n# Index data",DatFile,append=T)
 write(paste("# Number of cpue sets\n",NcpueDataSeries,sep=""),DatFile,append=T)
 write("# Type of index (1=weight;2=numbers)",DatFile,append=T)
 write(TheData$IndexType,DatFile,append=T,ncol=NcpueDataSeries)
 write("# Treatment of sigma",DatFile,append=T)
 write(TheData$FixSigmaCpue,DatFile,append=T,ncol=NcpueDataSeries)
 write("# Treatment of q",DatFile,append=T)
 write(TheData$TreatQcpue,DatFile,append=T,ncol=NcpueDataSeries)
 write("# Environmental index",DatFile,append=T)
 write(TheData$EnvIndCpue,DatFile,append=T,ncol=NcpueDataSeries)
 write(paste("# Minimum sigma\n",TheData$SigmaCpueOffset,sep=""),DatFile,append=T)
 
 write(paste("# The cpue data\n",TheData$Ncpue,sep=""),DatFile,append=T)
 CpueOut <- matrix(0,nrow=TheData$Ncpue,ncol=7)
 write("#CpueInd #Fleet	#Sex	#Year	#Step	#Index	        #CV",DatFile,append=T)
 for (Idata in 1:TheData$Ncpue)
  {
   CpueOut[Idata,1:5] <- TheData$IndexI[Idata,] 
   IdataSet <- TheData$FixSigmaCpue[TheData$IndexI[Idata,1]+1]
   SigmaUse <- Report$SigmaCpue[IdataSet+1]*TheData$IndexR[Idata,2]
   CpueOut[Idata,4] <- CpueOut[Idata,4] + Year1
   CpueOut[Idata,7] <- TheData$IndexR[Idata,2]
   if (Nsim == -1) CpueOut[Idata,6] <- TheData$IndexR[Idata,1]
   if (Nsim == 0) CpueOut[Idata,6] <- Report$PredCpue[Idata,1]
   if (Nsim > 0) CpueOut[Idata,6] <- Report$PredCpue[Idata,1]*exp(rnorm(1,0,SigmaUse)-SigmaUse^2/2.0)
  }
 write(t(CpueOut),DatFile,ncol=7,append=T)
 
 # Catch-in-numbers
 NcatchDataSeries <- TheData$NcatchDataSeries
 write("\n# Numbers data",DatFile,append=T)
 write(paste("# Number of catch data sets\n",NcatchDataSeries,sep=""),DatFile,append=T)
 write("# Treatment of catch in numbers series",DatFile,append=T)
 write(TheData$FixSigmaCatchN,DatFile,append=T,ncol=NcatchDataSeries)
 write(paste("# Minimum sigma\n",TheData$SigmaCatchNOffset,sep=""),DatFile,append=T)
 write(paste("# The numbers data\n",TheData$Nnumbers,sep=""),DatFile,append=T)
 NumbersOut <- matrix(0,nrow=TheData$Nnumbers,ncol=6)
 write("#Group	#Fleet	year	step	catch	CV",DatFile,append=T)
 for (Idata in 1:TheData$Nnumbers)
  {
   NumbersOut[Idata,1:4] <- TheData$NumbersI[Idata,]
   IdataSet <- TheData$FixSigmaCatchN[TheData$NumbersI[Idata,1]+1]
   SigmaUse <- Report$SigmaNumbers[IdataSet+1]*TheData$NumbersR[Idata,2]
   NumbersOut[Idata,3] <- NumbersOut[Idata,3] + Year1
   NumbersOut[Idata,6] <- TheData$NumbersR[Idata,2]
   if (Nsim == -1) NumbersOut[Idata,5] <- TheData$NumbersR[Idata,1]
   if (Nsim == 0) NumbersOut[Idata,5] <- Report$PredNumbers[Idata,1]
   if (Nsim > 0) NumbersOut[Idata,5] <- Report$PredNumbers[Idata,1]*exp(rnorm(1,0,SigmaUse)-SigmaUse^2/2.0)
  }
 write(t(NumbersOut),DatFile,ncol=6,append=T)
  
 # length data 
 NlenComp <- TheData$NlenComp
 write(paste("\n# Length compostion\n",NlenComp),DatFile,append=T)
 write("#Fleet	Sex	SEASON	StpLen	Ind",DatFile,append=T)
 LengthOut <- matrix("",nrow=NlenComp,ncol=5+MaxLen)
 for (Idata in 1:NlenComp)
  {
   LengthOut[Idata,1:4] <- TheData$LenCompI[Idata,]+c(0,0,Year1,0)
   LengthOut[Idata,5] <- TheData$Stage1W[Idata]
   NlenC <- Nlen[TheData$LenCompI[Idata,2]+1]
   Ifleet <- TheData$LenCompI[Idata,1]+1
   Istep <- TheData$LenCompI[Idata,4]+1
   Isex <- TheData$LenCompI[Idata,2]+1
   Overdisp <- TheData$LambdaLength2[Ifleet,Istep,Isex]
   if (Nsim == -1) LengthOut[Idata,6:(5+NlenC)] <- TheData$LenCompR[Idata,1:NlenC]*TheData$Stage1W[Idata]
   if (Nsim == 0) LengthOut[Idata,6:(5+NlenC)] <- Report$PredLengthComp[Idata,1:NlenC]*TheData$Stage1W[Idata]
   if (Nsim > 0)LengthOut[Idata,6:(5+NlenC)] <- as.numeric(rmultinom(1,round(TheData$Stage1W[Idata]*Overdisp),Report$PredLengthComp[Idata,1:NlenC]))/Overdisp
  }
 write(t(LengthOut),DatFile,ncol=5+MaxLen,append=T)
 
 # larval data
 write("\n# Larval index",DatFile,append=T)
 write(paste("# Likelihood for larval data (0=Lognormal; otherwise normal)\n",TheData$LarvalLikeOpt,sep=""),DatFile,append=T)
 write(paste("# Delay from purulus stage and entry to the model\n",TheData$Larval_Offset,sep=""),DatFile,append=T)
 write(paste("# Number of data points\n",TheData$NLarvalData,sep=""),DatFile,append=T)
 write("# Area Year Index SD",DatFile,append=T)
 LarvalOut <- matrix(0,nrow=TheData$NLarvalData,ncol=4)
 if (TheData$NLarvalData > 0)
  for (Idata in 1:TheData$NLarvalData)
   {
    LarvalOut[Idata,1:2] <- TheData$Lar_dataI[Idata,]
    LarvalOut[Idata,2] <- LarvalOut[Idata,2] + Year1 - TheData$BurnIn
    LarvalOut[Idata,4] <- TheData$Lar_dataR[Idata,2]
    SigmaUse <- TheData$Lar_dataR[Idata,2]
    if (Nsim == -1) LarvalOut[Idata,3] <- TheData$Lar_dataR[Idata,1]
    if (Nsim == 0) LarvalOut[Idata,3] <- Report$PredLarval[Idata,1]
    if (Nsim > 0) LarvalOut[Idata,3] <- Report$PredLarval[Idata,1]*exp(rnorm(1,0,SigmaUse)-SigmaUse^2/2.0)
   }
 write(t(LarvalOut),DatFile,ncol=4,append=T)
  
 # Envronmental data 
 write("\n# Environmental Data",DatFile,append=T)
 write(paste("# Number of environmental series\n",TheData$NenvSeries,sep=""),DatFile,append=T)
 write("# Years of data per series",DatFile,append=T)
 write(rep(Nyear+TheData$BurnIn+1,TheData$NenvSeries),DatFile,append=T,ncol=TheData$NenvSeries)
 EnvOut <- matrix(0,nrow=(Nyear+TheData$BurnIn+1)*TheData$NenvSeries,ncol=2)
 Ipnt <- 0
 if (TheData$NenvSeries > 0)
  for (Iseries in 1:TheData$NenvSeries)
   for (Iyear in 1:(Nyear+TheData$BurnIn+1))
    {
     Ipnt <- Ipnt + 1
     EnvOut[Ipnt,1] <- Iyear+Year1-TheData$BurnIn-1
     EnvOut[Ipnt,2] <- TheData$EnvData[Iyear,Iseries+1]
    }
 write(t(EnvOut),DatFile,ncol=2,append=T)
    
 
 
 
  
 
}
