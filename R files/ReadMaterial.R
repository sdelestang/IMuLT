#### Test whether index is correct - stop function if it is not an show line
testindex <- function(Data,Row,Col,Ind) {
  if(as.numeric(Data[Row,Col])!=(Ind-1)){
    stop(paste("Col",Col,"is wrong in selectivity pointer matrix, line",Row),call. = F) }}

## Fancy way to choose the model run
csemod <- function(x){
  mod <-  as.numeric(dlg_input(c('Choose a model:',paste(1:length(x), x, sep=(" : ") )), 1)$res)
  if (!length(mod)) {# The user clicked the 'cancel' button
    cat(paste("OK, the default model is",x[1],"\n"))
  } else {
    cat(paste("Model", x[mod], "has been chosen"), "\n")
  }
  return(x[mod])}

asnum <- function(x) suppressWarnings(!is.na(as.numeric(x)))

## Runs through data loking for NAs and data errors
isnafunc <- function(x){
  if(length(x)>0) {
    if(is.na(sum(x))){ 
      print(paste('There is a NA in', names(Data)[i])) }
  } else {print(paste('There is no data for', names(Data)[i]))   }
}

ReadStarterFile <- function(StarterFile)
 {
  ReturnObj <- NULL
  ReturnObj$DataFileName <- StarterFile[1,1]
  ReturnObj$ControlFileName <- StarterFile[2,1]
  ReturnObj$SelexFileName <- StarterFile[3,1]
  ReturnObj$RetainFileName <- StarterFile[4,1]
  ReturnObj$RecruitFileName <- StarterFile[5,1]
  ReturnObj$GrowthFileName <- StarterFile[6,1]
  ReturnObj$MoveFileName <- StarterFile[7,1]
  ReturnObj$TagFileName <- StarterFile[8,1]
  ReturnObj$PropFFileName <- StarterFile[9,1]
  ReturnObj$ProjectionsFileName <- StarterFile[10,1]
  Index <- MatchTable(StarterFile,Char2="#",Char3="Stop",Char4="after") 
  ReturnObj$MaxPhase <- as.numeric(StarterFile[Index,1])
  return(ReturnObj)
 }  

# ===================================================================================


ReadGeneralFile <- function(DataFile)
 {
  # Read in the data file
  Index <- MatchTable(DataFile,Char1="#",Char2="First",Char3="year"); Year1 <- as.numeric(DataFile[Index+1,1])
  Index <- MatchTable(DataFile,Char1="#",Char2="Last",Char3="year"); Year2 <- as.numeric(DataFile[Index+1,1])
  Nyear <- Year2-Year1+1
  write(paste("Year1",Year1),EchoFile,append=T)
  write(paste("Year2",Year2),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Maximum",Char3="projection"); MaxProjYr <- as.numeric(DataFile[Index+1,1])
  Index <- MatchTable(DataFile,Char1="#",Char2="Time",Char3="steps"); Nstep <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of time steps",Nstep),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="areas"); Narea <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of areas",Narea),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Burn-in"); BurnIn <- max(as.numeric(DataFile[Index+1,1:Narea]), na.rm=T)
  write(paste("Max length of burn-in",BurnIn),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Burn-in"); BurnInVec <- as.numeric(DataFile[Index+1,1:Narea]); BurnIn <- max(BurnInVec)
  write(paste("Area specific length of burn-in",BurnInVec),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="sexes"); Nsex <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of sexes",Nsex),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="ages"); Nage <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of ages",Nage),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="fleets"); Nfleet <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of fleets",Nfleet),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="size-classes"); Nlen <- as.numeric(DataFile[Index+1,c(1,2)]); MaxLen <- max(Nlen)
  write(paste("Number of size-classes",Nlen),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="The",Char3="Time"); TimeStepLen <- as.numeric(DataFile[Index+2,1:Nstep]);
  write("Length of each time-step",EchoFile,append=T)
  write(TimeStepLen,EchoFile,append=T,ncol=1000)
  TimeStepLenA <- TimeStepLen; for (Iyear in 1:(Nyear+MaxProjYr)) TimeStepLenA <- rbind(TimeStepLenA,TimeStepLen)
  write("Matrix of time steps",EchoFile,append=T)
  write(t(TimeStepLenA),EchoFile,append=T,ncol=Nstep)
  Index <- MatchTable(DataFile,Char1="#",Char2="Loop",Char3="counter"); Num_Iteration <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of iterations to set up the intial size-structure",Num_Iteration),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Years",Char3="over",Char4="which"); Tune_Years <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of years to tune over to set up the initial size-structure",Tune_Years),EchoFile,append=T)

  Index <- MatchTable(DataFile,Char1="#",Char2="Lower",Char3="Length"); 
  LowerLen <- matrix(0,Nsex,MaxLen+1)
  for (Isex in 1:Nsex)
    LowerLen[Isex,] <- as.numeric(DataFile[Index+Isex,1:(MaxLen+1)]);
  if (!is.na(DataFile[(Index+Isex),(MaxLen+2)])) { print("Error reading length-bin; too many inputs: Stopping"); AAA }
  write("Lower bounds on size-classes",EchoFile,append=T)
  write(t(LowerLen),EchoFile,append=T,ncol=MaxLen+1)
  
  MidLenBin <- matrix(0,Nsex,MaxLen+1)
  for (Isex in 1:Nsex)
    for (Ilen in 1:Nlen[Isex])  
      MidLenBin[Isex,Ilen] <- (LowerLen[Isex,Ilen]+LowerLen[Isex,Ilen+1])/2
  write("Midpoints of the size-classes",EchoFile,append=T)
  write(t(MidLenBin),EchoFile,append=T,ncol=MaxLen+1)
  print("READ IN THE GENERAL FILE")
  write("READ IN THE GENERAL FILE\n\n",EchoFile,append=T)
  
  ReturnObj <- NULL
  ReturnObj$Year1 <- Year1
  ReturnObj$Year2 <- Year2
  ReturnObj$MaxProjYr <- MaxProjYr
  ReturnObj$Nyear <- Nyear
  ReturnObj$Nstep <- Nstep
  ReturnObj$Narea <- Narea
  ReturnObj$Nage <- Nage
  ReturnObj$Nsex <- Nsex
  ReturnObj$Nfleet <- Nfleet
  ReturnObj$MaxLen <- MaxLen
  ReturnObj$Nlen <- Nlen
  ReturnObj$TimeStepLen <- TimeStepLenA
  ReturnObj$BurnIn <- BurnIn
  ReturnObj$BurnInVec <- BurnInVec
  ReturnObj$Num_Iteration <- Num_Iteration
  ReturnObj$Tune_Years <- Tune_Years
  ReturnObj$MidLenBin <- MidLenBin
  ReturnObj$LowLenBin <- LowerLen
  return(ReturnObj)

}

# ===================================================================================
ReadTagFile <- function(TagFile,PropFFile,GeneralSpecs,DataFile)
{
 # Tag loss rates
 Index <- MatchTable(TagFile,Char1="#Initial_tag_loss"); InitialLoss <- as.numeric(TagFile[Index+1,1]);
 Index <- MatchTable(TagFile,Char1="#Long-term_tag_loss"); TagLossRate <- as.numeric(TagFile[Index+1,1]);
 Index <- MatchTable(TagFile,Char1="#",Char2="Types",Char4="reporting"); NrepSplit <- as.numeric(TagFile[Index+1,1]);
 Index <- MatchTable(TagFile,Char1="#",Char2="Reporting",Char3="rates"); RepRate <- as.numeric(TagFile[Index+1,1:NrepSplit]);
 Index <- MatchTable(TagFile,Char1="#",Char2="Use",Char3="size"); FitTagSizes <- as.numeric(TagFile[Index+1,1:NrepSplit]);
 Index <- MatchTable(TagFile,Char1="#",Char2="Number",Char4="periods"); NtagLag <- as.numeric(TagFile[Index+1,1]);
 Index <- MatchTable(TagFile,Char1="#",Char2="Groups"); NtagGroups <- as.numeric(TagFile[Index+1,1]);
 Index <- MatchTable(TagFile,Char1="#",Char2="Group",Char3="years"); Year1Tag <- as.numeric(TagFile[Index+1,1:NtagGroups]); Year2Tag <- as.numeric(TagFile[Index+2,1:NtagGroups]);
 TagYr1 <- Year1Tag[1]; TagYr2 <- Year2Tag[NtagGroups]; NyearTags <- TagYr2-TagYr1+1;
 TagRel <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NyearTags,GeneralSpecs$Nstep,GeneralSpecs$MaxLen+1))
 Index <- MatchTable(TagFile,Char1="#",Char2="Releases");
 Ipnt <- 0
 for (Isex in 1:GeneralSpecs$Nsex)
  for (Igrp in 1:NtagGroups)
   for (Iarea in 1:GeneralSpecs$Narea)
    {
     Ipnt <- Ipnt +1
     for (Iyear in 1:NyearTags)
      for (Istep in 1:GeneralSpecs$Nstep)
       {
        Ipnt <- Ipnt + 1
        TagRel[Isex,Igrp,Iarea,Iyear,Istep,] <- as.numeric(TagFile[Index+Ipnt,6:(6+GeneralSpecs$MaxLen)]);
        # sum(TagRel[,,,,,])
       }
    }   
 #print(Index+Ipnt)
 Index <- MatchTable(TagFile,Char1="#",Char2="Recaptures")[1];
 TagRec <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NrepSplit,NyearTags,GeneralSpecs$Nstep,GeneralSpecs$MaxLen+1))
 Ipnt <- 0
 for (Isex in 1:GeneralSpecs$Nsex)
  for (Igrp in 1:NtagGroups)
   for (Iarea in 1:GeneralSpecs$Narea)
    for (Irep in 1:NrepSplit)
     {
      Ipnt <- Ipnt +1
      for (Iyear in 1:NyearTags)
       for (Istep in 1:GeneralSpecs$Nstep)
       {
        Ipnt <- Ipnt + 1
        TagRec[Isex,Igrp,Iarea,Irep,Iyear,Istep,] <- as.numeric(TagFile[Index+Ipnt,7:(7+GeneralSpecs$MaxLen)]);
        #print(TagRec[Isex,Igrp,Iarea,Irep,Iyear,Istep,])
       }
    }   
 # print(Index+Ipnt)
    
 Index <- MatchTable(TagFile,Char1="#",Char2="Recaptures")[2];
 RecapObs <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NrepSplit,NyearTags,GeneralSpecs$Nstep))  
 Ipnt <- 0
 for (Isex in 1:GeneralSpecs$Nsex)
  for (Igrp in 1:NtagGroups)
   for (Iarea in 1:GeneralSpecs$Narea)
    for (Irep in 1:NrepSplit)
     {
      Ipnt <- Ipnt +1
      for (Iyear in 1:NyearTags)
       {
        Ipnt <- Ipnt + 1
        RecapObs[Isex,Igrp,Iarea,Irep,Iyear,] <- as.numeric(TagFile[Index+Ipnt,6:(5+GeneralSpecs$Nstep)]);
        #print(TagRec[Isex,Igrp,Iarea,Irep,Iyear,Istep,])
       }
    }   
 #print(Index+Ipnt)
 # Totals 
 NrelTotal <- matrix(0,nrow=GeneralSpecs$Nsex,ncol=NtagGroups)
 NotReportedObs <- matrix(0,nrow=GeneralSpecs$Nsex,ncol=NtagGroups)
 for (Isex in 1:GeneralSpecs$Nsex)
  for (Igrp in 1:NtagGroups)
   for (Iarea in 1:GeneralSpecs$Narea)
    for (Iyear in 1:NyearTags)
     for (Istep in 1:GeneralSpecs$Nstep)
      {
       NrelTotal[Isex,Igrp] <- NrelTotal[Isex,Igrp] + TagRel[Isex,Igrp,Iarea,Iyear,Istep,1]
       for (Irep in 1:NrepSplit)
        NotReportedObs[Isex,Igrp] <- NotReportedObs[Isex,Igrp] + TagRec[Isex,Igrp,Iarea,Irep,Iyear,Istep,1]
      }
 #print(NrelTotal)
 #print(NotReportedObs)

 for (Isex in 1:GeneralSpecs$Nsex)
  for (Igrp in 1:NtagGroups)
   NotReportedObs[Isex,Igrp] <- (NrelTotal[Isex,Igrp]-NotReportedObs[Isex,Igrp])/NrelTotal[Isex,Igrp]
 #print(NotReportedObs)
 for (Isex in 1:GeneralSpecs$Nsex)
  for (Igrp in 1:NtagGroups)
   for (Iarea in 1:GeneralSpecs$Narea)
    for (Iyear in 1:NyearTags)
     for (Istep in 1:GeneralSpecs$Nstep)
      for (Irep in 1:NrepSplit)
       RecapObs[Isex,Igrp,Iarea,Irep,Iyear,Istep] = TagRec[Isex,Igrp,Iarea,Irep,Iyear,Istep,1]/NrelTotal[Isex,Igrp]

 PropRepSplit<-array(0,dim=c(NyearTags,GeneralSpecs$Nstep,GeneralSpecs$Narea,NrepSplit));
 Index <- 1
 Ipnt <- 0
 for (Iyear in 1:NyearTags)
  for (Istep in 1:GeneralSpecs$Nstep)
   for (Iarea in 1:GeneralSpecs$Narea)
    {
     Ipnt <- Ipnt + 1
     PropRepSplit[Iyear,Istep,Iarea,]  <- as.numeric(PropFFile[Index+Ipnt,4:(3+NrepSplit)]);
    }

 ReturnObj <- NULL 
 ReturnObj$InitialLoss <- InitialLoss
 ReturnObj$TagLossRate <- TagLossRate
 ReturnObj$NrepSplit <- NrepSplit
 ReturnObj$RepRate <- RepRate
 ReturnObj$FitTagSizes <- FitTagSizes
 ReturnObj$NtagLag <- NtagLag
 ReturnObj$NtagGroups <- NtagGroups
 ReturnObj$Year1Tag <- Year1Tag
 ReturnObj$Year2Tag <- Year2Tag
 ReturnObj$TagYr1 <- TagYr1
 ReturnObj$TagYr2 <- TagYr2
 ReturnObj$NyearTags <- NyearTags
 ReturnObj$TagRel <- TagRel
 ReturnObj$TagRec <- TagRec
 ReturnObj$RecapObs <- RecapObs
 ReturnObj$NrelTotal <-  NrelTotal
 ReturnObj$NotReportedObs <- NotReportedObs
 ReturnObj$PropRepSplit <- PropRepSplit
#  print(str(ReturnObj))
 return(ReturnObj)

}

# ===================================================================================


ReadDataFile <- function(DataFile,GeneralSpecs)
{
  
  # Catch data
  Index <- MatchTable(DataFile,Char1="#",Char2="Catch",Char3="data"); Ncatch  <- as.numeric(DataFile[Index+1,1]); Index <- Index + 2
  Catch <- array(0,dim=c(GeneralSpecs$Nyear+GeneralSpecs$MaxProjYr,GeneralSpecs$Nstep,GeneralSpecs$Nfleet))
  for (Icatch in 1:Ncatch)
   {
    Year <- as.numeric(DataFile[Index+Icatch,1])-GeneralSpecs$Year1+1;Step <- as.numeric(DataFile[Index+Icatch,2]);Fleet <-as.numeric(DataFile[Index+Icatch,3])   
    Catch[Year,Step,Fleet] <- as.numeric(DataFile[Index+Icatch,4])
  }  
  if (asnum((DataFile[Index+Icatch+1,4]))) { print("Error reading Catch data; too many inputs: Stopping"); AAA }
  write("Catch data by year fleet, step and fleet",EchoFile,append=T)
  write(Catch,EchoFile,append=T,ncol=GeneralSpecs$Nfleet)
  
  # Index data
  Index <- MatchTable(DataFile,Char1="#",Char2="Index",Char3="data"); 
  NcpueDataSeries <- as.numeric(DataFile[Index+2,1]);
  IndexType <- as.numeric(DataFile[Index+4,1:NcpueDataSeries]);
  FixSigmaCpue <- as.numeric(DataFile[Index+6,1:NcpueDataSeries]);
  TreatQcpue <- as.numeric(DataFile[Index+8,1:NcpueDataSeries]);
  EnvIndCpue  <- as.numeric(DataFile[Index+10,1:NcpueDataSeries]);
  EffCrIndCpue  <- as.numeric(DataFile[Index+12,1:NcpueDataSeries]);
  EffCrLag  <- as.numeric(DataFile[Index+14,1:max(EffCrIndCpue)]);
  SigmaCpueOffset <- as.numeric(DataFile[Index+16,1]);
  Index <- MatchTable(DataFile,Char1="#",Char2="The",Char3="cpue",Char4="data"); Ncpue  <- as.numeric(DataFile[Index+1,1]); Index <- Index + 2
  write(paste("Number of cpue points",Ncpue),EchoFile,append=T)
  IndexI <- matrix(0,nrow=Ncpue,ncol=5)
  IndexR <- matrix(0,nrow=Ncpue,ncol=2)
  Icpue <- 1
  for (Icpue in 1:Ncpue)
   {
    IndexI[Icpue,] <- as.numeric(DataFile[Index+Icpue,1:5]) - c(0,1,1,0,1)  ## This allows R indexing to be used in the input file
    IndexI[Icpue,4] <- IndexI[Icpue,4] - GeneralSpecs$Year1
    IndexR[Icpue,] <- as.numeric(DataFile[Index+Icpue,6:7])
  }  
  if (asnum(DataFile[Index+Icpue+1,6])) { print("Error reading CPUE data; too many inputs: Stopping"); AAA }
      
  write("Cpue data",EchoFile,append=T)
  write(t(cbind(IndexI,IndexR)),EchoFile,append=T,ncol=7)
  
  NQpars <- 0 
  for (IdataSet in 1:NcpueDataSeries)
   if (EnvIndCpue[IdataSet] != 0)  NQpars <- NQpars + 1
  write(paste("Number of Q-related parameters",NQpars),EchoFile,append=T)

  # Numbers index data
  Index <- MatchTable(DataFile,Char1="#",Char2="Numbers",Char3="data"); 
  NcatchDataSeries <- as.numeric(DataFile[Index+2,1]);
  if(NcatchDataSeries>0) {  FixSigmaCatchN <- as.numeric(DataFile[Index+4,1:NcatchDataSeries]);
  } else {FixSigmaCatchN <- 0 }
  SigmaCatchNOffset <- as.numeric(DataFile[Index+6,1]);
  Index <- MatchTable(DataFile,Char1="#",Char2="The",Char3="numbers",Char4="data"); 
  Nnumbers  <- as.numeric(DataFile[Index+1,1]); Index <- Index + 2
  write(paste("Number of numbers points",Nnumbers),EchoFile,append=T)
  NumbersI <- matrix(0,nrow=Nnumbers,ncol=4)
  NumbersR <- matrix(0,nrow=Nnumbers,ncol=2)
  if(NcatchDataSeries>0){
    for (Inumber in 1:Nnumbers) 
     {
      NumbersI[Inumber,] <- as.numeric(DataFile[Index+Inumber,1:4])- c(0,1,0,1)
      NumbersI[Inumber,3] <- NumbersI[Inumber,3] - GeneralSpecs$Year1
      NumbersR[Inumber,] <- as.numeric(DataFile[Index+Inumber,5:6])
    }
    if (asnum(DataFile[Index+Inumber+1,6])) { print("Error reading Numbers data; too many inputs: Stopping"); AAA }
    
     }
  write("Numbers data",EchoFile,append=T)
  write(t(cbind(NumbersI,NumbersR)),EchoFile,append=T,ncol=6)
  
  # Size-composition
  Index <- MatchTable(DataFile,Char1="#",Char2="Length",Char3="compostion"); NlenComp  <- as.numeric(DataFile[Index+1,1]); Index <- Index + 2
  write("\nLength data",EchoFile,append=T)
  write(paste("Number of lines of composition data",NlenComp),EchoFile,append=T)
  LenCompI <- matrix(0,nrow=NlenComp,ncol=4)
  Stage1W <- rep(0,length=NlenComp)
  LenCompR <- matrix(0,nrow=NlenComp,ncol=GeneralSpecs$MaxLen)
  for (IlenC in 1:NlenComp)
   {
    LenCompI[IlenC,] <- as.numeric(DataFile[Index+IlenC,1:4]) - c(1,1,0,1)  ## This allows for R indexing in the dat file
    LenCompI[IlenC,3] <- LenCompI[IlenC,3] - GeneralSpecs$Year1
    Stage1W[IlenC] <- as.numeric(DataFile[Index+IlenC,5])
    NlenC <- GeneralSpecs$Nlen[LenCompI[IlenC,2]+1]
    LenCompR[IlenC,1:NlenC] <- as.numeric(DataFile[Index+IlenC,6:(5+NlenC)])
    LenCompR[IlenC,] <- LenCompR[IlenC,]/sum(LenCompR[IlenC,])
  }
  if (!is.na(DataFile[Index+IlenC,(6+NlenC)])) { print("Error reading Length data; too many inputs: Stopping"); AAA }
  write("Size-composition dats",EchoFile,append=T)
  write(t(cbind(LenCompI,LenCompR)),EchoFile,append=T,ncol=GeneralSpecs$MaxLen+4)

  # Read in the larval data
  Index <- MatchTable(DataFile,Char1="#",Char2="Larval",Char3="index")
  LarvalLikeOpt <- as.numeric(DataFile[Index+2,1]);                        # 0 for log-normal; otherwise normal
  Larval_Offset <- as.numeric(DataFile[Index+4,1]);                        # Time-delay between settlement and recruitment
  NLarvalData <- as.numeric(DataFile[Index+6,1]);
  Lar_dataI <- matrix(0,nrow=NLarvalData,ncol=2);
  Lar_dataR <- matrix(0,nrow=NLarvalData,ncol=2);
  Index <- Index + 7
  if (NLarvalData > 0){
    for (Idata in 1:NLarvalData)
    {
     for (II in 1:2) Lar_dataI[Idata,II] <- as.numeric(DataFile[Index+Idata,II]) - c(1,0)[II]
     for (II in 1:2) Lar_dataR[Idata,II] <- as.numeric(DataFile[Index+Idata,II+2])
     Lar_dataI[Idata,2] <- Lar_dataI[Idata,2]  - GeneralSpecs$Year1 + max(GeneralSpecs$BurnIn)
    }
    if (asnum(DataFile[Index+Idata+1,II+2])) { print("Error reading Larval data; too many inputs: Stopping"); AAA }
  }
  # Read in the environmental data
  Index <- MatchTable(DataFile,Char1="#",Char2="Environmental",Char3="Data")
  NenvSeries  <- as.numeric(DataFile[Index+2,1]); 
  write("\nEnvironmental data",EchoFile,append=T)
  write(paste("Number of environmental series",NenvSeries),EchoFile,append=T)
  YearsPerSeries <- as.numeric(DataFile[Index+4,1:NenvSeries])
  EnvData <- array(0,dim=c(GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1,GeneralSpecs$Nstep,NenvSeries))
  Ipnt <- Index+5
  for (Iseries in 1:NenvSeries)
   for (Iyr in 1:YearsPerSeries[Iseries])  
    {
     Ipnt <- Ipnt + 1
     Iyear <- as.numeric(DataFile[Ipnt,1])-GeneralSpecs$Year1+max(GeneralSpecs$BurnIn)+1;
     Istep <- as.numeric(DataFile[Ipnt,2]);
     Ienv <- as.numeric(DataFile[Ipnt,3]);
     EnvData[Iyear,Istep,Iseries]<- Ienv
   } 
  for (Istep in 1:GeneralSpecs$Nstep)
   {
    EnvData2 <- cbind(c(GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)+0:(GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn))),rep(Istep-1,1+GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)),EnvData[,Istep,])
    write(t(EnvData2),ncol=NenvSeries+2,EchoFile,append=T)
   } 
  if (asnum(DataFile[Ipnt+1,3])) { print("Error reading Enviromental data; too many inputs: Stopping"); AAA }
  
  print("READ IN THE DATA FILE")
  write("READ IN THE DATA FILE\n\n",EchoFile,append=T)
  
  ReturnObj <- NULL
  ReturnObj$Catch <- Catch
  ReturnObj$Ncpue <- Ncpue
  ReturnObj$NcpueDataSeries <- NcpueDataSeries
  ReturnObj$IndexType <- IndexType
  ReturnObj$FixSigmaCpue <- FixSigmaCpue
  ReturnObj$SigmaCpueOffset <- SigmaCpueOffset
  ReturnObj$TreatQcpue <- TreatQcpue
  ReturnObj$EnvIndCpue <- EnvIndCpue
  ReturnObj$EffCrIndCpue <- EffCrIndCpue
  ReturnObj$EffCrLag <- EffCrLag
  ReturnObj$IndexI <- IndexI
  ReturnObj$IndexR <- IndexR
  ReturnObj$LarvalLikeOpt <- LarvalLikeOpt
  ReturnObj$Larval_Offset <- Larval_Offset
  ReturnObj$NLarvalData <- NLarvalData
  ReturnObj$Lar_dataI <- Lar_dataI
  ReturnObj$Lar_dataR <- Lar_dataR
  ReturnObj$NenvSeries <- NenvSeries
  ReturnObj$EnvData <- EnvData
  ReturnObj$NQparPass = NQpars
  
  ReturnObj$Nnumbers <- Nnumbers
  ReturnObj$NcatchDataSeries <- NcatchDataSeries
  ReturnObj$FixSigmaCatchN <- FixSigmaCatchN
  ReturnObj$SigmaCatchNOffset <-SigmaCatchNOffset
  ReturnObj$NumbersI <- NumbersI
  ReturnObj$NumbersR <- NumbersR
  ReturnObj$NlenComp <- NlenComp
  ReturnObj$LenCompI <- LenCompI
  ReturnObj$Stage1W <- Stage1W
  ReturnObj$LenCompR <- LenCompR
  return(ReturnObj)
  
}

# ===================================================================================

ReadControlFile <- function(ControlFile,GeneralSpecs,DataSpecs)
{
  write("READING IN THE CONTROL FILE",EchoFile,append=T)
  
  # Weight-length regression
  Index <- MatchTable(ControlFile,Char1="#",Char2="weight-at-length",Char3=NULL)+1;
  WeightLen <- matrix(0,2,GeneralSpecs$MaxLen)
  for (Isex in 1:2) 
  {
    for (Jlen in 1:GeneralSpecs$Nlen[Isex]) WeightLen[Isex,Jlen] <- as.numeric(ControlFile[Index,Jlen])
    if (!is.na(ControlFile[Index,GeneralSpecs$Nlen[Isex]+1])) { print("Error reading weight-at-length; too many inputs: Stopping"); AAA }
    Index <- Index + 1
  }    
  write("Weight-length regressions",EchoFile,append=T)
  write(t(WeightLen),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Maturity",Char3=NULL)+1;
  Maturity <- matrix(0,GeneralSpecs$Narea,GeneralSpecs$MaxLen)
  for (Iarea in 1:GeneralSpecs$Narea) 
  {
    for (Jlen in 1:GeneralSpecs$Nlen[2]) Maturity[Iarea,Jlen] <- as.numeric(ControlFile[Index,Jlen])
    if (!is.na(ControlFile[Index,GeneralSpecs$Nlen[2]+1])) { print("Error reading maturity-at-length; too many inputs: Stopping"); AAA }
    Index <- Index + 1
  }    
  write("Maturity",EchoFile,append=T)
  write(t(Maturity),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)

  Index <- MatchTable(ControlFile,Char1="#",Char2="Egg",Char3="production")+1;
  MatFem <- matrix(0,GeneralSpecs$Narea,GeneralSpecs$MaxLen)
  for (Iarea in 1:GeneralSpecs$Narea) 
  {
    for (Jlen in 1:GeneralSpecs$Nlen[2]) MatFem[Iarea,Jlen] <- as.numeric(ControlFile[Index,Jlen])
    if (!is.na(ControlFile[Index,GeneralSpecs$Nlen[2]+1])) { print("Error reading eggs-at-length; too many inputs: Stopping"); AAA }
    Index <- Index + 1
  }    
  write("Egg Production",EchoFile,append=T)
  write(t(MatFem),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Egg",Char3="time")+1;
  MatTimeStep <- as.numeric(ControlFile[Index,1])
  
  # Link between fleets and areas
  Index <- MatchTable(ControlFile,Char1="#",Char2="Fleet",Char3="Area");
  Fleet_area <- rep(0,GeneralSpecs$Nfleet)
  for (Ifleet in 1:GeneralSpecs$Nfleet) 
    Fleet_area[Ifleet] <- as.numeric(ControlFile[Index+Ifleet,2])
  if (asnum(ControlFile[Index+Ifleet+1,2])) { print("Error reading Fleet data; too many inputs: Stopping"); AAA }
  write("Area for each fleet",EchoFile,append=T)
  write(Fleet_area,EchoFile,append=T)
  
  Area_fleet<- matrix(0,nrow=GeneralSpecs$Narea,ncol=GeneralSpecs$Nfleet)
  Narea_fleet <- rep(0,GeneralSpecs$Narea)
  for (Iarea in 1:GeneralSpecs$Narea)
  {  
    for (Ifleet in 1:GeneralSpecs$Nfleet)  
    {
      if (Fleet_area[Ifleet] == Iarea-1)
        Area_fleet[Iarea,Ifleet] = 1
      else
        Area_fleet[Iarea,Ifleet] = 0
      Narea_fleet[Iarea] = Narea_fleet[Iarea] + Area_fleet[Iarea,Ifleet]
    }
  }  
  write("Number of areas and fleets",EchoFile,append=T)
  write(Area_fleet,EchoFile,append=T)
  write("Number of fleets by area",EchoFile,append=T)
  write(Narea_fleet,EchoFile,append=T)
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Number",Char4="Zones");
  Nzone <- as.numeric(ControlFile[Index+1,1])
  Index <- MatchTable(ControlFile,Char1="#",Char2="Areas",Char4="each",Char5="Zone");
  NareasPerZone <- as.numeric(ControlFile[Index+1,1:Nzone])
  AreasPerZone <- matrix(-1,nrow=Nzone,ncol=GeneralSpecs$Narea)
  Index <- MatchTable(ControlFile,Char1="#",Char2="The",Char3="Zones");
  for (Izone in 1:Nzone) {Index <- Index+1
    AreasPerZone[Izone,1:NareasPerZone[Izone]] <- as.numeric(ControlFile[Index,1:NareasPerZone[Izone]])
    if (asnum(ControlFile[Index,NareasPerZone[Izone]+1]) | is.na(ControlFile[Index,NareasPerZone[Izone]])) { print("Error reading Number of areas / Zone data; not the correct number of inputs: Stopping"); AAA }
  }

  write("Number of zones",EchoFile,append=T)
  write(Nzone,EchoFile,append=T)
  write( NareasPerZone,EchoFile,append=T)
  write("Links between zones and area",EchoFile,append=T)
  write(t(AreasPerZone),EchoFile,ncol=GeneralSpecs$Narea,append=T)
  
  write("discard mortality",EchoFile,append=T)
  Phi <- array(0,dim=c(GeneralSpecs$Nfleet,GeneralSpecs$Nage,GeneralSpecs$Nyear+GeneralSpecs$MaxProjYr,GeneralSpecs$Nstep))
  Index <- MatchTable(ControlFile,Char1="#",Char2="Discard",Char3="mortality")+1;
  Ipnt <- 0
  for (Iage in 1:(GeneralSpecs$Nage))   
    for (Ifleet in 1:GeneralSpecs$Nfleet) 
      for (Istep in 1:GeneralSpecs$Nstep)  
      {
        Ipnt <- Ipnt + 1
        pos <- Index+Ipnt
        ControlFile[498,]
        testindex(ControlFile,pos,1,Iage) 
        testindex(ControlFile,pos,2,Ifleet) 
        testindex(ControlFile,pos,3,Istep) 
        for (Iyear in 1:GeneralSpecs$Nyear) Phi[Ifleet,Iage,Iyear,Istep] <- as.numeric(ControlFile[Index+Ipnt,3+Iyear])
        if (!is.na(ControlFile[Index+Ipnt,3+Iyear+1]) | is.na(ControlFile[Index+Ipnt,3+Iyear])) { print("Error reading Discard Mortality data; not the correct number of inputs: Stopping"); AAA }
      }
  
  Nout <- GeneralSpecs$Nyear*GeneralSpecs$Nfleet*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=3+GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Iage in 1:(GeneralSpecs$Nage))
      for (Iyear in 1:GeneralSpecs$Nyear)
      {
        Ipnt <- Ipnt + 1
        OutM[Ipnt,1:3]  <- c(Ifleet-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
        OutM[Ipnt,(4:(3+GeneralSpecs$Nstep))] <- Phi[Ifleet,Iage,Iyear,]
      }
  write(t(OutM),EchoFile,append=T,ncol=3+GeneralSpecs$Nstep)
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Recruitment_deviations"); 
  RecYr1 <- as.numeric(ControlFile[Index+1,1]) - GeneralSpecs$Year1+max(GeneralSpecs$BurnIn)
  RecYr2 <- GeneralSpecs$Year2+GeneralSpecs$MaxProjYr - GeneralSpecs$Year1+max(GeneralSpecs$BurnIn)
  write(paste("First recruitment year",RecYr1),EchoFile,append=T)
  write(paste("Last recruitment year",RecYr2),EchoFile,append=T)
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Spatial_deviations_in_recruitment"); 
  RecSpatYr1 <- as.numeric(ControlFile[Index+1,1]) - GeneralSpecs$Year1+max(GeneralSpecs$BurnIn)
  RecSpatYr2 <- GeneralSpecs$Year2+GeneralSpecs$MaxProjYr  - GeneralSpecs$Year1+max(GeneralSpecs$BurnIn)
  write(paste("First spatial recruitment year",RecSpatYr1),EchoFile,append=T)
  write(paste("Last spatial recruitment year",RecSpatYr2),EchoFile,append=T)
  
  # Initial conditions
  Index <- MatchTable(ControlFile,Char1="#",Char2="Initial_dev_option"); 
  InitOpt <- as.numeric(ControlFile[Index+1,1])
  InitParSpec <- as.numeric(ControlFile[Index+2,1])
  
  # read the data weights
  Index <- MatchTable(ControlFile,Char1="#",Char2="Weights",Char3="on"); 
  LambdaCpue <- as.numeric(ControlFile[Index+1,1])
  LambdaNumbers <- as.numeric(ControlFile[Index+2,1])
  LambdaLength <- as.numeric(ControlFile[Index+3,1])
  LambdaLarval <- as.numeric(ControlFile[Index+4,1])
  LambdaTag1 <- as.numeric(ControlFile[Index+5,1])
  LambdaTag2 <- as.numeric(ControlFile[Index+6,1])
  WeightInitialN <- as.numeric(ControlFile[Index+7,1])
  WeightInit3 <- as.numeric(ControlFile[Index+8,1])
  
  LambdaCpue2 <- rep(1.0,Data$NcpueDataSeries)
  LambdaNumbers2 <- rep(1.0,Data$NcatchDataSeries)
  LambdaLength2 <- array(1.0,dim=c(GeneralSpecs$Nfleet,GeneralSpecs$Nstep,GeneralSpecs$Nsex))
  Index <- MatchTable(ControlFile,Char1="#",Char2="Weights",Char3="by"); 
  UseLambdaDetailed <- as.numeric(ControlFile[Index+2,1])
  if (UseLambdaDetailed > 0)
  {
    for (II in 1:UseLambdaDetailed)
    {
      TheSpec <-as.numeric(ControlFile[Index+2+II,1:5])  
      if (TheSpec[1]==1)
      {
        IdataSet <- TheSpec[2]+1; Wght <- TheSpec[5]
        LambdaCpue2[IdataSet] <- Wght
      }  
      if (TheSpec[1]==2)
      {
        IdataSet <- TheSpec[2]+1; Wght <- TheSpec[5]
        LambdaNumbers2[IdataSet] <- Wght
      }  
      if (TheSpec[1]==3)
      {
        Ifleet <- TheSpec[2]+1;Istep <- TheSpec[3]; Isex <- TheSpec[4]+1; Wght <- TheSpec[5]
        if (Istep == -1)
          LambdaLength2[Ifleet,,Isex] <- Wght
        else
          LambdaLength2[Ifleet,Istep+1,Isex] <- Wght
      }  
      
    }  
  }
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Number",Char4="variance"); 
  NvarTypes <- as.numeric(ControlFile[Index+1,1])
  if (NvarTypes > 0) VarTypes <- as.numeric(ControlFile[Index+3,1:NvarTypes])
  if (NvarTypes <= 0) VarTypes <- 0
  
  Index <- MatchTable(ControlFile,Char1="#",Char2="Efficiency",Char3="parameters"); 
  NefficPar <- as.numeric(ControlFile[Index+1,1])

  print("READ IN THE CONTROL FILE")
  write("READ IN THE CONTROL FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$Fleet_area <- Fleet_area
  ReturnObj$Narea_fleet <- Narea_fleet
  ReturnObj$Area_fleet <- Area_fleet
  ReturnObj$Nzone <- Nzone
  ReturnObj$NareasPerZone <- NareasPerZone
  ReturnObj$AreasPerZone <- AreasPerZone
  ReturnObj$NefficPar <- NefficPar
  ReturnObj$Phi1 <- Phi
  ReturnObj$WeightLen <- WeightLen
  ReturnObj$MatFem <- MatFem
  ReturnObj$MatTimeStep <- MatTimeStep
  ReturnObj$RecYr1 <- RecYr1
  ReturnObj$RecYr2 <- RecYr2
  ReturnObj$RecSpatYr1 <- RecSpatYr1 
  ReturnObj$RecSpatYr2 <- RecSpatYr2
  ReturnObj$InitOpt <- InitOpt
  ReturnObj$InitParSpec <- InitParSpec
  ReturnObj$LambdaCpue <- LambdaCpue
  ReturnObj$LambdaNumbers <- LambdaNumbers
  ReturnObj$LambdaLength <- LambdaLength
  ReturnObj$LambdaLarval <- LambdaLarval
  ReturnObj$LambdaTag1 <- LambdaTag1
  ReturnObj$LambdaTag2 <- LambdaTag2
  ReturnObj$LambdaCpue2 <- LambdaCpue2
  ReturnObj$LambdaNumbers2 <- LambdaNumbers2
  ReturnObj$LambdaLength2 <- LambdaLength2
  ReturnObj$WeightInitialN <- WeightInitialN
  ReturnObj$WeightInit3 <- WeightInit3
  ReturnObj$NvarTypes <- NvarTypes
  ReturnObj$VarTypes <- VarTypes
  return(ReturnObj)
  
}
# ===================================================================================

ReadMoveFile <- function(MoveFile,GeneralSpecs)
{
  Index <- MatchTable(MoveFile,Char1="#",Char2="Number",Char3="of",Char4="movement");
  NmovePatterns <- as.numeric(MoveFile[Index+1,1])
  write(paste("Number of movement patterns",NmovePatterns),EchoFile,append=T)
  MoveSpec<-matrix(0,NmovePatterns,4)
  for (Ipatt in 1:NmovePatterns){
    for (II in 1:4) {
      MoveSpec[Ipatt,II] <- as.numeric(MoveFile[Index+2+Ipatt,II])
      
    }
    }
  write("Specifications for movement",EchoFile,append=T)
  write(MoveSpec,EchoFile,append=T,ncol=4)
    
  Index <- Index + 2+NmovePatterns+1
  NmovePars = 0
  for (Ipatt in 1:NmovePatterns)
   {
    if (MoveSpec[Ipatt,2]==1) NmovePars = NmovePars + 1
    if (MoveSpec[Ipatt,2]==2) NmovePars = NmovePars + 2
    }  
  write(paste("Number of movement parameters",NmovePars),EchoFile,append=T)
  
  Index <- MatchTable(MoveFile,Char1="#",Char2="Movement",Char3="specifications")+2;
  MovePnt <- array(0,dim=c(GeneralSpecs$Narea,GeneralSpecs$Nage,GeneralSpecs$Nyear,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Iage in 1:(GeneralSpecs$Nage))
   for (Iarea in 1:GeneralSpecs$Narea) 
    for (Istep in 1:GeneralSpecs$Nstep)  
     {
      pos <- Index+Ipnt
      testindex(MoveFile,pos,1,Iage) 
      testindex(MoveFile,pos,2,Iarea) 
      testindex(MoveFile,pos,3,Istep) 
      for (Iyear in 1:GeneralSpecs$Nyear) MovePnt[Iarea,Iage,Iyear,Istep] <- as.numeric(MoveFile[Index+Ipnt,3+Iyear])
      Ipnt <- Ipnt + 1
    }
  write("Specifications for movement pointers",EchoFile,append=T)
  Nout <- GeneralSpecs$Nyear*GeneralSpecs$Narea*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=3+GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Iage in 1:(GeneralSpecs$Nage))
   for (Iarea in 1:GeneralSpecs$Narea)
    for (Iyear in 1:GeneralSpecs$Nyear)
      {
        Ipnt <- Ipnt + 1
        OutM[Ipnt,1:3]  <- c(Iarea-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
        OutM[Ipnt,(4:(3+GeneralSpecs$Nstep))] <- MovePnt[Iarea,Iage,Iyear,]
      }
  write(t(OutM),EchoFile,append=T,ncol=3+GeneralSpecs$Nstep)
  
    
  print("READ IN THE MOVEMENT FILE")
  write("READ IN THE MOVEMENT FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NmovePatterns <- NmovePatterns
  ReturnObj$MoveSpec <- MoveSpec
  ReturnObj$MovePnt <- MovePnt
  ReturnObj$NmovePars <- NmovePars
  return(ReturnObj)
 }


# ===================================================================================

ReadSelexFile <- function(SelexFile,GeneralSpecs)
{
  Index <- MatchTable(SelexFile,Char1="#",Char2="Number",Char3="Selex")+1;
  NselPatterns <- as.numeric(SelexFile[Index,1])
  write(paste("Number of selectivity patterns",NselPatterns),EchoFile,append=T)
  
  Index <- MatchTable(SelexFile,Char1="#",Char2="Pattern",Char3="Type")[1];
  SelSpec <- matrix(0,NselPatterns,5)
  for (Isel in 1:NselPatterns)
    for (Icol in 1:5) SelSpec[Isel,Icol] <- as.numeric(SelexFile[Index+Isel,Icol]) 
  write("Specifications for selectivity",EchoFile,append=T)
  write(t(SelSpec),EchoFile,append=T,ncol=5)

  NselPars <- 0
  for (Isel in 1:NselPatterns) 
  {
    if (SelSpec[Isel,2]== 2) NselPars <- NselPars + SelSpec[Isel,4]
    if (SelSpec[Isel,2]== 3) NselPars <- NselPars + 2
  }
  write(paste("Number of selectivity parameters",NselPars),EchoFile,append=T)
  
  Index <- MatchTable(SelexFile,Char1="#",Char2="Specifications",Char3="for",Char4="selectivity")+2;
  SelPnt <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nfleet,GeneralSpecs$Nyear,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Ifleet in 1:GeneralSpecs$Nfleet) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          pos <- Index+Ipnt
          testindex(SelexFile,pos,1,Isex) 
          testindex(SelexFile,pos,2,Iage) 
          testindex(SelexFile,pos,3,Ifleet) 
          testindex(SelexFile,pos,4,Istep) 
          for (Iyear in 1:GeneralSpecs$Nyear) SelPnt[Isex,Iage,Ifleet,Iyear,Istep] <- as.numeric(SelexFile[pos,4+Iyear])
          Ipnt <- Ipnt + 1
        } 
  write("Specifications for selectivity pointers",EchoFile,append=T)
  Nout <- GeneralSpecs$Nyear*GeneralSpecs$Nfleet*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Ifleet in 1:GeneralSpecs$Nfleet)
   for (Isex in 1:(GeneralSpecs$Nsex))
    for (Iage in 1:(GeneralSpecs$Nage))
     for (Iyear in 1:GeneralSpecs$Nyear)
      {
       Ipnt <- Ipnt + 1
       OutM[Ipnt,1:4]  <- c(Ifleet-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
       OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- SelPnt[Isex,Iage,Ifleet,Iyear,]
      }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)

  # Fixed selectivity
  Index <- MatchTable(SelexFile,Char1="#",Char2="selectivity",Char3=NULL)+1;
  NfixedSelex <- as.numeric(SelexFile[Index,1])
  write(paste("\nNumber of fixed selectivity patterns",NfixedSelex),EchoFile,append=T)
  SelexFI <- matrix(0,NfixedSelex,GeneralSpecs$MaxLen)
  for (Ifleet in 1:NfixedSelex) 
    for (Jlen in 1:GeneralSpecs$Nlen[1])SelexFI[Ifleet,Jlen] <- as.numeric(SelexFile[Index+Ifleet,Jlen])
  write("Fixed selectivity patterns",EchoFile,append=T)
  write(t(SelexFI),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)

  # Legal patterns
  Index <- MatchTable(SelexFile,Char1="#",Char2="Number",Char3="Legal")+1;
  NlegalPatterns <- as.numeric(SelexFile[Index,1])
  write(paste("Number of Legal patterns",NlegalPatterns),EchoFile,append=T)
  
  Index <- MatchTable(SelexFile,Char1="#",Char2="Pattern",Char3="Type")[2];
  LegalSpec <- matrix(0,NlegalPatterns,4)
  for (Isel in 1:NlegalPatterns)
    for (Icol in 1:4) LegalSpec[Isel,Icol] <- as.numeric(SelexFile[Index+Isel,Icol]) 
  write("Specifications for legal",EchoFile,append=T)
  write(t(LegalSpec),EchoFile,append=T,ncol=4)
  
  # Fixed legal
  Index <- MatchTable(SelexFile,Char1="#",Char2="legal",Char3=NULL)+1;
  NfixedLegal <- as.numeric(SelexFile[Index,1])
  write(paste("\nNumber of fixed legal patterns",NfixedLegal),EchoFile,append=T)
  LegalFI <- matrix(0,NfixedLegal,GeneralSpecs$MaxLen)
  for (Ifleet in 1:NfixedLegal) 
    for (Jlen in 1:GeneralSpecs$Nlen[1]) LegalFI[Ifleet,Jlen] <- as.numeric(SelexFile[Index+Ifleet,Jlen])
  write("Fixed legal patterns",EchoFile,append=T)
  write(t(LegalFI),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)

  Index <- MatchTable(SelexFile,Char1="#",Char2="Specifications",Char3="for",Char4="Fleet")+2;
  LegalFleetPnt <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nfleet,GeneralSpecs$Nyear,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Ifleet in 1:GeneralSpecs$Nfleet) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          pos <- Index+Ipnt
          testindex(SelexFile,pos,1,Isex) 
          testindex(SelexFile,pos,2,Iage) 
          testindex(SelexFile,pos,3,Ifleet) 
          testindex(SelexFile,pos,4,Istep) 
          for (Iyear in 1:GeneralSpecs$Nyear) LegalFleetPnt[Isex,Iage,Ifleet,Iyear,Istep] <- as.numeric(SelexFile[Index+Ipnt,4+Iyear])
          Ipnt <- Ipnt + 1
        } 
  write("Specifications for fleet legal pointers",EchoFile,append=T)
  Nout <- GeneralSpecs$Nyear*GeneralSpecs$Nfleet*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Isex in 1:(GeneralSpecs$Nsex))
      for (Iage in 1:(GeneralSpecs$Nage))
        for (Iyear in 1:GeneralSpecs$Nyear)
        {
          Ipnt <- Ipnt + 1
          OutM[Ipnt,1:4]  <- c(Ifleet-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
          OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- LegalFleetPnt[Isex,Iage,Ifleet,Iyear,]
        }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)

  Index <- MatchTable(SelexFile,Char1="#",Char2="Specifications",Char3="for",Char4="legal")+2;
  LegalPnt <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Narea,GeneralSpecs$Nyear,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Iarea in 1:GeneralSpecs$Narea) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
         pos <- Index+Ipnt
         testindex(SelexFile,pos,1,Isex) 
         testindex(SelexFile,pos,2,Iage) 
         testindex(SelexFile,pos,3,Iarea) 
         testindex(SelexFile,pos,4,Istep) 
         for (Iyear in 1:GeneralSpecs$Nyear) LegalPnt[Isex,Iage,Iarea,Iyear,Istep] <- as.numeric(SelexFile[Index+Ipnt,4+Iyear])
         Ipnt <- Ipnt + 1
        } 
  write("Specifications for legal pointers",EchoFile,append=T)
  Nout <- GeneralSpecs$Nyear*GeneralSpecs$Narea*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:(GeneralSpecs$Nsex))
      for (Iage in 1:(GeneralSpecs$Nage))
        for (Iyear in 1:GeneralSpecs$Nyear)
         {
          Ipnt <- Ipnt + 1
          OutM[Ipnt,1:4]  <- c(Iarea-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
          OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- LegalPnt[Isex,Iage,Iarea,Iyear,]
         }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)
  
  Index <- MatchTable(SelexFile,Char1="#",Char2="Reference",Char3="selectivity",Char4="pattern");
  LegalRef <- rep(0,GeneralSpecs$MaxLen)
  for (Jlen in 1:GeneralSpecs$Nlen[1]) LegalRef[Jlen] <- as.numeric(SelexFile[Index+1,Jlen])
  write("Reference legal pattern",EchoFile,append=T)
  write(LegalRef,EchoFile,append=T,ncol=GeneralSpecs$MaxLen)
  
  Index <- MatchTable(SelexFile,Char1="#",Char2="IsRed",Char3="specifications",Char4="-")+2;
  IsRed <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Narea,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Iarea in 1:GeneralSpecs$Narea) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          pos <- Index+Ipnt
          testindex(SelexFile,pos,1,Isex) 
          testindex(SelexFile,pos,2,Iage) 
          testindex(SelexFile,pos,3,Iarea) 
          testindex(SelexFile,pos,4,Istep) 
          IsRed[Isex,Iage,Iarea,Istep] <- as.numeric(SelexFile[Index+Ipnt,5])
          Ipnt <- Ipnt + 1
        } 
  write("Specifications for IsRed",EchoFile,append=T)
  Nout <- GeneralSpecs$Narea*GeneralSpecs$Nsex*(GeneralSpecs$Nage)*(GeneralSpecs$Nstep)
  OutM <- matrix(0,nrow=Nout,ncol=5)
  Ipnt <- 0
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:(GeneralSpecs$Nsex))
      for (Iage in 1:(GeneralSpecs$Nage))
        for (Istep in 1:GeneralSpecs$Nstep)
        {
          Ipnt <- Ipnt + 1
          OutM[Ipnt,1:4]  <- c(Iarea-1,Isex-1,Iage-1,Istep-1)
          OutM[Ipnt,5] <- IsRed[Isex,Iage,Iarea,Istep]
        }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)

  print("READ IN THE SELEX FILE")
  write("READ IN THE SELEX FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NselPatterns <- NselPatterns
  ReturnObj$NlegalPatterns <- NlegalPatterns
  ReturnObj$SelSpec <- SelSpec
  ReturnObj$SelPnt <- SelPnt
  ReturnObj$SelexFI <- SelexFI
  ReturnObj$NselPars <- NselPars
  ReturnObj$NfixedSelex <- NfixedSelex
  ReturnObj$LegalSpec <- LegalSpec
  ReturnObj$LegalPnt <- LegalPnt
  ReturnObj$LegalFleetPnt <- LegalFleetPnt     
  ReturnObj$LegalFI <- LegalFI
  ReturnObj$LegalRef <- LegalRef
  ReturnObj$NlegalSelex <- NfixedLegal
  ReturnObj$IsRed <- IsRed
  return(ReturnObj)
}

# ===================================================================================

ReadRetenFile <- function(RetenFile,GeneralSpecs)
{
  Index <- MatchTable(RetenFile,Char1="#",Char2="Number",Char3="Retain")+1;
  NretPatterns <- as.numeric(RetenFile[Index,1])
  write(paste("Number of retention patterns",NretPatterns),EchoFile,append=T)

    Index <- MatchTable(RetenFile,Char1="#",Char2="Pattern",Char3="Type");
  RetSpec <- matrix(0,NretPatterns,5)
  for (Iret in 1:NretPatterns)
    for (Icol in 1:5) RetSpec[Iret,Icol] <- as.numeric(RetenFile[Index+Iret,Icol]) 
  write("Specifications for retention",EchoFile,append=T)
  write(RetSpec,EchoFile,append=T,ncol=5)
  NretPars <- 0
  for (Iret in 1:NretPatterns) 
  {
    if (RetSpec[Iret,2]== 2) NretPars <- NretPars + RetSpec[Iret,4]
    if (RetSpec[Iret,2]== 3) NretPars <- NretPars + 2
    if (RetSpec[Iret,2]== 4) NretPars <- NretPars + 1
    if (RetSpec[Iret,2]== 6) NretPars <- NretPars + 3
    if (RetSpec[Iret,2]== 7) NretPars <- NretPars + 2
    if (RetSpec[Iret,2]== 8) NretPars <- NretPars + 1
  }
  write(paste("Number of retention parameters",NretPars),EchoFile,append=T)
  
  Index <- MatchTable(RetenFile,Char1="#",Char2="Specifications",Char3="for")+2;
  RetPnt <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nfleet,GeneralSpecs$Nyear,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Ifleet in 1:GeneralSpecs$Nfleet) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          pos <- Index+Ipnt
          testindex(RetenFile,pos,1,Isex) 
          testindex(RetenFile,pos,2,Iage) 
          testindex(RetenFile,pos,3,Ifleet) 
          testindex(RetenFile,pos,4,Istep) 
          for (Iyear in 1:GeneralSpecs$Nyear) RetPnt[Isex,Iage,Ifleet,Iyear,Istep] <- as.numeric(RetenFile[Index+Ipnt,4+Iyear])
          Ipnt <- Ipnt + 1
        } 
  write("Specifications for retention pointers",EchoFile,append=T)
  Nout <- GeneralSpecs$Nyear*GeneralSpecs$Nfleet*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Isex in 1:(GeneralSpecs$Nsex))
      for (Iage in 1:(GeneralSpecs$Nage))
        for (Iyear in 1:GeneralSpecs$Nyear)
        {
          Ipnt <- Ipnt + 1
          OutM[Ipnt,1:4]  <- c(Ifleet-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
          OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- RetPnt[Isex,Iage,Ifleet,Iyear,]
        }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)

  # Fixed retention
  Index <- MatchTable(RetenFile,Char1="#",Char2="retention",Char3=NULL)+1;
  NfixedReten <- as.numeric(RetenFile[Index,1])
  write(paste("Number of fixed retention patterns",NfixedReten),EchoFile,append=T)
  RetenFI <- matrix(0,NfixedReten,GeneralSpecs$MaxLen)
  for (Ifleet in 1:NfixedReten) 
    for (Jlen in 1:GeneralSpecs$MaxLen) RetenFI[Ifleet,Jlen] <- as.numeric(RetenFile[Index+Ifleet,Jlen])
  write("Fixed retention patterns",EchoFile,append=T)
  write(RetenFI,EchoFile,append=T,ncol=GeneralSpecs$MaxLen)

  print("READ IN THE RETAIN FILE")
  write("READ IN THE RETAIN FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj <- NULL
  ReturnObj$NretPatterns <- NretPatterns
  ReturnObj$RetSpec <- RetSpec
  ReturnObj$RetPnt <- RetPnt
  ReturnObj$RetenFI <- RetenFI
  ReturnObj$NretPars <- NretPars
  ReturnObj$NfixedReten <- NfixedReten
  return(ReturnObj)
  
}

# ===================================================================================

ReadRecruitFile <- function(RecruitFile,GeneralSpecs)
 {
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Number",Char3="of",Char4="sex_area_allocation")+1;
  NrecruitPatternsA <- as.numeric(RecruitFile[Index,1])
  write(paste("Number of recruitment patterns",NrecruitPatternsA),EchoFile,append=T)

  RecruitSpecsA<-matrix (0,nrow=NrecruitPatternsA,ncol=3)                             
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Allocate_yearxarea")+1;
  for (Ipat in 1:NrecruitPatternsA)
    RecruitSpecsA[Ipat,]<- as.numeric(RecruitFile[Index+Ipat,1:3])
  write("Specifications for recuitment",EchoFile,append=T)
  write(t(RecruitSpecsA),EchoFile,append=T,ncol=2)

  Index <- MatchTable(RecruitFile,Char1="#",Char2="Number",Char3="of",Char4="length_allocation")+1;
  NrecruitPatternsB <- as.numeric(RecruitFile[Index,1])
  write(paste("Number of rectuitment patterns",NrecruitPatternsB),EchoFile,append=T)
  
  RecruitSpecsB<-matrix (0,nrow=NrecruitPatternsB,ncol=1+2*GeneralSpecs$Nsex)                             
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Allocate_length")+1;
  for (Ipat in 1:NrecruitPatternsB)
    RecruitSpecsB[Ipat,]<- as.numeric(RecruitFile[Index+Ipat,1:(1+2*GeneralSpecs$Nsex)])
  write("Specifications for recuitment",EchoFile,append=T)
  write(t(RecruitSpecsB),EchoFile,append=T,ncol=1+2*GeneralSpecs$Nsex)
  
  NrecruitPars = 0
  for (IrecPat in 1:NrecruitPatternsA)
   {
    if (RecruitSpecsA[IrecPat,2]==0) NrecruitPars <- NrecruitPars + (GeneralSpecs$Nsex+GeneralSpecs$Narea-2)
    if (RecruitSpecsA[IrecPat,2]==1) NrecruitPars <- NrecruitPars + (GeneralSpecs$Nsex*GeneralSpecs$Narea-1)
   } 
  for (Ipat in 1:NrecruitPatternsB)
    for (Isex in 1:GeneralSpecs$Nsex)
      if (RecruitSpecsB[Ipat,3+Isex] >=0) NrecruitPars <- NrecruitPars+ RecruitSpecs[Ipat,3+GeneralSpecs$Nsex+Isex]
  write(paste("Number of recuitment parameters",NrecruitPars),EchoFile,append=T,ncol=3+GeneralSpecs$Nsex)
  
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Recruit",Char3="by",Char4="year")+2;
  RecruitPnt <- matrix(0,nrow=GeneralSpecs$Nyear+GeneralSpecs$MaxProjYr,ncol=GeneralSpecs$Nstep)
  Ipnt <- 0
  for (Istep in 1:GeneralSpecs$Nstep)  
  {
    for (Iyear in 1:(GeneralSpecs$Nyear+GeneralSpecs$MaxProjYr))  RecruitPnt[Iyear,Istep] <- as.numeric(RecruitFile[Index+Ipnt,1+Iyear]) 
    Ipnt <- Ipnt + 1
  } 
  write("Specifications for recruitment pointers",EchoFile,append=T)
  write(t(RecruitPnt),EchoFile,append=T,ncol=GeneralSpecs$Nstep)
  
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Recruit",Char3="by",Char4="area")+1;
  RecruitLenPnt <- rep(0,length=GeneralSpecs$Narea)
  for (Iarea in 1:GeneralSpecs$Narea) RecruitLenPnt[Iarea] <- as.numeric(RecruitFile[Index,Iarea])   
  
  # Recruitment proportions
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Recruitment",Char3="fractions",Char4="use")+1;
  CalcRecruitFrac <- as.numeric(RecruitFile[Index,1]) 
  write(paste("Use pre-specified recruitment fractions",CalcRecruitFrac),EchoFile,append=T)
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Number",Char3="of",Char4="pre-specified")+1;
  NfixedRecruits <- as.numeric(RecruitFile[Index,1])
  write(paste("Number of fixed recruitment patterns",NfixedRecruits),EchoFile,append=T)
  # AEP 
  RecruitFrac <- matrix(1,NfixedRecruits,GeneralSpecs$MaxLen)  # make blank index
  #RecruitFrac <- matrix(1,GeneralSpecs$Nsex,GeneralSpecs$MaxLen)
  if(CalcRecruitFrac==0){
    for (Isex in 1:NfixedRecruits) 
      for (Jlen in 1:GeneralSpecs$MaxLen) RecruitFrac[Isex,Jlen] <- as.numeric(RecruitFile[Index+Isex+1,Jlen])
    }
   write("Fixed recruitment patterns",EchoFile,append=T)
   write(t(RecruitFrac),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)
    
   
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Bias",Char3="ramp")+1;
  Bias_Ramp_Yr1 <- as.numeric(RecruitFile[Index,1])-GeneralSpecs$Year1;
  Bias_Ramp_Yr2 <- as.numeric(RecruitFile[Index,2])-GeneralSpecs$Year1;
  Bias_Ramp_Yr3 <- as.numeric(RecruitFile[Index,3])-GeneralSpecs$Year1;
  Bias_Ramp_Yr4 <- as.numeric(RecruitFile[Index,4])-GeneralSpecs$Year1;
  

  print("READ IN THE RECRUIT FILE")
  write("READ IN THE RECRUIT FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NrecruitPatternsA <- NrecruitPatternsA
  ReturnObj$RecruitSpecsA <- RecruitSpecsA
  ReturnObj$NrecruitPatternsB <- NrecruitPatternsB
  ReturnObj$RecruitSpecsB <- RecruitSpecsB
  ReturnObj$RecruitPnt <- RecruitPnt
  ReturnObj$RecruitLenPnt <- RecruitLenPnt
  ReturnObj$RecruitFrac <- RecruitFrac
  ReturnObj$NrecruitPars <- NrecruitPars
  ReturnObj$NfixedRecruits <- NfixedRecruits
  ReturnObj$CalcRecruitFrac <- CalcRecruitFrac
  ReturnObj$Bias_Ramp_Yr1  <- Bias_Ramp_Yr1;
  ReturnObj$Bias_Ramp_Yr2  <- Bias_Ramp_Yr2;
  ReturnObj$Bias_Ramp_Yr3  <- Bias_Ramp_Yr3;
  ReturnObj$Bias_Ramp_Yr4  <- Bias_Ramp_Yr4;
  return(ReturnObj)
}

# ===================================================================================

ReadGrowthFile <- function(GrowthFile,GeneralSpecs)
{
  Index <- MatchTable(GrowthFile,Char1="#",Char2="Number",Char3="of",Char4="growth")+1;
  NgrowthPatterns <- as.numeric(GrowthFile[Index,1])
  write(paste("Number of growth patterns",NgrowthPatterns),EchoFile,append=T)
  
  GrowthSpecs<-matrix (0,nrow=NgrowthPatterns,ncol=6)
  for (Ipat in 1:NgrowthPatterns)
    GrowthSpecs[Ipat,]<- as.numeric(GrowthFile[Index+1+Ipat,1:6])
  write("Specifications for growth",EchoFile,append=T)
  write(t(GrowthSpecs),EchoFile,append=T,ncol=6)
  
  NgrowthPars <- 0
  write(paste("Number of growth parameters",NgrowthPars),EchoFile,append=T,ncol=3+GeneralSpecs$Nsex)
  
  Index <- MatchTable(GrowthFile,Char1="#",Char2="Specifications",Char3="for")+2;  #  Changed from 3 to 2 as there was some erroneous text in original growth file
  GrowthPnt <- array(0,dim=c(GeneralSpecs$Narea,GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nyear,GeneralSpecs$Nstep))
  Ipnt <- 0
  for (Isex in 1:GeneralSpecs$Nsex)
   for (Iage in 1:(GeneralSpecs$Nage))  
     for (Iarea in 1:GeneralSpecs$Narea)
       for (Istep in 1:GeneralSpecs$Nstep)  
      {
       pos <- Index+Ipnt
       testindex(GrowthFile,pos,1,Isex) 
       testindex(GrowthFile,pos,2,Iage) 
       testindex(GrowthFile,pos,3,Iarea) 
       testindex(GrowthFile,pos,4,Istep) 
       for (Iyear in 1:GeneralSpecs$Nyear) GrowthPnt[Iarea,Isex,Iage,Iyear,Istep] <- as.numeric(GrowthFile[Index+Ipnt,4+Iyear])
       Ipnt <- Ipnt + 1
      } 
  write("Specifications for growth pointers",EchoFile,append=T)
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:GeneralSpecs$Nsex)
      for (Iage in 1:(GeneralSpecs$Nage))  
        for (Iyear in 1:GeneralSpecs$Nyear) 
          write(GrowthPnt[Iarea,Isex,Iage,Iyear,],EchoFile,append=T,ncol=GeneralSpecs$Nstep)
  
  # Size-transition matrix
  Index <- MatchTable(GrowthFile,Char1="#",Char2="Number",Char3="prespecified")+1;
  NfixedGrowth = as.numeric(GrowthFile[Index,1])
  write(paste("Number of fixed growth patterns",NfixedGrowth),EchoFile,append=T)
  NfixedGrowthSex = as.numeric(GrowthFile[Index+2,1:NfixedGrowth])
  write("Number of fixed growth patterns by sex",EchoFile,append=T)
  write(NfixedGrowthSex,EchoFile,append=T,ncol=NfixedGrowth)
  
  TransInp <- array(0,dim=c(NfixedGrowth,GeneralSpecs$MaxLen,GeneralSpecs$MaxLen))
  for (Igrow in 1:NfixedGrowth) 
  {
    Isex <- NfixedGrowthSex[Igrow]
    Index <- Index + 1
    for (Ilen in 1:GeneralSpecs$Nlen[Isex+1]) 
    {  
      for (Jlen in 1:GeneralSpecs$Nlen[Isex+1]) TransInp[Igrow,Ilen,Jlen] <- as.numeric(GrowthFile[Index+4,Jlen])
      Index <- Index + 1
    }    
  }  
  for (Igrow in 1:NfixedGrowth) 
   write(TransInp[Igrow,,],EchoFile,append=T,ncol=GeneralSpecs$MaxLen)
  
  print("READ IN THE GROWTH FILE")
  write("READ IN THE GROWTH FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NgrowthPatterns <- NgrowthPatterns
  ReturnObj$GrowthSpecs <- GrowthSpecs
  ReturnObj$GrowthPnt <- GrowthPnt
  ReturnObj$TransInp <-TransInp
  ReturnObj$NfixedGrowth <- NfixedGrowth
  ReturnObj$NfixedGrowthSex <- NfixedGrowthSex
  ReturnObj$NgrowthPars <- NgrowthPars
  return(ReturnObj)
  
}

# ===================================================================================

ReadProjFile <- function(ProjFile,GeneralSpecs,Phi1)
{
  Index <- MatchTable(ProjFile,Char1="#",Char2="Number",Char3="of",Char4="projection")+1; 
  Nproj <- as.numeric(ProjFile[Index,1])
  print(Nproj)

  Index <- MatchTable(ProjFile,Char1="#",Char2="Specifications",Char3="for",Char4="selectivity")+2;
  SelPntFut <- array(-1,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nfleet,GeneralSpecs$MaxProjYr,GeneralSpecs$Nstep))
  Ipnt <- 0
  if (Nproj > 0)
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Ifleet in 1:GeneralSpecs$Nfleet) 
        for (Istep in 1:GeneralSpecs$Nstep)  
         {
          pos <- Index+Ipnt
          testindex(ProjFile,pos,1,Isex) 
          testindex(ProjFile,pos,2,Iage) 
          testindex(ProjFile,pos,3,Ifleet) 
          testindex(ProjFile,pos,4,Istep) 
          for (Iyear in 1:Nproj) SelPntFut[Isex,Iage,Ifleet,Iyear,Istep] <- as.numeric(ProjFile[pos,4+Iyear])
          Ipnt <- Ipnt + 1
         } 
  write("Specifications for selectivity pointers",EchoFile,append=T)
  Nout <- Nproj*GeneralSpecs$Nfleet*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  if (Nproj > 0)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
   for (Isex in 1:(GeneralSpecs$Nsex))
    for (Iage in 1:(GeneralSpecs$Nage))
     for (Iyear in 1:Nproj)
      {
       Ipnt <- Ipnt + 1
       OutM[Ipnt,1:4]  <- c(Ifleet-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Nyear+GeneralSpecs$Year1-1)
       OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- SelPntFut[Isex,Iage,Ifleet,Iyear,]
      }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)
    
  Index <- MatchTable(ProjFile,Char1="#",Char2="Specifications",Char3="for",Char4="retention.")+2;
  RetPntFut <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nfleet,GeneralSpecs$MaxProjYr,GeneralSpecs$Nstep))
  Ipnt <- 0
  if (Nproj > 0)
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Ifleet in 1:GeneralSpecs$Nfleet) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          pos <- Index+Ipnt
          testindex(ProjFile,pos,1,Isex) 
          testindex(ProjFile,pos,2,Iage) 
          testindex(ProjFile,pos,3,Ifleet) 
          testindex(ProjFile,pos,4,Istep) 
          for (Iyear in 1:Nproj) RetPntFut[Isex,Iage,Ifleet,Iyear,Istep] <- as.numeric(ProjFile[Index+Ipnt,4+Iyear])
          Ipnt <- Ipnt + 1
        } 
  write("Specifications for retention pointers",EchoFile,append=T)
  Nout <- Nproj*GeneralSpecs$Nfleet*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  if (Nproj > 0)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Isex in 1:(GeneralSpecs$Nsex))
      for (Iage in 1:(GeneralSpecs$Nage))
        for (Iyear in 1:Nproj)
        {
          Ipnt <- Ipnt + 1
          OutM[Ipnt,1:4]  <- c(Ifleet-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Nyear+GeneralSpecs$Year1-1)
          OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- RetPntFut[Isex,Iage,Ifleet,Iyear,]
        }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)

  Index <- MatchTable(ProjFile,Char1="#",Char2="Specifications",Char3="for",Char4="Fleet")+2;
  LegalFleetPntFut <- array(0,dim=c(GeneralSpecs$Nsex,GeneralSpecs$Nage,GeneralSpecs$Nfleet,GeneralSpecs$MaxProjYr,GeneralSpecs$Nstep))
  Ipnt <- 0
  if (Nproj > 0)
  for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))  
      for (Ifleet in 1:GeneralSpecs$Nfleet) 
        for (Istep in 1:GeneralSpecs$Nstep)  
        {
          pos <- Index+Ipnt
          testindex(ProjFile,pos,1,Isex) 
          testindex(ProjFile,pos,2,Iage) 
          testindex(ProjFile,pos,3,Ifleet) 
          testindex(ProjFile,pos,4,Istep) 
          for (Iyear in 1:Nproj) LegalFleetPntFut[Isex,Iage,Ifleet,Iyear,Istep] <- as.numeric(ProjFile[Index+Ipnt,4+Iyear])
          Ipnt <- Ipnt + 1
        } 
  write("Specifications for fleet legal pointers",EchoFile,append=T)
  Nout <- Nproj*GeneralSpecs$Nfleet*GeneralSpecs$Nsex*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=4+GeneralSpecs$Nstep)
  Ipnt <- 0
  if (Nproj > 0)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Isex in 1:(GeneralSpecs$Nsex))
      for (Iage in 1:(GeneralSpecs$Nage))
        for (Iyear in 1:Nproj)
        {
          Ipnt <- Ipnt + 1
          OutM[Ipnt,1:4]  <- c(Ifleet-1,Isex-1,Iage-1,Iyear+GeneralSpecs$Nyear++GeneralSpecs$Year1-1)
          OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- LegalFleetPntFut[Isex,Iage,Ifleet,Iyear,]
        }
  write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)

  write("discard mortality",EchoFile,append=T)
  Phi <- Phi1
  Index <- MatchTable(ProjFile,Char1="#",Char2="Discard",Char3="mortality")+1;
  Ipnt <- 0
  if (Nproj > 0)
  for (Iage in 1:(GeneralSpecs$Nage))   
    for (Ifleet in 1:GeneralSpecs$Nfleet) 
      for (Istep in 1:GeneralSpecs$Nstep)  
      {
        Ipnt <- Ipnt + 1
        pos <- Index+Ipnt
        ControlFile[498,]
        testindex(ProjFile,pos,1,Iage) 
        testindex(ProjFile,pos,2,Ifleet) 
        testindex(ProjFile,pos,3,Istep) 
        for (Iyear in 1:Nproj) Phi[Ifleet,Iage,Iyear+GeneralSpecs$Nyear,Istep] <- as.numeric(ProjFile[Index+Ipnt,3+Iyear])
      }

  Nout <- (Nproj+GeneralSpecs$Nyear)*GeneralSpecs$Nfleet*(GeneralSpecs$Nage)
  OutM <- matrix(0,nrow=Nout,ncol=3+GeneralSpecs$Nstep)
  Ipnt <- 0
  if (Nproj > 0)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
    for (Iage in 1:(GeneralSpecs$Nage))
      for (Iyear in 1:(GeneralSpecs$Nyear+Nproj))
      {
        Ipnt <- Ipnt + 1
        OutM[Ipnt,1:3]  <- c(Ifleet-1,Iage-1,Iyear+GeneralSpecs$Year1-1)
        OutM[Ipnt,(4:(3+GeneralSpecs$Nstep))] <- Phi[Ifleet,Iage,Iyear,]
      }
  write(t(OutM),EchoFile,append=T,ncol=3+GeneralSpecs$Nstep)

  
  print("READ IN THE PROJECTION FILE")
  write("READ IN THE PROJECTION FILE\n\n",EchoFile,append=T)

  ReturnObj <- NULL
  ReturnObj$Nproj <- Nproj
  ReturnObj$SelPntFut <- SelPntFut
  ReturnObj$RetPntFut <- RetPntFut
  ReturnObj$LegalFleetPntFut <- LegalFleetPntFut
  ReturnObj$Phi <- Phi

  return(ReturnObj)
  
}

# ===================================================================================

ReadInitialValues <- function(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,
                                 GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
{
 # Main parameters
 # R0, M-bar, M-at-age-offset, WjotyesScaleM, RedsScaleQ, SigmaR
 OK <- 1
 NmainPars = 7+(GeneralSpecs$Nage)+2*GeneralSpecs$Narea;                       
 Index <- MatchTable(ControlFile,Char1="#",Char2="Basic",Char3="parameters");
 MainPars <- rep(0,NmainPars)
 MainBnd <- matrix(0,nrow=NmainPars,ncol=2)
 MainPhase <- rep(NA,NmainPars)
 for (Ipar in 1:NmainPars)
  {
   MainPars[Ipar] <- as.numeric(ControlFile[Index+Ipar,3])
   MainBnd[Ipar,1] <- as.numeric(ControlFile[Index+Ipar,1])
   MainBnd[Ipar,2] <- as.numeric(ControlFile[Index+Ipar,2])
   MainPhase[Ipar] <- as.numeric(ControlFile[Index+Ipar,4])
 }
 write("Starting values for main parameters",EchoFile,append=T)
 write(MainPars,EchoFile,append=T)
 if(is.na(sum(MainPars))) { print('There is a NA in MainPars'); OK <- 0   }  
 
 # Recruitment estimation
 Index <- MatchTable(ControlFile,Char1="#",Char2="Recruitment_deviations"); 
 RecPhase <- as.numeric(ControlFile[Index+3,1])
 write(paste("Phase for recruitment estimates",RecPhase),EchoFile,append=T)
 NrecDev <- ControlSpecs$RecYr2-ControlSpecs$RecYr1+1
 RecDevBnd <- matrix(0,nrow=NrecDev,ncol=2)
 RecDevBnd[,1] <- -15
 RecDevBnd[,2] <-  15
 RecDevPhase <- rep(RecPhase,NrecDev)
 
 Index <- MatchTable(ControlFile,Char1="#",Char2="Spatial_deviations_in_recruitment"); 
 RecSpatPhase <- as.numeric(ControlFile[Index+3,1])
 write(paste("Phase for spatial recruitment estimates",RecSpatPhase),EchoFile,append=T)
 NrecSpatDev <- (ControlSpecs$RecSpatYr2-ControlSpecs$RecSpatYr1+1)*(GeneralSpecs$Narea-1)
 RecSpatDevBnd <- matrix(0,nrow=NrecSpatDev,ncol=2)
 RecSpatDevBnd[,1] <- -15
 RecSpatDevBnd[,2] <-  15
 RecSpatDevPhase <- rep(RecSpatPhase,NrecSpatDev)
 
 Index <- MatchTable(ControlFile,Char1="#",Char2="Prespecify_rec_devs"); 
 PreSpecifyRecDevs <- as.numeric(ControlFile[Index+1,1])
 RecDevPars <- rep(0,NrecDev)
 if (PreSpecifyRecDevs == 1)
  for (Dyr in 1: NrecDev) RecDevPars[Dyr] <-  as.numeric(ControlFile[Index+1+Dyr,1])
 #print(RecDevPars);
 if(is.na(sum(RecDevPars))) { print('There is a NA in RecDevPars'); OK <- 0   }  
 
 Index <- MatchTable(ControlFile,Char1="#",Char2="Prespecify_spatial_rec_devs"); 
 PreSpecifySpatRecDevs <- as.numeric(ControlFile[Index+1,1])
 RecSpatDevPars <- rep(0,NrecSpatDev)
 if (PreSpecifySpatRecDevs == 1){
   for (Dyr in 1: NrecSpatDev) RecSpatDevPars[Dyr] <-  as.numeric(ControlFile[Index+1+Dyr,1])}
 if (NrecSpatDev==0) { NrecSpatDev <- 1; RecSpatDevPars = 0; RecSpatDevBnd <- matrix(c(-15,15),ncol=2,nrow=2); RecSpatDevPhase <- -1; } 
 if(is.na(sum(PreSpecifySpatRecDevs))) { print('There is a NA in PreSpecifySpatRecDevs'); OK <- 0   }  
 
 Index <- MatchTable(SelexFile,Char1="#",Char2="Selectivity",Char3="parameters")+1;
 NlenSel <- max(1,SelexSpecs$NselPars)
 SelPars <- rep(0,NlenSel)
 SelBnd <- matrix(0,NlenSel,ncol=2)
 SelPhase <- rep(NA,NlenSel)
 if (SelexSpecs$NselPars > 0)
  {  
   for (Ipar in 1:SelexSpecs$NselPars)
    {
     SelPars[Ipar] <- as.numeric(SelexFile[Index+Ipar,3])
     SelBnd[Ipar,1] <- as.numeric(SelexFile[Index+Ipar,1])
     SelBnd[Ipar,2] <- as.numeric(SelexFile[Index+Ipar,2])
     SelPhase[Ipar] <- as.numeric(SelexFile[Index+Ipar,4])
    }
  } else { SelPhase[1] <- -100  }
 write("Initial selectivity parameters",EchoFile,append=T)
 if (SelexSpecs$NselPars>0) write(SelPars,EchoFile,append=T,ncol=SelexSpecs$NselPars)
 if(is.na(sum(SelPars))) { print('There is a NA in SelPars'); OK <- 0   }  
 
 Index <- MatchTable(RetenFile,Char1="#",Char2="Retention",Char3="parameters")+1;
 NlenRet <- max(1,RetenSpecs$NretPars)
 RetPars <- rep(0,NlenRet)
 RetBnd <- matrix(0,nrow=RetenSpecs$NretPars,ncol=2)
 RetPhase <- rep(NA,NlenRet)
 if (RetenSpecs$NretPars > 0)
 {  
   for (Ipar in 1:RetenSpecs$NretPars)
    {
     RetPars[Ipar] <- as.numeric(RetenFile[Index+Ipar,3])
     RetBnd[Ipar,1] <- as.numeric(RetenFile[Index+Ipar,1])
     RetBnd[Ipar,2] <- as.numeric(RetenFile[Index+Ipar,2])
     RetPhase[Ipar] <- as.numeric(RetenFile[Index+Ipar,4])
   }
 } else {RetPhase[1] <- -100}
 write("Initial retention parameters",EchoFile,append=T)
 if (RetenSpecs$NretPars>0) write(RetPars,EchoFile,append=T,ncol=RetenSpecs$NretPars)
 if(is.na(sum(RetPars))) { print('There is a NA in RetentionPars'); OK <- 0   }  
 
 Index <- MatchTable(RecruitFile,Char1="#",Char2="Recuitment1",Char3="parameters")+1;
 NRecPars <- ifelse(RecruitSpecs$CalcRecruitFrac==0, RecruitSpecs$NrecruitPars, RecruitSpecs$NrecruitPars+RecruitSpecs$NfixedRecruits*2)
 RecruitPars <- rep(0,NRecPars)
 RecruitBnd <- matrix(0,nrow=NRecPars,ncol=2)
 RecruitPhase <- rep(NA,NRecPars)
 for (Ipar in 1:NRecPars)
  {
   RecruitPars[Ipar] <- as.numeric(RecruitFile[Index+Ipar,3])
   RecruitBnd[Ipar,1] <- as.numeric(RecruitFile[Index+Ipar,1])
   RecruitBnd[Ipar,2] <- as.numeric(RecruitFile[Index+Ipar,2])
   RecruitPhase[Ipar] <- as.numeric(RecruitFile[Index+Ipar,4])
  }
 write("Initial recuitment parameters",EchoFile,append=T)
 write(RecruitPars,EchoFile,append=T,ncol=NRecPars)
 if(is.na(sum(RecruitPars))) { print('There is a NA in RecruitmentPars'); OK <- 0   }  
 
 Index <- MatchTable(RecruitFile,Char1="#",Char2="Puerulus",Char3="Power")+1;
 NPuerPow <- as.numeric(RecruitFile[Index,1])
 PuerPowPars <- rep(0,NPuerPow)
 PuerPowBnd <- matrix(0,nrow=NPuerPow,ncol=2)
 PuerPowPhase <- rep(NA,NPuerPow)
 for (Ipar in 1:NPuerPow)
 {
   PuerPowPars[Ipar] <- as.numeric(RecruitFile[Index+Ipar+1,3])
   PuerPowBnd[Ipar,1] <- as.numeric(RecruitFile[Index+Ipar+1,1])
   PuerPowBnd[Ipar,2] <- as.numeric(RecruitFile[Index+Ipar+1,2])
   PuerPowPhase[Ipar] <- as.numeric(RecruitFile[Index+Ipar+1,4])
 }
 write("Puerulus power parameters",EchoFile,append=T)
 write(PuerPowPars,EchoFile,append=T,ncol=NPuerPow)
 if(is.na(sum(PuerPowPars))) { print('There is a NA in PuerlusPowerPars'); OK <- 0   }  
 
 NlenGr <- max(1,GrowthSpecs$NgrowthPars)
 GrowthPars <- rep(0,NlenGr)
 GrowthBnd <- matrix(0,nrow=NlenGr,ncol=2)
 GrowthPhase <- rep(NA,NlenGr)
 if (GrowthSpecs$NgrowthPars >0)
  {
   for (Ipar in 1:GrowthSpecs$NgrowthPars)
    {
     GrowthPars[Ipar] <- as.numeric(GrowthFile[Index+Ipar,3])
     GrowthBnd[Ipar,1] <- as.numeric(GrowthFile[Index+Ipar,1])
     GrowthBnd[Ipar,2] <- as.numeric(GrowthFile[Index+Ipar,2])
     GrowthPhase[Ipar] <- as.numeric(GrowthFile[Index+Ipar,4])
   }
  }  
 else
  GrowthPhase[1] <- -100
 if(is.na(sum(GrowthPars))) { print('There is a NA in GrowthPars'); OK <- 0   }  
 
 Index <- MatchTable(MoveFile,Char1="#",Char2="Movement",Char3="parameters")+1;
 MovePars <- rep(0,MoveSpecs$NmovePars)
 MoveBnd <- matrix(0,nrow=MoveSpecs$NmovePars,ncol=2)
 MovePhase <- rep(NA,MoveSpecs$NmovePars)
 for (Ipar in 1:MoveSpecs$NmovePars)
  {
   MovePars[Ipar] <- as.numeric(MoveFile[Index+Ipar,3])
   MoveBnd[Ipar,1] <- as.numeric(MoveFile[Index+Ipar,1])
   MoveBnd[Ipar,2] <- as.numeric(MoveFile[Index+Ipar,2])
   MovePhase[Ipar] <- as.numeric(MoveFile[Index+Ipar,4])
  }
 write("Initial movement parameters",EchoFile,append=T)
 write(MovePars,EchoFile,append=T)
 if(is.na(sum(MovePars))) { print('There is a NA in MovePars'); OK <- 0   }  
 
 
 if (ControlSpecs$InitOpt==0||ControlSpecs$InitOpt==1||ControlSpecs$InitOpt==3||ControlSpecs$InitOpt==5) NInitPar <- GeneralSpecs$Narea*GeneralSpecs$Nage*(GeneralSpecs$Nlen[1]+GeneralSpecs$Nlen[2]);
 if (ControlSpecs$InitOpt==2||ControlSpecs$InitOpt==4) NInitPar <- GeneralSpecs$Narea*(GeneralSpecs$Nlen[1]+GeneralSpecs$Nlen[2]);

# Phase for log(initial_N)
 if (ControlSpecs$InitOpt %in% c(0,1,2))   {Ipar = (1+6+(GeneralSpecs$Narea+GeneralSpecs$Nage)):(6+(GeneralSpecs$Narea*2+GeneralSpecs$Nage)); MainPhase[Ipar] <- -1; }
 # Phase for log(initial F)
 if (ControlSpecs$InitOpt %in% c(0,1,2,3,5))   {MainPhase[length(MainPhase)] <- -1; }

 #print(NInitPar)
 Index <- MatchTable(ControlFile,Char1="#",Char2="Initial",Char3="size",Char4="parameters");
 InitPars <- rep(0,NInitPar)
 InitParsBnd <- matrix(0,nrow=NInitPar,ncol=2)
 InitParsPhase <- rep(NA,NInitPar)
 for (Ipar in 1:NInitPar)
   {
    if (ControlSpecs$InitParSpec==0) { Jpar <- Ipar; } else { Jpar <- 1 }
    InitPars[Ipar] <- as.numeric(ControlFile[Index+Jpar,3])
    InitParsBnd[Ipar,1] <- as.numeric(ControlFile[Index+Jpar,1])
    InitParsBnd[Ipar,2] <- as.numeric(ControlFile[Index+Jpar,2])
    InitParsPhase[Ipar] <- as.numeric(ControlFile[Index+Jpar,4])
    if (ControlSpecs$InitOpt==0 || ControlSpecs$InitOpt==4) InitParsPhase[Ipar] <- -1
   } 
 if(is.na(sum(InitPars))) { print('There is a NA in InitialPars'); OK <- 0   }  
 
 Index <- MatchTable(ControlFile,Char1="#",Char2="Q",Char3="parameters");
 NQ <- max(1,GeneralSpecs$NQpars)
 QPars <- rep(0,NQ)
 QBnd <- matrix(0,nrow=NQ,ncol=2)
 QPhase <- rep(NA,NQ)
 if (GeneralSpecs$NQpars >0)
  {
   for (Ipar in 1:GeneralSpecs$NQpars)
    {
     QPars[Ipar] <- as.numeric(ControlFile[Index+Ipar,3])
     QBnd[Ipar,1] <- as.numeric(ControlFile[Index+Ipar,1])
     QBnd[Ipar,2] <- as.numeric(ControlFile[Index+Ipar,2])
     QPhase[Ipar] <- as.numeric(ControlFile[Index+Ipar,4])
    }
  }   else { QPhase[1] <- -100}
 if(is.na(sum(QPars))) { print('There is a NA in QPars'); OK <- 0   }  
 
 Index <- MatchTable(ControlFile,Char1="#",Char2="Efficiency",Char3="parameters") +1;
 Nef <- max(1,Data$NefficPar)
 efPars <- rep(0,Nef)
 efBnd <- matrix(0,nrow=Nef,ncol=2)
 efPhase <- rep(-10,Nef)
 if (Nef>0)  {
   for (Ipar in 1:Nef)    {
     efPars[Ipar] <- as.numeric(ControlFile[Index+Ipar,3])
     efBnd[Ipar,1] <- as.numeric(ControlFile[Index+Ipar,1])
     efBnd[Ipar,2] <- as.numeric(ControlFile[Index+Ipar,2])
     efPhase[Ipar] <- as.numeric(ControlFile[Index+Ipar,4])
   }
 }  
 else { efPhase[1] <- -100 }
 if(is.na(sum(efPars))) { print('There is a NA in Efficiency Pars'); OK <- 0   } 

 if(OK==1) { print("All initial parameters are numerics no NAs were found")}
 ReturnObj <- NULL
 ReturnObj$MainPars$Initial <- MainPars
 ReturnObj$MainPars$Bnd <- MainBnd
 ReturnObj$MainPars$Phase <- MainPhase
 ReturnObj$RecruitPars$Initial <- RecruitPars
 ReturnObj$RecruitPars$Bnd <- RecruitBnd
 ReturnObj$RecruitPars$Phase <- RecruitPhase
 ReturnObj$PuerPowPars$Initial <- PuerPowPars
 ReturnObj$PuerPowPars$Bnd <- PuerPowBnd
 ReturnObj$PuerPowPars$Phase <- PuerPowPhase
 ReturnObj$SelPars$Initial <- SelPars
 ReturnObj$SelPars$Bnd <- SelBnd
 ReturnObj$SelPars$Phase <- SelPhase
 ReturnObj$RetPars$Initial <- RetPars
 ReturnObj$RetPars$Bnd <- RetBnd
 ReturnObj$RetPars$Phase <- RetPhase
 ReturnObj$RecDevs$Initial <- RecDevPars
 ReturnObj$RecDevs$Bnd <- RecDevBnd
 ReturnObj$RecDevs$Phase <- RecDevPhase
 ReturnObj$RecSpatDevs$Initial <- RecSpatDevPars
 ReturnObj$RecSpatDevs$Bnd <- RecSpatDevBnd
 ReturnObj$RecSpatDevs$Phase <- RecSpatDevPhase
 ReturnObj$Qpars$Initial <- QPars
 ReturnObj$Qpars$Bnd <- QBnd
 ReturnObj$Qpars$Phase <- QPhase
 ReturnObj$efpars$Initial <- efPars
 ReturnObj$efpars$Bnd <- efBnd
 ReturnObj$efpars$Phase <- efPhase
 ReturnObj$InitPars$Initial <- InitPars
 ReturnObj$InitPars$Bnd <- InitParsBnd
 ReturnObj$InitPars$Phase <- InitParsPhase
 ReturnObj$MovePars$Initial <- MovePars
 ReturnObj$MovePars$Bnd <- MoveBnd
 ReturnObj$MovePars$Phase <- MovePhase
 ReturnObj$GrowthPars$Initial <- GrowthPars
 ReturnObj$GrowthPars$Bnd <- GrowthBnd
 ReturnObj$GrowthPars$Phase <- GrowthPhase
 return(ReturnObj)
  
}


Parssolved <- function(InitialVars){
  suppressWarnings(rm(parsolve, pos=1 ))
  MaxPhase <<- 1
  for(n in 1:length(unique(names(InitialVars)))){
    par <- unique(names(InitialVars))[n]
    tmp_n <- InitialVars[par][[1]]$Phase
    pos <- which(tmp_n>0)
    if(max(tmp_n)>MaxPhase) MaxPhase <<- max(tmp_n)
    if(length(pos)>0) {
      if(par!="MainPars") pos <-1
      if(!exists('parsolve')) {parsolve <- data.frame(Parameter=par,Position=pos)} else parsolve <- rbind(parsolve,data.frame(Parameter=par,Position=pos))
    }}
  if(!exists('parsolve')) {parsolve <- NA}
  write.table(parsolve,'Output/Parameters_solved.txt',quote = F, sep='\t',row.names = F)
  print(parsolve)
}


toXL <- function(x, rnames=FALSE){
  write.table(x, "clipboard", sep="\t", row.names = rnames)
}

