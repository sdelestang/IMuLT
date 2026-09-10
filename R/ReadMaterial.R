
#' Test Index Value in Data Matrix
#'
#' Internal function to validate index values in pointer matrices during data parsing.
#'
#' @param Data Data frame or matrix to check
#' @param Row Row index to check
#' @param Col Column index to check
#' @param Ind Expected index value (will be compared to Ind-1)
#'
#' @return NULL. Stops with error if validation fails.
#' @keywords internal
testindex <- function(Data,Row,Col,Ind) {
  if(as.numeric(Data[Row,Col])!=(Ind-1)){
    stop(paste("Col",Col,"is wrong in selectivity pointer matrix, line",Row),call. = F) }}

#' Test if Value Can Be Converted to Numeric
#'
#' Internal helper to check if a value is numeric or can be coerced to numeric.
#'
#' @param x Value to test
#' @return Logical vector indicating whether each element can be converted to numeric
#' @keywords internal
asnum <- function(x) suppressWarnings(!is.na(as.numeric(x)))

#' Check Data Object for Missing Values
#'
#' Internal function that runs through data checking for NAs and reports errors.
#' Used during data loading to validate input files.
#'
#' @param x Data object to check (vector, list, data frame)
#' @param i Index in Data list (for error messages)
#' @return NULL. Prints messages if NAs found or data missing.
#' @keywords internal
isnafunc <- function(x,i){
  if(length(x)>0) {
    if(is.na(sum(x))){
      print(paste('There is a NA in', names(Data)[i])) }
  } else {print(paste('There is no data for', names(Data)[i]))   }
}

isnafunc2 <- function(){
  natmp <- NULL; blanktmp <- NULL
  for(i in 1:length(Data)){
    x <- Data[[i]]
  if(length(x)>0) {
    if(is.na(sum(x))){ natmp <- c(natmp, names(Data)[i])   }
  } else {  blanktmp <- c(blanktmp, names(Data)[i])   }
  }
  return(list(natmp=natmp, blanktmp=blanktmp))
  }

#' Parse STARTER.DAT File
#'
#' Internal function to read and parse the STARTER.DAT file which contains
#' file paths for all other model input files.
#'
#' @param StarterFile Data frame from read.table() of STARTER.DAT
#' @return List containing:
#' \itemize{
#'   \item DataFileName - Path to DATA.DAT
#'   \item ControlFileName - Path to CONTROL.DAT
#'   \item SelexFileName - Path to SELEXSPEC.DAT
#'   \item RetainFileName - Path to RETAINSPEC.DAT
#'   \item RecruitFileName - Path to RECRUITSPEC.DAT
#'   \item GrowthFileName - Path to GROWTHSPEC.DAT
#'   \item MoveFileName - Path to MOVESPEC.DAT
#'   \item TagFileName - Path to TAGFILE.TXT
#'   \item PropFFileName - Path to PROPORTIONS.TXT
#'   \item ProjectionsFileName - Path to PROJECTIONS.DAT
#'   \item MaxPhase - Maximum estimation phase
#' }
#' @keywords internal
ReadStarterFile <- function(StarterFile)
{
  ReturnObj <- NULL
  ReturnObj$DataFileName <- StarterFile[1,1]
  ReturnObj$ControlFileName <- StarterFile[2,1]
  ReturnObj$SelexFileName <- StarterFile[3,1]
  ReturnObj$RetainFileName <- StarterFile[4,1]
  ReturnObj$RecruitFileName <- StarterFile[5,1]
  ReturnObj$ReproFileName <- StarterFile[6,1]
  ReturnObj$GrowthFileName <- StarterFile[7,1]
  ReturnObj$MoveFileName <- StarterFile[8,1]
  ReturnObj$TagFileName <- StarterFile[9,1]
  ReturnObj$PropFFileName <- StarterFile[10,1]
  ReturnObj$ProjectionsFileName <- StarterFile[11,1]
  Index <- MatchTable(StarterFile,Char2="#",Char3="Stop",Char4="after")
  ReturnObj$MaxPhase <- as.numeric(StarterFile[Index,1])
  return(ReturnObj)
}


#' Parse General Model Specifications from DATA.DAT
#'
#' Internal function to extract general model structure and dimensions from the
#' DATA.DAT file. Reads year range, spatial structure, temporal structure, and
#' size class definitions.
#'
#' @param DataFile Data frame from read.table() of DATA.DAT
#'
#' @return List containing model structure specifications:
#' \itemize{
#'   \item Year1, Year2 - First and last years of assessment
#'   \item MaxProjYr - Maximum projection years
#'   \item Nyear - Number of years in assessment
#'   \item Nstep - Number of time steps per year
#'   \item Narea - Number of spatial areas
#'   \item Nage - Number of age classes
#'   \item Nsex - Number of sexes modeled
#'   \item Nfleet - Number of fleets
#'   \item MaxLen - Maximum number of size classes across sexes
#'   \item Nlen - Vector of size classes per sex
#'   \item TimeStepLen - Matrix of time step proportions (year × step)
#'   \item BurnIn - Maximum burn-in period across areas > BurnIn - BurnInVec are years with No F but moving to equlibrium population
#'   \item BurnInVec - Vector of burn-in periods by area that experiance F
#'   \item Num_Iteration - Iterations for initial size structure
#'   \item Tune_Years - Years to tune initial conditions over
#'   \item MidLenBin - Matrix of midpoints of size classes (sex × size)
#'   \item LowLenBin - Matrix of lower bounds of size classes (sex × size)
#' }
#'
#' @details
#' The function uses MatchTable() to locate specific keywords in the data file
#' and extract corresponding values. All extracted specifications are written
#' to the Echo.out file for verification. Midpoints of size classes are
#' calculated as averages of adjacent lower bounds.
#'
#' @keywords internal
ReadGeneralFile <- function(DataFile)
 {
  # Read in the data file
  print("READ IN THE GENERAL FILE")
  Index <- MatchTable(DataFile,Char1="#",Char2="First",Char3="year"); Year1 <- as.numeric(DataFile[Index+1,1])
  Index <- MatchTable(DataFile,Char1="#",Char2="Last",Char3="year"); Year2 <- as.numeric(DataFile[Index+1,1])
  Nyear <- Year2-Year1+1
  write(paste("Year1",Year1),EchoFile,append=T)
  write(paste("Year2",Year2),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Maximum",Char3="projection"); MaxProjYr <- as.numeric(DataFile[Index+1,1])
  write(paste("MaxProjYr",MaxProjYr),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Time",Char3="steps"); Nstep <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of time steps",Nstep),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="areas"); Narea <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of areas",Narea),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Burn-in",Char3="whole"); BurnIn <- as.numeric(DataFile[Index+1,1])
  write(paste("Burn-in full model",BurnIn),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Burn-in",Char3="for"); BurnInVec <- as.numeric(DataFile[Index+1,1:Narea]);
  write(paste("Area specific length of burn-in for F to be applied",BurnInVec),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="sexes"); Nsex <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of sexes",Nsex),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="ages"); Nage <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of ages",Nage),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="fleets"); Nfleet <- as.numeric(DataFile[Index+1,1])
  write(paste("Number of fleets",Nfleet),EchoFile,append=T)
  Index <- MatchTable(DataFile,Char1="#",Char2="Number",Char3="of",Char4="size-classes"); Nlen <- as.numeric(DataFile[Index+1,c(1:Nsex)]); MaxLen <- max(Nlen, na.rm=T)
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


#' Parse Tag-Recapture Data from TAGFILE.TXT
#'
#' Internal function to read and parse tag-recapture data including tag releases,
#' recaptures, reporting rates, and tag loss parameters.
#'
#' @param TagFile Data frame from read.table() of TAGFILE.TXT
#' @param PropFFile Data frame from read.table() of PROPORTIONS.TXT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#' @param DataFile Data frame from read.table() of DATA.DAT (currently unused)
#'
#' @return List containing tag-recapture specifications:
#' \itemize{
#'   \item InitialLoss - Initial tag loss rate at tagging
#'   \item TagLossRate - Long-term chronic tag loss rate
#'   \item NrepSplit - Number of reporting rate categories
#'   \item RepRate - Vector of reporting rates by category
#'   \item FitTagSizes - Vector indicating whether to fit size data
#'   \item NtagLag - Number of tag lag periods
#'   \item NtagGroups - Number of tagging groups/cohorts
#'   \item Year1Tag, Year2Tag - First and last year of tags by group
#'   \item TagYr1, TagYr2 - Overall first and last tag years
#'   \item NyearTags - Total number of years with tags
#'   \item TagRel - Array of tag releases (sex × group × area × year × step × size)
#'   \item TagRec - Array of tag recaptures (sex × group × area × reporting × year × step × size)
#'   \item RecapObs - Array of observed recapture proportions (sex × group × area × reporting × year × step)
#'   \item NrelTotal - Matrix of total releases by sex and group
#'   \item NotReportedObs - Matrix of proportion not reported by sex and group
#'   \item PropRepSplit - Array of reporting category proportions (year × step × area × reporting)
#' }
#'
#' @details
#' The function reads multi-dimensional arrays of tag release and recapture data,
#' accounting for spatial structure, temporal dynamics, size structure, and
#' multiple reporting rate categories. Calculates total releases and proportions
#' not reported for use in model fitting.
#'
#' @note This function is currently not used in LoadData() but may be activated
#'   for models incorporating tag-recapture data.
#'
#' @keywords internal
ReadTagFile <- function(TagFile,PropFFile,GeneralSpecs,DataFile)
{
  print("READ IN THE TAGGING FILE")
  # Tag loss rates
  Index <- MatchTable(TagFile,Char1="#",Char2="IsTagData"); IsTagData <- as.numeric(TagFile[Index+1,1]);
  if(IsTagData==0){
    InitialLoss <- 0;TagLossRate <- 0;NrepSplit <- 1;RepRate <- 0;FitTagSizes <- 0;NtagLag <- 0;NtagGroups <- 1; Year1Tag <- GeneralSpecs$Year1
    TagYr1 <- GeneralSpecs$Year1; TagYr2 <- GeneralSpecs$Year1+2; NyearTags <- TagYr2-TagYr1+1;
    TagRel <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NyearTags,GeneralSpecs$Nstep,GeneralSpecs$MaxLen+1))
    TagRec <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NrepSplit,NyearTags,GeneralSpecs$Nstep,GeneralSpecs$MaxLen+1))
    RecapObs <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NrepSplit,NyearTags,GeneralSpecs$Nstep))
    NrelTotal <- matrix(0,nrow=GeneralSpecs$Nsex,ncol=NtagGroups)
    NotReportedObs <- matrix(0,nrow=GeneralSpecs$Nsex,ncol=NtagGroups)
    PropRepSplit<-array(0,dim=c(NyearTags,GeneralSpecs$Nstep,GeneralSpecs$Narea,NrepSplit));  }

  if(IsTagData==1){
    Index <- MatchTable(TagFile,Char1="#",Char2="Initial",Char3="tagloss"); InitialLoss <- as.numeric(TagFile[Index+1,1]);
    Index <- MatchTable(TagFile,Char1="#",Char2="Longterm",Char3="tagloss"); TagLossRate <- as.numeric(TagFile[Index+1,1]);
    Index <- MatchTable(TagFile,Char1="#",Char2="Number",Char4="reporting"); NrepSplit <- as.numeric(TagFile[Index+1,1]);
    Index <- MatchTable(TagFile,Char1="#",Char2="Reporting",Char3="rates"); RepRate <- as.numeric(TagFile[Index+1,1:NrepSplit]);
    Index <- MatchTable(TagFile,Char1="#",Char2="Use",Char3="size"); FitTagSizes <- as.numeric(TagFile[Index+1,1:NrepSplit]);
    Index <- MatchTable(TagFile,Char1="#",Char2="Number",Char4="tsteps"); NtagLag <- as.numeric(TagFile[Index+1,1]);
    Index <- MatchTable(TagFile,Char1="#",Char2="Release",Char3="areas"); NtagGroups <- as.numeric(TagFile[Index+1,1]);
    Index <- MatchTable(TagFile,Char1="#",Char2="First",Char3="release"); Year1Tag <- as.numeric(TagFile[Index+1,1:NtagGroups])

    TagYr1 <- min(Year1Tag); TagYr2 <- GeneralSpecs$Year2; NyearTags <- TagYr2-TagYr1+1;
    TagRel <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NyearTags,GeneralSpecs$Nstep,GeneralSpecs$MaxLen+1))
    Index <- MatchTable(TagFile,Char1="#",Char2="Release",Char4="lbin");
    NRtags <- as.numeric(TagFile[Index+1,1])
    for (i in 1:NRtags) {
      Isex <- as.numeric(TagFile[Index+2+i,1])
      Igrp <- as.numeric(TagFile[Index+2+i,2])
      Iarea <- as.numeric(TagFile[Index+2+i,2])
      Iyear <- as.numeric(TagFile[Index+2+i,3])-TagYr1+1
      Istep <- as.numeric(TagFile[Index+2+i,4])
      TagRel[Isex,Igrp,Iarea,Iyear,Istep,] <- as.numeric(TagFile[Index+2+i,5:(5+GeneralSpecs$MaxLen)]);
  }

    Index <- MatchTable(TagFile,Char1="#",Char2="Recaptures",Char4="lbin");
    TagRec <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NrepSplit,NyearTags,GeneralSpecs$Nstep,GeneralSpecs$MaxLen+1))
    NRtags <- as.numeric(TagFile[Index+1,1])
    for (i in 1:NRtags) {
      Isex <- as.numeric(TagFile[Index+2+i,1])
      Igrp <- as.numeric(TagFile[Index+2+i,2])
      Iarea <- as.numeric(TagFile[Index+2+i,3])
      Itype <- as.numeric(TagFile[Index+2+i,4])
      Iyear <- as.numeric(TagFile[Index+2+i,5])-TagYr1+1
      Istep <- as.numeric(TagFile[Index+2+i,6])
      TagRec[Isex,Igrp,Iarea,Itype,Iyear,Istep,] <- as.numeric(TagFile[Index+2+i,7:(7+GeneralSpecs$MaxLen)]);
    }

    Index <- MatchTable(TagFile,Char1="#",Char2="Recaptures",Char4="timestep");
    RecapObs <- array(0,dim=c(GeneralSpecs$Nsex,NtagGroups,GeneralSpecs$Narea,NrepSplit,NyearTags,GeneralSpecs$Nstep))
    NRtags <- as.numeric(TagFile[Index+1,1])
    for (i in 1:NRtags) {
      Isex <- as.numeric(TagFile[Index+2+i,1])
      Igrp <- as.numeric(TagFile[Index+2+i,2])
      Iarea <- as.numeric(TagFile[Index+2+i,3])
      Itype <- as.numeric(TagFile[Index+2+i,4])
      Iyear <- as.numeric(TagFile[Index+2+i,5])-TagYr1+1
      RecapObs[Isex,Igrp,Iarea,Itype,Iyear,] <- as.numeric(TagFile[Index+2+i,6:(5+GeneralSpecs$Nstep)]);
    }

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

  for (Isex in 1:GeneralSpecs$Nsex){
    for (Igrp in 1:NtagGroups){
      NotReportedObs[Isex,Igrp] <- (NrelTotal[Isex,Igrp]-NotReportedObs[Isex,Igrp])/NrelTotal[Isex,Igrp]}}

  for (Isex in 1:GeneralSpecs$Nsex){
    for (Igrp in 1:NtagGroups){
      for (Iarea in 1:GeneralSpecs$Narea){
        for (Iyear in 1:NyearTags){
          for (Istep in 1:GeneralSpecs$Nstep){
            for (Irep in 1:NrepSplit){
              RecapObs[Isex,Igrp,Iarea,Irep,Iyear,Istep] = TagRec[Isex,Igrp,Iarea,Irep,Iyear,Istep,1]/NrelTotal[Isex,Igrp]
              }}}}}}

  PropRepSplit<-array(0,dim=c(NyearTags,GeneralSpecs$Nstep,GeneralSpecs$Narea,NrepSplit));
  Nobs <-  as.numeric(PropnFile[1,1])
  for (i in 1:Nobs) {
    Iyear <- as.numeric(PropnFile[i+2,1])-TagYr1+1
    Istep <- as.numeric(PropnFile[i+2,2])
    Iarea <- as.numeric(PropnFile[i+2,3])
    PropRepSplit[Iyear,Istep,Iarea,]  <- as.numeric(PropnFile[i+2,4:(3+NrepSplit)]);
      }

  write("READ IN THE TAG FILE\n\n",EchoFile,append=T)
  }
ReturnObj <- NULL
ReturnObj$IsTagData <- IsTagData
ReturnObj$InitialLoss <- InitialLoss
ReturnObj$TagLossRate <- TagLossRate
ReturnObj$NrepSplit <- NrepSplit
ReturnObj$RepRate <- RepRate
ReturnObj$FitTagSizes <- FitTagSizes
ReturnObj$NtagLag <- NtagLag
ReturnObj$NtagGroups <- NtagGroups
ReturnObj$Year1Tag <- Year1Tag
ReturnObj$Year2Tag <- GeneralSpecs$Year2
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

#' Parse Fishery Data from DATA.DAT
#'
#' Internal function to read and parse all fishery-dependent and fishery-independent
#' data from DATA.DAT including catch, CPUE indices, length compositions, larval
#' settlement data, and environmental covariates.
#'
#' @param DataFile Data frame from read.table() of DATA.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#'
#' @return List containing all data specifications:
#' \itemize{
#'   \item Catch - Array of catches (year × step × fleet)
#'   \item Ncpue - Number of CPUE observations
#'   \item NcpueDataSeries - Number of CPUE data series
#'   \item IndexType - Vector indicating index type (1=weight, 2=numbers)
#'   \item FixSigmaCpue - Vector for sigma treatment by series
#'   \item SigmaCpueOffset - Minimum sigma for CPUE
#'   \item TreatQcpue - Vector for catchability treatment
#'   \item EnvIndCpue - Vector linking series to environmental indices
#'   \item EffCrIndCpue - Vector for efficiency creep indices
#'   \item EffCrLag - Temporal lag structure for efficiency creep
#'   \item IndexI - Matrix of CPUE integer specifications (series, fleet, sex, year, step)
#'   \item IndexR - Matrix of CPUE real values (index, CV)
#'   \item LarvalLikeOpt - Likelihood for larval data (0=lognormal, else normal)
#'   \item Larval_Offset - Time delay from settlement to recruitment (years)
#'   \item NLarvalData - Number of larval/puerulus observations
#'   \item Lar_dataI - Matrix of larval integer specs (area, year)
#'   \item Lar_dataR - Matrix of larval real values (index, SD)
#'   \item NenvSeries - Number of environmental data series
#'   \item EnvData - Array of environmental covariates (year × step × series)
#'   \item NQparPass - Number of Q-related parameters
#'   \item Nnumbers - Number of catch-in-numbers observations
#'   \item NcatchDataSeries - Number of catch-in-numbers series
#'   \item FixSigmaCatchN - Sigma treatment for numbers data
#'   \item SigmaCatchNOffset - Minimum sigma for numbers
#'   \item NumbersI - Matrix of numbers integer specs (series, fleet, year, step)
#'   \item NumbersR - Matrix of numbers real values (catch, CV)
#'   \item NlenComp - Number of length composition samples
#'   \item LenCompI - Matrix of length comp integer specs (fleet, sex, year, step)
#'   \item Stage1W - Vector of stage 1 weights for length comps
#'   \item LenCompR - Matrix of length composition proportions (sample × size)
#' }
#'
#' @details
#' Reads all observational data used for model fitting. Data are parsed from
#' specific sections of DATA.DAT identified by keyword matching. All data are
#' written to Echo.out for verification. Length compositions are normalized to
#' sum to 1.0 within each sample.
#'
#' @keywords internal
ReadDataFile <- function(DataFile,GeneralSpecs)
{
  print("READ IN THE DATA FILE")
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
  write(Catch,EchoFile,append=T,ncolumns=(GeneralSpecs$Nyear+GeneralSpecs$MaxProjYr))

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
  SigmaCpueCeiling <- as.numeric(DataFile[Index+18,1]);
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
  write("Size-composition data",EchoFile,append=T)
  write(t(cbind(LenCompI,LenCompR)),EchoFile,append=T,ncol=GeneralSpecs$MaxLen+4)

  # Read in the larval data
  Index <- MatchTable(DataFile,Char1="#",Char2="Larval",Char3="index")
  LarvalLikeOpt <- as.numeric(DataFile[Index+2,1]);                        # 0 for log-normal; otherwise normal
  Larval_Offset <- as.numeric(DataFile[Index+4,1]);                        # Time-delay between settlement and recruitment
  NLarvalData <- as.numeric(DataFile[Index+6,1]);
  Lar_dataI <- matrix(0,nrow=NLarvalData,ncol=2);
  Lar_dataR <- matrix(0,nrow=NLarvalData,ncol=2);
  Index <- Index + 7
  write("Puerulus data",EchoFile,append=T)
  if (NLarvalData > 0){
    for (Idata in 1:NLarvalData)
    {
      for (II in 1:2) Lar_dataI[Idata,II] <- as.numeric(DataFile[Index+Idata,II]) - c(1,0)[II]
      for (II in 1:2) Lar_dataR[Idata,II] <- as.numeric(DataFile[Index+Idata,II+2])
      Lar_dataI[Idata,2] <- Lar_dataI[Idata,2]  - GeneralSpecs$Year1 + max(GeneralSpecs$BurnIn)
    }
    if (asnum(DataFile[Index+Idata+1,II+2])) { print("Error reading Larval data; too many inputs: Stopping"); AAA }

    write(t(cbind(Lar_dataI,Lar_dataR)),EchoFile,append=T,ncol=4)

  }

  # Read in the environmental data
  Index <- MatchTable(DataFile,Char1="#",Char2="Environmental",Char3="Data")
  NenvSeries <- as.numeric(DataFile[Index+2,1]); if(NenvSeries==0) NenvSeries <- 1
  write("\nEnvironmental data",EchoFile,append=T)
  write(paste("Number of environmental series",NenvSeries),EchoFile,append=T)
  YearsPerSeries <- as.numeric(DataFile[Index+4,1:NenvSeries])
  EnvData <- array(0,dim=c(GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1,GeneralSpecs$Nstep,NenvSeries))
  Ipnt <- Index+5
  if(as.numeric(DataFile[Index+2,1])>0){
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
    } }
  if (asnum(DataFile[Ipnt+1,3])) { print("Error reading Enviromental data; too many inputs: Stopping"); AAA }

  write("READ IN THE DATA FILE\n\n",EchoFile,append=T)

  ReturnObj <- NULL
  ReturnObj$Catch <- Catch
  ReturnObj$Ncpue <- Ncpue
  ReturnObj$NcpueDataSeries <- NcpueDataSeries
  ReturnObj$IndexType <- IndexType
  ReturnObj$FixSigmaCpue <- FixSigmaCpue
  ReturnObj$SigmaCpueOffset <- SigmaCpueOffset
  ReturnObj$SigmaCpueCeiling <- SigmaCpueCeiling
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


#' Parse Control File Specifications from CONTROL.DAT
#'
#' Internal function to read and parse model control parameters, biological
#' specifications, fleet-area linkages, data weights, and parameter linking
#' options from CONTROL.DAT.
#'
#' @param ControlFile Data frame from read.table() of CONTROL.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#' @param DataSpecs List from ReadDataFile() containing data specifications
#'
#' @return List containing control specifications:
#' \itemize{
#'   \item Fleet_area - Vector linking fleets to areas
#'   \item Narea_fleet - Number of fleets per area
#'   \item Area_fleet - Matrix of fleet-area associations (area × fleet)
#'   \item Nzone - Number of management zones
#'   \item NareasPerZone - Number of areas per zone
#'   \item AreasPerZone - Matrix of area-zone associations (zone × area)
#'   \item NefficPar - Number of efficiency creep parameters
#'   \item Phi1 - Array of discard mortality rates (fleet × age × year × step)
#'   \item WeightLen - Matrix of weight-at-length (sex × size)
#'   \item MatFem - Matrix of egg production at length (area × size)
#'   \item MatAge - Vector of age at maturity by area
#'   \item MatTimeStep - Time step for egg production calculation
#'   \item BioTimeStep - Time step for biomass calculation
#'   \item RecYr1, RecYr2 - First and last years for recruitment deviations
#'   \item RecSpatYr1, RecSpatYr2 - First and last years for spatial recruitment deviations
#'   \item MparsLink - Vector of parameter linking specifications
#'   \item MparsPrior - Matrix of prior specifications (mean, SD, type)
#'   \item InitOpt - Option for initial conditions calculation
#'   \item InitParSpec - Option for initial parameter values
#'   \item LambdaCpue - Global weight on CPUE data
#'   \item LambdaNumbers - Global weight on catch numbers data
#'   \item LambdaLength - Global weight on length composition data
#'   \item LambdaLarval - Global weight on larval data
#'   \item LambdaTag1, LambdaTag2 - Global weights on tag data
#'   \item LambdaCpue2 - Vector of CPUE weights by series
#'   \item LambdaNumbers2 - Vector of numbers weights by series
#'   \item LambdaLength2 - Array of length comp weights (fleet × step × sex)
#'   \item WeightInitialN - Weight on initial numbers penalty
#'   \item WeightInit3 - Additional initial conditions weight
#'   \item NvarTypes - Number of variance parameter types
#'   \item VarTypes - Vector indicating which variance parameters to estimate
#' }
#'
#' @details
#' Reads all model control parameters including biological parameters (growth,
#' maturity, fecundity), fleet specifications, management zones, discard mortality
#' rates, data weighting schemes, and parameter estimation controls. Supports
#' both global and detailed (fleet/sex/step specific) data weighting. All
#' specifications are written to Echo.out for verification.
#'
#' @keywords internal
ReadControlFile <- function(ControlFile,GeneralSpecs,DataSpecs)
{
  print("READ IN THE CONTROL FILE")
  write("READING IN THE CONTROL FILE",EchoFile,append=T)

  # Weight-length regression
  Index <- MatchTable(ControlFile,Char1="#",Char2="weight-at-length",Char3=NULL)+1;
  WeightLen <- matrix(0,GeneralSpecs$Nsex,GeneralSpecs$MaxLen)
  for (Isex in 1:GeneralSpecs$Nsex)
  {
    for (Jlen in 1:GeneralSpecs$Nlen[Isex]) WeightLen[Isex,Jlen] <- as.numeric(ControlFile[Index,Jlen])
    if (!is.na(ControlFile[Index,GeneralSpecs$Nlen[Isex]+1]) & ControlFile[Index,GeneralSpecs$Nlen[Isex]+1]!="") { print("Error reading weight-at-length; too many inputs: Stopping"); AAA }
    Index <- Index + 1
  }
  write("Weight-length regressions",EchoFile,append=T)
  write(t(WeightLen),EchoFile,append=T,ncol=GeneralSpecs$MaxLen)

  Index <- MatchTable(ControlFile,Char1="#",Char2="Egg",Char3="time")+1;
  MatTimeStep <- as.numeric(ControlFile[Index,1])

  Index <- MatchTable(ControlFile,Char1="#",Char2="Biomass",Char3="time")+1;
  BioTimeStep <- as.numeric(ControlFile[Index,1])

  # Link between fleets and areas
  Index <- MatchTable(ControlFile,Char1="#",Char2="Fleet",Char3="Area");
  Fleet_area <- rep(0,GeneralSpecs$Nfleet)
  Fleet_name <- rep(0,GeneralSpecs$Nfleet)
  for (Ifleet in 1:GeneralSpecs$Nfleet) {
    Fleet_area[Ifleet] <- as.numeric(ControlFile[Index+Ifleet,2])
    Fleet_name[Ifleet] <- ControlFile[Index+Ifleet,3]}
  Fleet_name[is.na(Fleet_name)] <- 'Not_provided'
  if (asnum(ControlFile[Index+Ifleet+1,2])) { print("Error reading Fleet data; too many inputs: Stopping"); AAA }
  write("Area for each fleet",EchoFile,append=T)
  write(Fleet_area,EchoFile,append=T)
  write("fleet name",EchoFile,append=T)
  write(Fleet_name,EchoFile,append=T)

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
        testindex(ControlFile,pos,1,Iage)
        testindex(ControlFile,pos,2,Ifleet)
        testindex(ControlFile,pos,3,Istep)
        ControlFile[Index+Ipnt,ControlFile[Index+Ipnt,]==""|is.na(ControlFile[Index+Ipnt,])] <- NA
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

  # Main parameters linking conditions
  Index <- MatchTable(ControlFile,Char1="#",Char2="Basic");
  npars <- 1*GeneralSpecs$Narea+GeneralSpecs$Nage+4
  MparsLink <- as.numeric(ControlFile[(Index+1):(Index+npars),5])
  MparsPrior <- apply(as.matrix(ControlFile[(Index+1):(Index+npars),6:8]),2,as.numeric)

  # Efficiency parameters linking conditions
  Index <- MatchTable(ControlFile,Char1="#",Char2="Efficiency",Char3="parameters");
  npars <- as.numeric(ControlFile[(Index+1),1])
  EffparsLink <- as.numeric(ControlFile[(Index+2):(Index+npars+1),5])
  EffparsPrior <- apply(as.matrix(ControlFile[(Index+2):(Index+npars+1),6:8]),2,as.numeric)


  # read the data weights
  Index <- MatchTable(ControlFile,Char1="#",Char2="Weights",Char3="on");
  LambdaCpue <- as.numeric(ControlFile[Index+1,1])
  LambdaNumbers <- as.numeric(ControlFile[Index+2,1])
  LambdaLength <- as.numeric(ControlFile[Index+3,1])
  LambdaLarval <- as.numeric(ControlFile[Index+4,1])
  LambdaTag1 <- as.numeric(ControlFile[Index+5,1])
  LambdaTag2 <- as.numeric(ControlFile[Index+6,1])
  # WeightInitialN <- as.numeric(ControlFile[Index+7,1])
  # WeightInit3 <- as.numeric(ControlFile[Index+8,1])

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
        Fleet <- TheSpec[2]+1; IdataSet <- unique(Data$IndexI[Data$IndexI[,2]==Fleet,1]); Wght <- TheSpec[5]
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
  #ReturnObj$MatFem <- MatFem
  #ReturnObj$MatAge <- MatAge
  ReturnObj$MatTimeStep <- MatTimeStep
  ReturnObj$BioTimeStep <- BioTimeStep
  ReturnObj$RecYr1 <- RecYr1
  ReturnObj$RecYr2 <- RecYr2
  ReturnObj$RecSpatYr1 <- RecSpatYr1
  ReturnObj$RecSpatYr2 <- RecSpatYr2
  ReturnObj$MparsLink <- MparsLink
  ReturnObj$MparsPrior <- MparsPrior
  ReturnObj$EffparsLink <- EffparsLink
  ReturnObj$EffparsPrior <- EffparsPrior
  #ReturnObj$InitOpt <- InitOpt
  #ReturnObj$InitParSpec <- InitParSpec
  ReturnObj$LambdaCpue <- LambdaCpue
  ReturnObj$LambdaNumbers <- LambdaNumbers
  ReturnObj$LambdaLength <- LambdaLength
  ReturnObj$LambdaLarval <- LambdaLarval
  ReturnObj$LambdaTag1 <- LambdaTag1
  ReturnObj$LambdaTag2 <- LambdaTag2
  ReturnObj$LambdaCpue2 <- LambdaCpue2
  ReturnObj$LambdaNumbers2 <- LambdaNumbers2
  ReturnObj$LambdaLength2 <- LambdaLength2
  # ReturnObj$WeightInitialN <- WeightInitialN
  # ReturnObj$WeightInit3 <- WeightInit3
  ReturnObj$NvarTypes <- NvarTypes
  ReturnObj$VarTypes <- VarTypes
  return(ReturnObj)

}

#' Parse Reproduction Specifications from REPROD.DAT
#'
#' Internal function to read and parse reproductive biology parameters from
#' REPROD.DAT, including maturity, multiple spawning, and fecundity specifications.
#' Combines these components to construct an egg production array across age,
#' area, year, and length.
#'
#' @param ReprodFile Data frame from read.table() of REPROD.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#' @param DataSpecs List from ReadDataFile() containing data specifications
#'
#' @return List containing reproduction specifications:
#' \itemize{
#'   \item MatFem - Array of egg production (age × area × year × length bin)
#'     combining maturity, multiple spawning, and fecundity
#'   \item MatAge - Numeric vector of maturity-at-age values by area
#' }
#'
#' @details
#' Three reproductive components are read and combined into the egg production array:
#'
#' \strong{Maturity}: Modelled as a logistic function of length with parameters
#' \code{Par_a} (L50) and \code{Par_b} (slope). A pattern index of -1 sets
#' maturity to 1 for all lengths.
#'
#' \strong{Multiple spawning}: Modelled as a scaled logistic function with parameters
#' \code{Par_a} (inflection point), \code{Par_b} (slope), and \code{Par_c} (scalar).
#' Captures variation in spawning frequency across lengths.
#'
#' \strong{Fecundity}: Modelled as a power function of length (\code{Par_a * L ^ Par_b}).
#' A pattern index of -1 sets fecundity to 1 for all lengths.
#'
#' For each component, the input file defines parameter sets and a specification
#' matrix that maps each age × area × year combination to a parameter set index
#' (0-indexed). The final \code{MatFem} array is the product of all three components
#' evaluated at mid-length bin values, and is written to Echo.out.
#'
#' @keywords internal
ReadReprodFile <- function(ReprodFile,GeneralSpecs,DataSpecs)
{
  print("READ IN THE REPRODUCTION FILE")
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Age",Char4="maturity")+1;
  MatAge <- as.integer(ReprodFile[Index,1:GeneralSpecs$Narea])
  write("Maturity at age",EchoFile,append=T)
  write(MatAge,EchoFile,append=T,ncol=GeneralSpecs$Narea)

  ## Get Maturity
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Number",Char4="maturity")+1;
  Nrows <- as.numeric(ReprodFile[Index,1])
  Names <- ReprodFile[(Index+1),2:4]
  Matpars <- (ReprodFile[(Index+2):(Index+1+Nrows),1:3])
  Matpars[, 1:3] <- lapply(Matpars[, 1:3], as.numeric)
  colnames(Matpars) <- Names
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Specifications",Char4="maturity")+1;
  Names <- ReprodFile[Index,2:(GeneralSpecs$Nyear+3)];
  NageArea <- GeneralSpecs$Nage*GeneralSpecs$Narea
  MatSpec <-  ReprodFile[(Index+1):(Index+NageArea),1:(GeneralSpecs$Nyear+2)];

  ## Get MultipleSpawn
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Number",Char4="multiple")+1;
  Nrows <- as.numeric(ReprodFile[Index,1])
  Names <- ReprodFile[(Index+1),2:5]
  Mulpars <- (ReprodFile[(Index+2):(Index+1+Nrows),1:4])
  Mulpars[, 1:4] <- lapply(Mulpars[, 1:4], as.numeric)
  colnames(Mulpars) <- Names
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Specifications",Char4="multiple")+1;
  MulSpec <-  ReprodFile[(Index+1):(Index+NageArea),1:(GeneralSpecs$Nyear+2)];

  ## Get Fecundity
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Number",Char4="fecundity")+1;
  Nrows <- as.numeric(ReprodFile[Index,1])
  Names <- ReprodFile[(Index+1),2:4]
  Fecpars <- (ReprodFile[(Index+2):(Index+1+Nrows),1:3])
  Fecpars[, 1:3] <- lapply(Fecpars[, 1:3], as.numeric)
  colnames(Fecpars) <- Names
  Index <- MatchTable(ReprodFile,Char1="#",Char2="Specifications",Char4="fecundity")+1;
  FecSpec <-  ReprodFile[(Index+1):(Index+NageArea),1:(GeneralSpecs$Nyear+2)];

  nyears <- GeneralSpecs$BurnIn + GeneralSpecs$Nyear + GeneralSpecs$MaxProjYr + 1
  MatFem <- array(1, dim = c(GeneralSpecs$Nage, GeneralSpecs$Narea, nyears, GeneralSpecs$MaxLen))
  Mlbin <- GeneralSpecs$MidLenBin[1,][1:GeneralSpecs$MaxLen]

  for (Iage in 0:(GeneralSpecs$Nage - 1)) {
    for (Iarea in 0:(GeneralSpecs$Narea - 1)) {
      for (IyearLong in 1:nyears) {

        if (IyearLong <= GeneralSpecs$BurnIn) {
          Iyear <- 1
        } else if (IyearLong > GeneralSpecs$BurnIn + GeneralSpecs$Nyear) {
          Iyear <- GeneralSpecs$Nyear
        } else {
          Iyear <- IyearLong - GeneralSpecs$BurnIn
        }

        MatPoint <- as.numeric(MatSpec[MatSpec[,1] == Iage & MatSpec[,2] == Iarea, Iyear + 2])
        Maturity <- rep(1, GeneralSpecs$Nlen[1])
        if (MatPoint >= 0) Maturity <- 1 / (1 + exp((Mlbin - Matpars$Par_a[MatPoint + 1]) / Matpars$Par_b[MatPoint + 1]))

        MulPoint <- as.numeric(MulSpec[MulSpec[,1] == Iage & MulSpec[,2] == Iarea, Iyear + 2])
        Multiple <- rep(1, GeneralSpecs$Nlen[1])
        if (MulPoint >= 0) Multiple <- Mulpars$Par_c[MulPoint + 1] / (1 + exp((Mlbin - Mulpars$Par_a[MulPoint + 1]) / Mulpars$Par_b[MulPoint + 1]))

        FecPoint <- as.numeric(FecSpec[FecSpec[,1] == Iage & FecSpec[,2] == Iarea, Iyear + 2])
        Fecundity <- rep(1, GeneralSpecs$Nlen[1])
        if (FecPoint >= 0) Fecundity <- Fecpars$Par_a[FecPoint + 1] * Mlbin ^ Fecpars$Par_b[FecPoint + 1]

        MatFem[Iage + 1, Iarea + 1, IyearLong, ] <- Maturity * Multiple * Fecundity
      }
    }
  }

  write("Egg Production\nAge Area Year Lbins\n",EchoFile,append=T)
  for (Iage in 1:(GeneralSpecs$Nage)){
    for (Iarea in 1:(GeneralSpecs$Narea)){
      MatFemEcho <- data.frame(age=rep(Iage,nyears),area=rep(Iarea,nyears),year=(GeneralSpecs$Year1-GeneralSpecs$BurnIn):(GeneralSpecs$Year2+GeneralSpecs$MaxProjYr+1))
      write(t(cbind(MatFemEcho,MatFem[Iage,Iarea,,])),EchoFile,append=T,ncol=GeneralSpecs$MaxLen+3)
    }}

  write("READ IN THE REPRODUCTION FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$MatFem <- MatFem
  ReturnObj$MatAge <- MatAge
  return(ReturnObj)
}


#' Parse Movement Specifications from MOVESPEC.DAT
#'
#' Internal function to read and parse spatial movement/migration patterns and
#' parameters from MOVESPEC.DAT. Defines how individuals move between areas as
#' a function of age, year, and season.
#'
#' @param MoveFile Data frame from read.table() of MOVESPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#'
#' @return List containing movement specifications:
#' \itemize{
#'   \item NmovePatterns - Number of movement patterns defined
#'   \item MoveSpec - Matrix of movement pattern specifications (pattern × 4)
#'     \itemize{
#'       \item Column 1: Pattern ID
#'       \item Column 2: Type (0=none, 1=constant, 2=knife-edged)
#'       \item Column 3: Destination area
#'       \item Column 4: Extra specification
#'     }
#'   \item MovePnt - Array of movement pattern pointers (area × age × year × step)
#'     indicating which movement pattern applies for each combination
#'   \item NmovePars - Total number of movement parameters to estimate
#' }
#'
#' @details
#' Movement patterns define spatial redistribution of population between areas.
#' Each pattern can be:
#' \itemize{
#'   \item Type 0: No movement
#'   \item Type 1: Constant movement (1 parameter: proportion moving)
#'   \item Type 2: Knife-edged movement (2 parameters: size threshold and proportion)
#' }
#'
#' Movement pointers link each age/area/year/step combination to a specific
#' movement pattern. The number of parameters is calculated based on the types
#' of movement patterns specified. All specifications are written to Echo.out.
#'
#' @keywords internal
ReadMoveFile <- function(MoveFile,GeneralSpecs)
{
  print("READ IN THE MOVEMENT FILE")
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

  # Movement parameters linking conditions
  Index <- MatchTable(MoveFile,Char1="#",Char2="Movement",Char3="parameters")+2;
  if(NmovePars>0){
    MoveparsLink <- as.numeric(MoveFile[(Index):(Index+NmovePars-1),5])
    MoveparsPrior <- matrix(apply(as.matrix(MoveFile[(Index):(Index+NmovePars-1),6:8]),2,as.numeric), ncol=3)
  } else { MoveparsLink <- 0; MoveparsPrior <- matrix(c(0,0,0),ncol=3) }

  write("READ IN THE MOVEMENT FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NmovePatterns <- NmovePatterns
  ReturnObj$MoveSpec <- MoveSpec
  ReturnObj$MovePnt <- MovePnt
  ReturnObj$NmovePars <- NmovePars
  ReturnObj$MoveparsLink <- MoveparsLink
  ReturnObj$MoveparsPrior <- MoveparsPrior
  return(ReturnObj)
}

#' Parse Selectivity and Legal Size Specifications from SELEXSPEC.DAT
#'
#' Internal function to read and parse selectivity patterns, legal size regulations,
#' and retention specifications from SELEXSPEC.DAT. Defines how fishery selectivity
#' and legal size limits vary by fleet, area, age, sex, year, and season.
#'
#' @param SelexFile Data frame from read.table() of SELEXSPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#'
#' @return List containing selectivity and legal size specifications:
#' \itemize{
#'   \item NselPatterns - Number of selectivity patterns defined
#'   \item SelSpec - Matrix of selectivity pattern specifications (pattern × 5)
#'     \itemize{
#'       \item Column 1: Pattern ID
#'       \item Column 2: Type (2=logistic, 3=double-normal, 4=general)
#'       \item Column 3: First parameter specification
#'       \item Column 4: Number of parameters
#'       \item Column 5: Additional specification
#'     }
#'   \item SelPnt - Array of selectivity pattern pointers (sex × age × fleet × year × step)
#'   \item NselPars - Total number of selectivity parameters to estimate
#'   \item SelparsLink - Vector linking selectivity parameters (for shared parameters)
#'   \item SelexFI - Matrix of fixed selectivity-at-length patterns (pattern × length)
#'   \item NfixedSelex - Number of fixed selectivity patterns
#'   \item NlegalPatterns - Number of legal size patterns defined
#'   \item LegalSpec - Matrix of legal pattern specifications (pattern × 4)
#'   \item LegalPnt - Array of legal pattern pointers by area (sex × age × area × year × step)
#'   \item LegalFleetPnt - Array of legal pattern pointers by fleet (sex × age × fleet × year × step)
#'   \item LegalFI - Matrix of fixed legal-at-length patterns (pattern × length)
#'   \item LegalRef - Vector of reference legal selectivity pattern (length vector)
#'   \item NlegalSelex - Number of fixed legal patterns
#'   \item IsRed - Array indicating setose/red lobster retention (sex × age × area × step)
#' }
#'
#' @details
#' This function handles two main types of specifications:
#'
#' **Selectivity patterns** define the probability of capture given encounter with
#' fishing gear. Patterns can be:
#' \itemize{
#'   \item Type 2: Logistic selectivity (2 parameters: size-at-50%, slope)
#'   \item Type 3: Double-normal selectivity (multiple parameters)
#'   \item Type 4: General parametric forms
#' }
#'
#' **Legal size patterns** define retention probability for legal-sized animals,
#' accounting for minimum size limits and other regulations. These can vary by
#' both fleet (enforcement differences) and area (regional regulations).
#'
#' Fixed patterns provide length-specific selectivity/legal values that are not
#' estimated but held constant. The IsRed specifications handle special retention
#' rules for setose (red) lobsters which may be protected in certain areas/seasons.
#'
#' Pointers link each combination of biological/spatial/temporal factors to specific
#' selectivity or legal patterns, allowing flexible representation of how fishing
#' selectivity changes across the population structure. Parameter linking allows
#' multiple patterns to share the same underlying parameter values.
#'
#' All specifications are written to Echo.out for model verification.
#'
#' @keywords internal
ReadSelexFile <- function(SelexFile,GeneralSpecs)
{
  print("READ IN THE SELEX FILE")
  Index <- MatchTable(SelexFile,Char1="#",Char2="Number",Char3="Selex")+1;
  NselPatterns <- as.numeric(SelexFile[Index,1])
  write(paste("Number of selectivity patterns",NselPatterns),EchoFile,append=T)

  Index <- MatchTable(SelexFile,Char1="#",Char2="Pattern",Char3="Type")[1];
  SelSpec <- matrix(0,NselPatterns,5)
  for (Isel in 1:NselPatterns)
    for (Icol in 1:5) SelSpec[Isel,Icol] <- as.numeric(SelexFile[Index+Isel,Icol])
  write("Specifications for selectivity",EchoFile,append=T)
  write(t(SelSpec),EchoFile,append=T,ncol=5)

  NselPars <- sum(SelSpec[,4])
  # for (Isel in 1:NselPatterns)
  # {
  #   if (SelSpec[Isel,2]%in%c(2,4)) NselPars <- NselPars + SelSpec[Isel,4]
  #   if (SelSpec[Isel,2]== 3) NselPars <- NselPars + 2
  # }
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

  # Selectivity parameters linking conditions
  Index <- MatchTable(SelexFile,Char1="#",Char2="Selectivity",Char3="Parameters")+2;
  SelparsLink <- as.numeric(SelexFile[(Index):(Index+NselPars-1),5])
  SelparsPrior <- apply(as.matrix(SelexFile[(Index):(Index+NselPars-1),6:8]),2,as.numeric)

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
  LegalRef <- matrix(0,nrow=GeneralSpecs$Nsex, ncol=GeneralSpecs$MaxLen)
  for (Jsex in 1:GeneralSpecs$Nsex) {
    for (Jlen in 1:GeneralSpecs$Nlen[1]) {
      LegalRef[Jsex,Jlen] <- as.numeric(SelexFile[Index+Jsex,Jlen])}}
  write("Reference legal pattern",EchoFile,append=T)
  write(LegalRef,EchoFile,append=T,ncol=GeneralSpecs$MaxLen)

  Index <- MatchTable(SelexFile,Char1="#",Char2="IsMorph",Char3="specifications",Char4="-")+2;
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
  write("Specifications for IsMorph0",EchoFile,append=T)
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

  write("READ IN THE SELEX FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NselPatterns <- NselPatterns
  ReturnObj$NlegalPatterns <- NlegalPatterns
  ReturnObj$SelSpec <- SelSpec
  ReturnObj$SelPnt <- SelPnt
  ReturnObj$SelparsLink <- SelparsLink
  ReturnObj$SelparsPrior <- SelparsPrior
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

#' Parse Retention Specifications from RETENSPEC.DAT
#'
#' Internal function to read and parse retention/discard patterns from RETENSPEC.DAT.
#' Defines the probability that captured animals are retained (versus released/discarded)
#' as a function of size, varying by fleet, sex, age, year, and season.
#'
#' @param RetenFile Data frame from read.table() of RETENSPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#'
#' @return List containing retention specifications:
#' \itemize{
#'   \item NretPatterns - Number of retention patterns defined
#'   \item RetSpec - Matrix of retention pattern specifications (pattern × 5)
#'     \itemize{
#'       \item Column 1: Pattern ID
#'       \item Column 2: Type (2,3,4,6,7,8 for different functional forms)
#'       \item Column 3: First parameter specification
#'       \item Column 4: Number of parameters
#'       \item Column 5: Additional specification
#'     }
#'   \item RetPnt - Array of retention pattern pointers (sex × age × fleet × year × step)
#'     indicating which retention pattern applies for each combination
#'   \item NretPars - Total number of retention parameters to estimate
#'   \item RetenFI - Matrix of fixed retention-at-length patterns (pattern × length)
#'   \item NfixedReten - Number of fixed retention patterns
#' }
#'
#' @details
#' Retention patterns define the probability that a captured animal is kept rather
#' than discarded or released. This is distinct from selectivity (probability of
#' capture given encounter) and legal size regulations (regulatory retention rules).
#'
#' Retention pattern types include:
#' \itemize{
#'   \item Type 2: General parametric form (number of parameters specified in column 4)
#'   \item Type 3: Simple form with 2 parameters
#'   \item Type 4: Single parameter form
#'   \item Type 6: Three-parameter form
#'   \item Type 7: Two-parameter form
#'   \item Type 8: Single parameter form
#' }
#'
#' Fixed retention patterns provide length-specific retention probabilities that
#' are held constant (not estimated). Retention pointers link each combination of
#' sex, age, fleet, year, and season to a specific retention pattern, allowing
#' flexible representation of discard practices across the fishery.
#'
#' The total number of parameters to estimate is calculated by summing across all
#' patterns based on their types. All specifications are written to Echo.out.
#'
#' @keywords internal
ReadRetenFile <- function(RetenFile,GeneralSpecs)
{
  print("READ IN THE RETENTION FILE")
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


#' Parse Recruitment Specifications from RECRUITSPEC.DAT
#'
#' Internal function to read and parse recruitment patterns from RECRUITSPEC.DAT.
#' Defines recruitment allocation across sexes and areas, size distribution at
#' recruitment, and bias adjustment parameters for recruitment deviations.
#'
#' @param RecruitFile Data frame from read.table() of RECRUITSPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#'
#' @return List containing recruitment specifications:
#' \itemize{
#'   \item NrecruitPatternsA - Number of sex/area allocation patterns
#'   \item RecruitSpecsA - Matrix of sex/area allocation specifications (pattern × 3)
#'     \itemize{
#'       \item Column 1: Pattern ID
#'       \item Column 2: Allocation type (0=sex split + areas, 1=sex × area)
#'       \item Column 3: Additional specification
#'     }
#'   \item NrecruitPatternsB - Number of length allocation patterns
#'   \item RecruitSpecsB - Matrix of length allocation specifications (pattern × `[1+2×Nsex]`)
#'   \item RecruitPnt - Matrix of recruitment pattern pointers by year and step (year × step)
#'   \item RecruitLenPnt - Vector of length pattern pointers by area (length Narea)
#'   \item NrecruitPars - Total number of recruitment parameters to estimate
#'   \item RecruitFrac - Matrix of fixed recruitment size distributions (pattern × length)
#'   \item NfixedRecruits - Number of pre-specified recruitment patterns
#'   \item CalcRecruitFrac - Flag indicating whether to calculate (1) or use fixed (0) recruitment fractions
#'   \item Bias_Ramp_Yr1 - First year of bias ramp (relative to Year1)
#'   \item Bias_Ramp_Yr2 - Second year of bias ramp (relative to Year1)
#'   \item Bias_Ramp_Yr3 - Third year of bias ramp (relative to Year1)
#'   \item Bias_Ramp_Yr4 - Fourth year of bias ramp (relative to Year1)
#' }
#'
#' @details
#' This function handles three main components of recruitment specification:
#'
#' **Sex and Area Allocation (Pattern A)**: Defines how total recruitment is
#' distributed across sexes and spatial areas. Two allocation types are supported:
#' \itemize{
#'   \item Type 0: Estimates sex ratio plus independent area proportions (Narea parameters)
#'   \item Type 1: Estimates sex × area proportions independently (Nsex × `[Narea-1]` parameters)
#' }
#'
#' **Length Distribution (Pattern B)**: Defines the size distribution of recruits
#' entering the population. Can be estimated or fixed for each sex. The number of
#' parameters depends on the distributional form specified.
#'
#' **Bias Adjustment Ramp**: Four years defining the "bias ramp" for adjusting
#' recruitment deviation variance. This accounts for retrospective patterns in
#' recruitment estimation, with full variance applied between years 2-3 and
#' ramping at the edges.
#'
#' Pointers link years/steps to specific recruitment patterns, allowing temporal
#' variation in recruitment processes. If CalcRecruitFrac = 0, pre-specified
#' length distributions (RecruitFrac) are used instead of estimating them.
#'
#' All specifications are written to Echo.out for model verification.
#'
#' @keywords internal
ReadRecruitFile <- function(RecruitFile,GeneralSpecs)
{
  print("READ IN THE RECRUIT FILE")
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

  RecruitSpecsB<-matrix (0,nrow=NrecruitPatternsB,ncol=1+2+GeneralSpecs$Nsex)
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Allocate_length")+1;
  for (Ipat in 1:NrecruitPatternsB) RecruitSpecsB[Ipat,] <- as.numeric(RecruitFile[Index+Ipat,1:(1+2+GeneralSpecs$Nsex)])
  write("Specifications for recuitment",EchoFile,append=T)
  write(t(RecruitSpecsB),EchoFile,append=T,ncol=1+2*GeneralSpecs$Nsex)

  NrecruitPars = 0
  for (IrecPat in 1:NrecruitPatternsA)   {
    if (RecruitSpecsA[IrecPat,2]==0) NrecruitPars <- NrecruitPars + (GeneralSpecs$Narea)   # Sex split + Nareas-1
    if (RecruitSpecsA[IrecPat,2]==1) NrecruitPars <- NrecruitPars + (GeneralSpecs$Nsex*(GeneralSpecs$Narea-1))
  }
  for (Ipat in 1:NrecruitPatternsB)
    for (Isex in 1:GeneralSpecs$Nsex)
      if (RecruitSpecsB[Ipat,3+Isex]>=0) NrecruitPars <- NrecruitPars+RecruitSpecs[Ipat,3+GeneralSpecs$Nsex+Isex]
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

  # REcruitment Parameters linking conditions
  Index <- MatchTable(RecruitFile,Char1="#",Char2="Recuitment1");
  npars <- GeneralSpecs$Narea+(NfixedRecruits*2)
  RecparsLink <- as.numeric(RecruitFile[(Index+2):(Index+npars+1),5])
  RecparsPrior <- apply(as.matrix(RecruitFile[(Index+2):(Index+npars+1),6:8]),2,as.numeric)

  Index <- MatchTable(RecruitFile,Char1="#",Char2="Bias",Char3="ramp")+1;
  Bias_Ramp_Yr1 <- as.numeric(RecruitFile[Index,1])-GeneralSpecs$Year1;
  Bias_Ramp_Yr2 <- as.numeric(RecruitFile[Index,2])-GeneralSpecs$Year1;
  Bias_Ramp_Yr3 <- as.numeric(RecruitFile[Index,3])-GeneralSpecs$Year1;
  Bias_Ramp_Yr4 <- as.numeric(RecruitFile[Index,4])-GeneralSpecs$Year1;


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
  ReturnObj$RecparsLink <- RecparsLink
  ReturnObj$RecparsPrior <- RecparsPrior
  ReturnObj$Bias_Ramp_Yr1  <- Bias_Ramp_Yr1;
  ReturnObj$Bias_Ramp_Yr2  <- Bias_Ramp_Yr2;
  ReturnObj$Bias_Ramp_Yr3  <- Bias_Ramp_Yr3;
  ReturnObj$Bias_Ramp_Yr4  <- Bias_Ramp_Yr4;
  return(ReturnObj)
}


#' Parse Growth Specifications from GROWTHSPEC.DAT
#'
#' Internal function to read and parse growth patterns and size-transition matrices
#' from GROWTHSPEC.DAT. Defines how individuals transition between size classes
#' (molting patterns) as a function of area, sex, age, year, and season.
#'
#' @param GrowthFile Data frame from read.table() of GROWTHSPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#'
#' @return List containing growth specifications:
#' \itemize{
#'   \item NgrowthPatterns - Number of growth patterns defined
#'   \item GrowthSpecs - Matrix of growth pattern specifications (pattern × 6)
#'     \itemize{
#'       \item Column 1: Pattern ID
#'       \item Column 2: Growth type
#'       \item Columns 3-6: Pattern-specific specifications
#'     }
#'   \item GrowthPnt - Array of growth pattern pointers (area × sex × age × year × step)
#'     indicating which growth pattern applies for each combination
#'   \item TransInp - Array of pre-specified size-transition matrices (pattern × length × length)
#'   \item NfixedGrowth - Number of pre-specified growth/transition matrices
#'   \item NfixedGrowthSex - Vector indicating which sex each fixed pattern applies to
#'   \item NgrowthPars - Total number of growth parameters to estimate (currently 0)
#' }
#'
#' @details
#' Growth in size-structured models is represented through size-transition matrices
#' that define the probability of transitioning from one size class to another during
#' a molt event. This function handles both parametric growth patterns (which would
#' be estimated) and pre-specified transition matrices.
#'
#' **Pre-specified transition matrices** (TransInp) are length × length matrices where
#' element `[i,j]` represents the probability of an individual in size class i molting
#' to size class j. Each row should sum to 1.0. These matrices can differ by sex,
#' allowing sex-specific growth patterns.
#'
#' **Growth pointers** link each combination of area, sex, age, year, and season to
#' a specific growth pattern, providing flexibility to represent spatial and temporal
#' variation in growth rates (e.g., different growth in warm vs. cool waters, or
#' seasonal growth patterns).
#'
#' Currently, the model uses pre-specified transition matrices exclusively
#' (NgrowthPars = 0), but the structure supports future implementation of estimated
#' parametric growth functions.
#'
#' All specifications are written to Echo.out for model verification.
#'
#' @keywords internal
ReadGrowthFile <- function(GrowthFile,GeneralSpecs)
{
  print("READ IN THE GROWTH FILE")
  Index <- MatchTable(GrowthFile,Char1="#",Char2="Number",Char3="of",Char4="growth")+1;
  NgrowthPatterns <- as.numeric(GrowthFile[Index,1])
  write(paste("Number of growth patterns",NgrowthPatterns),EchoFile,append=T)

  GrowthSpecs<-matrix (0,nrow=NgrowthPatterns,ncol=6)
  for (Ipat in 1:NgrowthPatterns)
    GrowthSpecs[Ipat,]<- as.numeric(GrowthFile[Index+1+Ipat,1:6])
  write("Specifications for growth",EchoFile,append=T)
  write(t(GrowthSpecs),EchoFile,append=T,ncol=6)

  Index <- MatchTable(GrowthFile,Char1="#",Char2="Growth",Char3="parameters")+2;
  ## CHeck if there are growthpars
  NgrowthPars <- ifelse(sum(is.na(as.numeric(GrowthFile[Index:(Index+7),1])))==0, NgrowthPatterns*8, 0)
  write(paste("Number of growth parameters",NgrowthPars),EchoFile,append=T,ncol=3+GeneralSpecs$Nsex)

  # Growth parameters linking conditions
  if(NgrowthPars>0) {
    GrowparsLink <- as.numeric(GrowthFile[(Index):(Index+NgrowthPars-1),5])
    GrowparsPrior <- apply(as.matrix(GrowthFile[(Index):(Index+NgrowthPars-1),6:8]),2,as.numeric)
  } else {GrowparsLink <- 0; GrowparsPrior <- c(0,0)}

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

  write("READ IN THE GROWTH FILE\n\n",EchoFile,append=T)
  ReturnObj <- NULL
  ReturnObj$NgrowthPatterns <- NgrowthPatterns
  ReturnObj$GrowthSpecs <- GrowthSpecs
  ReturnObj$GrowthPnt <- GrowthPnt
  ReturnObj$GrowparsLink <- GrowparsLink
  ReturnObj$GrowparsPrior <- GrowparsPrior
  ReturnObj$TransInp <-TransInp
  ReturnObj$NfixedGrowth <- NfixedGrowth
  ReturnObj$NfixedGrowthSex <- NfixedGrowthSex
  ReturnObj$NgrowthPars <- NgrowthPars
  return(ReturnObj)
}


#' Parse Projection Specifications from PROJSPEC.DAT
#'
#' Internal function to read and parse fishery projection specifications from
#' PROJSPEC.DAT. Defines selectivity, retention, legal size regulations, and
#' discard mortality for future projection years beyond the assessment period.
#'
#' @param ProjFile Data frame from read.table() of PROJSPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#' @param Phi1 Array of discard mortality rates from assessment period to be
#'   extended into projection years (fleet × age × year × step)
#' @param Catch1 Array of observed catch from ReadDataFile() (year × step ×
#'   fleet), already dimensioned `Nyear+MaxProjYr` on the year axis with the
#'   projection-year rows at zero. Overlaid with the projected catch column
#'   from PROJECTIONS.DAT.
#'
#' @return List containing projection specifications:
#' \itemize{
#'   \item Nproj - Number of projection years
#'   \item SelPntFut - Array of future selectivity pattern pointers
#'     (sex × age × fleet × projection_year × step)
#'   \item RetPntFut - Array of future retention pattern pointers
#'     (sex × age × fleet × projection_year × step)
#'   \item LegalFleetPntFut - Array of future legal size pattern pointers
#'     (sex × age × fleet × projection_year × step)
#'   \item Phi - Extended discard mortality array including projection years
#'     (fleet × age × `[assessment_years + projection_years]` × step)
#'   \item ProjType - 1 = catch-based projection, 2 = harvest-rate/effort-based
#'   \item Catch - Catch array (year × step × fleet) with the projection-year
#'     rows filled in from the "catch" column of PROJECTIONS.DAT (unchanged
#'     from \code{Catch1} elsewhere)
#'   \item ProjHarvestRate - Array of imposed harvest rates
#'     (projection_year × step × fleet), filled in from the "Hrate" column of
#'     PROJECTIONS.DAT. Populated regardless of ProjType (harmless when
#'     ProjType == 1, since the compiled model then never reads it) so the
#'     scenario can be switched with a one-line change to PROJECTIONS.DAT
#'     without re-deriving either schedule.
#' }
#'
#' @details
#' This function extends the stock assessment model into future projection years
#' for forecasting and management strategy evaluation. It specifies how fishery
#' characteristics will be configured in future years:
#'
#' **Selectivity patterns** define future gear selectivity, typically assumed
#' constant at recent levels or following specified management scenarios.
#'
#' **Retention patterns** define future discard practices, which may change
#' under different management strategies or behavioral assumptions.
#'
#' **Legal size patterns** specify future size limit regulations that may be
#' modified as part of management scenarios being evaluated.
#'
#' **Discard mortality** (Phi) extends the survival rate of discarded animals
#' into projection years. The function takes the existing Phi array from the
#' assessment period and populates projection year values based on specifications
#' in the projection file.
#'
#' All pointer arrays have dimensions matching the assessment period structure
#' but indexed to projection years. Calendar years in projections are calculated
#' as: Nyear + Year1 + projection_year - 1.
#'
#' **Projection type and schedule** ("# Specifications for projections
#' (1=Catch;2=Harvestrate)" followed by a data block of
#' Year/Step/Fleet/Catch/Hrate rows). Both the catch and harvest-rate columns
#' are always parsed and stored (into Catch and ProjHarvestRate respectively)
#' regardless of ProjType, so which scenario actually runs is controlled
#' entirely by the ProjType flag passed through to the compiled model -- e.g.
#' switching ProjType from 1 to 2 to compare a catch-based projection against
#' a harvest-rate-based one doesn't require re-deriving either schedule.
#'
#' All specifications are written to Echo.out for verification. If Nproj = 0,
#' no projections are performed and arrays remain at default values.
#'
#' @keywords internal
ReadProjFile <- function(ProjFile, GeneralSpecs, Phi1, Catch1, Echo = TRUE)
{
  print("READ IN THE PROJECTION FILE")
  Index <- MatchTable(ProjFile,Char1="#",Char2="Number",Char3="of",Char4="projection")+1;
  Nproj <- as.numeric(ProjFile[Index,1])

  Index <- MatchTable(ProjFile,Char1="#",Char2="Specifications",Char3="for",Char4="gear")+2;
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
  if (Echo) {
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
  }

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
  if (Echo) {
    write("Projected Specifications for retention pointers",EchoFile,append=T)
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
  }

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
  if (Echo) {
    write("Projected Specifications for fleet legal pointers",EchoFile,append=T)
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
              OutM[Ipnt,(5:(4+GeneralSpecs$Nstep))] <- LegalFleetPntFut[Isex,Iage,Ifleet,Iyear,]
            }
    write(t(OutM),EchoFile,append=T,ncol=4+GeneralSpecs$Nstep)
  }

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
          testindex(ProjFile,pos,1,Iage)
          testindex(ProjFile,pos,2,Ifleet)
          testindex(ProjFile,pos,3,Istep)
          for (Iyear in 1:Nproj) Phi[Ifleet,Iage,Iyear+GeneralSpecs$Nyear,Istep] <- as.numeric(ProjFile[Index+Ipnt,3+Iyear])
        }
  if (Echo) {
    write("discard mortality",EchoFile,append=T)
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
  }

  # ── Projection type + catch / harvest-rate schedule ────────────────────
  Index <- MatchTable(ProjFile,Char1="#",Char2="Specifications",Char3="for",Char4="projections")+1
  ProjType <- as.numeric(ProjFile[Index,1])
  if (!ProjType %in% c(1,2))
    stop("ReadProjFile: unrecognised ProjType ", ProjType,
         " (expected 1 = catch or 2 = harvest rate/effort).", call. = FALSE)

  Catch <- Catch1
  ProjHarvestRate <- array(0,dim=c(GeneralSpecs$MaxProjYr,GeneralSpecs$Nstep,GeneralSpecs$Nfleet))

  if (Nproj > 0)
  {
    Index <- MatchTable(ProjFile,Char1="#",Char3="data")
    Nobs  <- as.numeric(ProjFile[Index+1,1])
    Index <- Index + 2

    if (Nobs > 0)
      for (Iobs in 1:Nobs)
      {
        YearRow  <- as.numeric(ProjFile[Index+Iobs,1]) - GeneralSpecs$Year1 + 1
        Step     <- as.numeric(ProjFile[Index+Iobs,2])
        Fleet    <- as.numeric(ProjFile[Index+Iobs,3])
        CatchVal <- as.numeric(ProjFile[Index+Iobs,4])
        HrateVal <- as.numeric(ProjFile[Index+Iobs,5])

        Catch[YearRow,Step,Fleet] <- CatchVal
        ProjHarvestRate[YearRow-GeneralSpecs$Nyear,Step,Fleet] <- HrateVal
      }

    if (Echo) {
      write(paste("Projection type (1=catch, 2=harvest rate):", ProjType),EchoFile,append=T)
      write("Projected catch / harvest rate schedule (Year Step Fleet Catch Hrate)",EchoFile,append=T)
      if (Nobs > 0) write(t(as.matrix(ProjFile[(Index+1):(Index+Nobs),1:5])),EchoFile,append=T,ncolumns=5)
    }
  }

  if (Echo) write("READ IN THE PROJECTION FILE\n\n",EchoFile,append=T)

  ReturnObj <- NULL
  ReturnObj$Nproj <- Nproj
  ReturnObj$SelPntFut <- SelPntFut
  ReturnObj$RetPntFut <- RetPntFut
  ReturnObj$LegalFleetPntFut <- LegalFleetPntFut
  ReturnObj$Phi <- Phi
  ReturnObj$ProjType <- ProjType
  ReturnObj$Catch <- Catch
  ReturnObj$ProjHarvestRate <- ProjHarvestRate

  return(ReturnObj)
}

#' Read Initial Parameter Values, Bounds, and Estimation Phases
#'
#' Internal function to extract starting values, bounds, and estimation phases for
#' all model parameters from control and specification files. Organizes parameters
#' by type and validates that all values are numeric (no missing values).
#'
#' @param ControlFile Data frame from read.table() of CONTROL.DAT
#' @param SelexFile Data frame from read.table() of SELEXSPEC.DAT
#' @param RetainFile Data frame from read.table() of RETENSPEC.DAT
#' @param RecruitFile Data frame from read.table() of RECRUITSPEC.DAT
#' @param GrowthFile Data frame from read.table() of GROWTHSPEC.DAT
#' @param MoveFile Data frame from read.table() of MOVESPEC.DAT
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#' @param ControlSpecs List from ReadControlFile() containing control specifications
#' @param SelexSpecs List from ReadSelexFile() containing selectivity specifications
#' @param RetenSpecs List from ReadRetenFile() containing retention specifications
#' @param GrowthSpecs List from ReadGrowthFile() containing growth specifications
#' @param MoveSpecs List from ReadMoveFile() containing movement specifications
#'
#' @return Nested list containing parameter specifications organized by type:
#' \itemize{
#'   \item MainPars - Unfished recruitment (R0), natural mortality (M), mortality offsets,
#'     scale factors, recruitment variability (sigmaR), and initial abundances by area
#'     \itemize{
#'       \item Initial - Vector of starting values
#'       \item Bnd - Matrix of bounds (parameter × 2: lower, upper)
#'       \item Phase - Vector of estimation phases (negative = fixed)
#'     }
#'   \item RecruitPars - Recruitment distribution parameters (sex ratio, area allocation, size distribution)
#'   \item PuerPowPars - Puerulus settlement power parameters
#'   \item SelPars - Selectivity parameters
#'   \item RetPars - Retention parameters
#'   \item RecDevs - Annual recruitment deviations
#'   \item RecSpatDevs - Spatial recruitment deviations
#'   \item Qpars - Catchability coefficients
#'   \item efpars - Fishing efficiency/power parameters
#'   \item InitPars - Initial population size-structure parameters
#'   \item MovePars - Movement/migration parameters
#'   \item GrowthPars - Growth parameters
#' }
#'
#' @details
#' This function consolidates all model parameters from multiple specification files
#' into a single organized structure for model initialization. For each parameter type,
#' it extracts three components:
#'
#' **Initial values**: Starting values for optimization or fixed values if not estimated
#'
#' **Bounds**: Lower and upper limits constraining parameter space during estimation
#'
#' **Phase**: Estimation phase controlling when parameters are estimated. Negative
#' phases indicate fixed parameters. Parameters are estimated in order of ascending
#' phase number (e.g., phase 1 parameters before phase 2).
#'
#' The function handles several special initialization options (InitOpt) that determine
#' how the initial population structure is specified:
#' \itemize{
#'   \item InitOpt 0: Unfished equilibrium (no initial size parameters estimated)
#'   \item InitOpt 1,3,5: Full age-length structure (area × age × length parameters)
#'   \item InitOpt 2,4: Length structure only (area × length parameters)
#' }
#'
#' Phase assignments are automatically adjusted based on InitOpt settings to ensure
#' appropriate parameters are fixed or estimated.
#'
#' The function validates all parameter values to ensure they are numeric and reports
#' any NA values that would cause model failure. Pre-specified recruitment deviations
#' can be provided rather than estimated.
#'
#' All parameter specifications are written to Echo.out for verification.
#'
#' @keywords internal
ReadInitialValues <- function(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,
                              GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
{
  # Main parameters
  # R0, M-bar, M-at-age-offset, WjotyesScaleM, RedsScaleQ, SigmaR
  OK <- 1
  NmainPars = 4+(GeneralSpecs$Nage)+GeneralSpecs$Narea;                       #// 5 is no virgin M
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
  if(is.na(sum(MainPars))) { warning("\nThere are NA's in Main Pars\n", call. = FALSE); OK <- 0   }

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
  if(is.na(sum(RecDevPars))) { warning("\nThere are NA's in Recruitment deviations Pars\n", call. = FALSE); OK <- 0   }

  Index <- MatchTable(ControlFile,Char1="#",Char2="Prespecify_spatial_rec_devs");
  PreSpecifySpatRecDevs <- as.numeric(ControlFile[Index+1,1])
  RecSpatDevPars <- rep(0,NrecSpatDev)
  if (PreSpecifySpatRecDevs == 1){
    for (Dyr in 1: NrecSpatDev) RecSpatDevPars[Dyr] <-  as.numeric(ControlFile[Index+1+Dyr,1])}
  if (NrecSpatDev==0) { NrecSpatDev <- 1; RecSpatDevPars = 0; RecSpatDevBnd <- matrix(c(-15,15),ncol=2,nrow=2); RecSpatDevPhase <- -1; }
  if(is.na(sum(PreSpecifySpatRecDevs))) { warning("\nThere are NA's in Spatial Recruitment Pars\n", call. = FALSE); OK <- 0   }

  Index <- MatchTable(SelexFile,Char1="#",Char2="Selectivity",Char3="Parameters")+1;
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
  if(is.na(sum(SelPars))) { warning("\nThere are NA's in Selectivity Pars\n", call. = FALSE); OK <- 0   }

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
  if(is.na(sum(RetPars))) { warning("\nThere are NA's in Retention Pars\n", call. = FALSE); OK <- 0   }

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
  if(is.na(sum(RecruitPars))) { warning("\nThere are NA's in Recruitment Pars\n", call. = FALSE); OK <- 0   }

  Index <- MatchTable(RecruitFile,Char1="#",Char2="Puerulus",Char3="Power")+1;
  NPuerPow <- as.numeric(RecruitFile[Index,1])
  PuerPowPars <- rep(0,NPuerPow)
  PuerPowBnd <- matrix(0,nrow=NPuerPow,ncol=2)
  PuerPowPhase <- rep(-99,NPuerPow)
  if(NPuerPow>0){
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
  }
  if(is.na(sum(PuerPowPars))) { warning("\nThere are NA's in Puerulus Pars\n", call. = FALSE); OK <- 0   }

  Index <- MatchTable(GrowthFile,Char1="#",Char2="Growth",Char3="parameters")+1;
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
  }  else  GrowthPhase[1] <- -100
  if(is.na(sum(GrowthPars))) { warning("\nThere are NA's in Growth Pars\n", call. = FALSE); OK <- 0   }

  Index <- MatchTable(MoveFile,Char1="#",Char2="Movement",Char3="parameters")+1;
  MovePars <- rep(0,MoveSpecs$NmovePars)
  MoveBnd <- matrix(0,nrow=MoveSpecs$NmovePars,ncol=2)
  MovePhase <- rep(-99,MoveSpecs$NmovePars)
  if(MoveSpecs$NmovePars>0){
    for (Ipar in 1:MoveSpecs$NmovePars)
    {
      MovePars[Ipar] <- as.numeric(MoveFile[Index+Ipar,3])
      MoveBnd[Ipar,1] <- as.numeric(MoveFile[Index+Ipar,1])
      MoveBnd[Ipar,2] <- as.numeric(MoveFile[Index+Ipar,2])
      MovePhase[Ipar] <- as.numeric(MoveFile[Index+Ipar,4])
    }
    write("Initial movement parameters",EchoFile,append=T)
    write(MovePars,EchoFile,append=T)}
  if(is.na(sum(MovePars))) { warning("\nThere are NA's in Move Pars\n", call. = FALSE); OK <- 0   }

  NInitPar <- 1;

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
  if(is.na(sum(QPars))) { warning("\nThere are NA's in Q Pars\n", call. = FALSE); OK <- 0   }

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
  if(is.na(sum(efPars))) { warning("\nThere are NA's in Efficiency Pars\n", call. = FALSE); OK <- 0   }

  if(OK==1) { message("All initial parameters are provided")}
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
 ReturnObj$MovePars$Initial <- MovePars
 ReturnObj$MovePars$Bnd <- MoveBnd
 ReturnObj$MovePars$Phase <- MovePhase
 ReturnObj$GrowthPars$Initial <- GrowthPars
 ReturnObj$GrowthPars$Bnd <- GrowthBnd
 ReturnObj$GrowthPars$Phase <- GrowthPhase
 return(ReturnObj)

}

#' Identify Parameters to be Estimated and Estimation Phases
#'
#' Internal function that analyzes the parameter structure to determine which
#' parameters will be estimated versus held fixed, and organizes them by
#' estimation phase. Creates a summary table for model documentation.
#'
#' @param InitialVars Nested list from ReadInitialValues() containing parameter
#'   specifications with Initial, Bnd, and Phase components for each parameter type
#'
#' @return Invisibly returns NULL. Side effects include:
#' \itemize{
#'   \item Creates global variable MaxPhase - highest estimation phase number
#'   \item Writes 'Output/Parameters_solved.txt' - table summarizing estimated parameters
#'   \item Prints summary table to console
#' }
#'
#' @details
#' This function processes the parameter structure created by ReadInitialValues()
#' to identify which parameters will be actively estimated during model fitting.
#' Parameters with positive phase values are estimated, while those with negative
#' or zero phases remain fixed.
#'
#' **Estimation phases** control the order and grouping of parameter estimation:
#' \itemize{
#'   \item Phase < 0: Parameter held fixed at initial value
#'   \item Phase 1: Estimated first (typically core parameters like R0, M)
#'   \item Phase 2+: Estimated after earlier phases converge
#' }
#'
#' The phased estimation approach improves optimization by:
#' \itemize{
#'   \item Allowing core parameters to stabilize before adding complexity
#'   \item Reducing correlations between parameter groups
#'   \item Improving convergence reliability
#'   \item Making it easier to diagnose estimation problems
#' }
#'
#' The function creates a summary table with three columns:
#' \itemize{
#'   \item Parameter - Parameter type name (e.g., "MainPars", "RecDevs")
#'   \item Number - Count or indices of parameters to be estimated
#'   \item Phases - Space-separated list of phase numbers for this parameter type
#' }
#'
#' For MainPars specifically, the Number column shows the actual parameter indices
#' that will be estimated (to show which specific main parameters are active). For
#' other parameter types, it defaults to 1 to indicate the parameter type is active.
#'
#' The global variable MaxPhase is set to support phase-limited optimization runs
#' and progress monitoring during model fitting.
#'
#' @examples
#' \dontrun{
#' InitVals <- ReadInitialValues(...)
#' Parssolved(InitVals)
#' # Creates Output/Parameters_solved.txt and sets MaxPhase
#' }
#'
#' @keywords internal
Parssolved <- function(InitialVars){
  suppressWarnings(rm(parsolve, pos=1 ))
  MaxPhase <<- 1
  for(n in 1:length(unique(names(InitialVars)))){
    par <- unique(names(InitialVars))[n]
    tmp_n <- InitialVars[par][[1]]$Phase
    tmp_nsum <- paste(unique(tmp_n[tmp_n>0]),collapse = ' ')
    pos <- which(tmp_n>0)
    npos <- length(pos)
    if(length(tmp_n)>0) if(max(tmp_n)>MaxPhase) MaxPhase <<- max(tmp_n)
    if(length(pos)>0) {
      if(par!="MainPars") pos <-1
      if(!exists('parsolve')) { parsolve <- data.frame(Parameter=par,Number=npos,Phases=tmp_nsum)} else parsolve <- rbind(parsolve,data.frame(Parameter=par,Number=npos,Phases=tmp_nsum))
    }
  }
  if(!exists('parsolve')) {parsolve <- NA}
  write.table(parsolve,'Output/Parameters_solved.txt',quote = F, sep='\t',row.names = F)
}

#' Copy Data to Clipboard for Excel Pasting
#'
#' Convenience function to copy R objects (data frames, matrices, vectors) to the
#' system clipboard in tab-delimited format suitable for pasting directly into
#' Excel or other spreadsheet applications.
#'
#' @param x Data object to copy (data.frame, matrix, or vector)
#' @param rnames Logical indicating whether to include row names in output.
#'   Default is FALSE
#'
#' @return Invisibly returns NULL. Side effect is copying data to clipboard.
#'
#' @details
#' This is a utility function for interactive data exploration and presentation.
#' It writes the data to the system clipboard using tab separators, which Excel
#' and similar programs recognize for proper column alignment.
#'
#' After running this function, you can paste (Ctrl+V or Cmd+V) the data directly
#' into Excel, maintaining the row and column structure.
#'
#' **Note**: This function uses the "clipboard" connection which is primarily
#' supported on Windows systems. On Mac/Linux, you may need to use alternative
#' approaches like pbcopy/pbpaste or xclip.
#'
#' @examples
#' \dontrun{
#' # Copy a data frame to clipboard
#' results <- data.frame(Year = 2020:2024, Catch = c(1000, 1200, 1100, 1300, 1250))
#' toXL(results)
#' # Now paste into Excel
#'
#' # Include row names
#' toXL(mtcars[1:5, 1:3], rnames = TRUE)
#' }
#'
#' @export
toXL <- function(x, rnames=FALSE){
  write.table(x, "clipboard", sep="\t", row.names = rnames)
}

#' Load and Parse All IMuLT Model Input Files
#'
#' Reads all model specification files (.DAT files) from the current directory,
#' parses them into structured data objects, and prepares them for model estimation.
#' This function must be called after navigating to a model run directory and before
#' loading parameters or running the model.
#'
#' @return NULL. Creates multiple global objects in the parent environment:
#' \itemize{
#'   \item Data - Master list containing all model specifications
#'   \item StarterFile, Starter - Starter file and parsed specifications
#'   \item DataFile, GeneralSpecs, TheData - General model structure and data
#'   \item ControlFile, ControlSpecs - Control parameters and specifications
#'   \item SelexFile, SelexSpecs - Selectivity specifications
#'   \item RetenFile, RetenSpecs - Retention/discard specifications
#'   \item RecruitFile, RecruitSpecs - Recruitment specifications
#'   \item GrowthFile, GrowthSpecs - Growth transition specifications
#'   \item MoveFile, MoveSpecs - Movement/migration specifications
#'   \item ProjFile, ProjectSpecs - Projection specifications
#'   \item EchoFile - Path to echo output file
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads STARTER.DAT to get names of all other input files
#'   \item Loads each specification file (.DAT format)
#'   \item Parses each file into structured R objects
#'   \item Combines all specifications into the master Data list
#'   \item Creates Output/Echo.out for model diagnostics
#' }
#'
#' Required input files in working directory:
#' \itemize{
#'   \item STARTER.DAT - File paths and basic settings
#'   \item DATA.DAT - Catch, CPUE, and biological data
#'   \item CONTROL.DAT - Model parameters and specifications
#'   \item SELEXSPEC.DAT - Selectivity specifications
#'   \item RETAINSPEC.DAT - Retention specifications
#'   \item RECRUITSPEC.DAT - Recruitment specifications
#'   \item GROWTHSPEC.DAT - Growth specifications
#'   \item MOVESPEC.DAT - Movement specifications
#'   \item PROJECTIONS.DAT - Projection specifications
#' }
#'
#' @note
#' \itemize{
#'   \item Must be run from within a model run directory containing all .DAT files
#'   \item Creates global variables - use in interactive sessions or scripts
#'   \item An Output/ directory must exist (or will be created) for Echo.out
#' }
#'
#' @examples
#' \dontrun{
#' # Standard workflow
#' choose_model()  # Navigate to model directory
#' LoadData()      # Load all input files
#' LoadPars()      # Load parameter specifications
#' SolveModelNew() # Run model estimation
#' }
#'
#' @seealso
#' \code{\link{choose_model}} for selecting model directory,
#' \code{\link{LoadPars}} for loading parameters after data,
#' \code{\link{BuildInputFiles}} for creating input files
#'
#' @export
LoadData <- function() {
  # read in files
  StarterFile <<- read.table("Starter.dat",comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  Starter <<- ReadStarterFile(StarterFile)
  DataFile <<- read.table(Starter$DataFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  ControlFile <<- read.table(Starter$ControlFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  SelexFile <<- read.table(Starter$SelexFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  RetenFile <<- read.table(Starter$RetainFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  RecruitFile <<- read.table(Starter$RecruitFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  ReprodFile <<- read.table(Starter$ReproFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  GrowthFile <<- read.table(Starter$GrowthFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  MoveFile <<- read.table(Starter$MoveFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  TagFile <<- read.table(Starter$TagFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)
  PropnFile <<- read.table(Starter$PropFFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:100)
  ProjFile <<- read.table(Starter$ProjectionsFileName,comment.char = "?",fill=T,blank.lines.skip=T,stringsAsFactors=F,col.names=1:200)

  EchoFile <<- "Output/Echo.out"
  write(paste(rep(".", 300),collapse=' ') ,EchoFile)
  write("The is an echo file",EchoFile)

  # Read in the data file
  GeneralSpecs <<- ReadGeneralFile(DataFile)
  Data <<- GeneralSpecs

  # Read in the data file
  TheData <<- ReadDataFile(DataFile,GeneralSpecs)
  Data <<- append(Data,TheData)
  GeneralSpecs$NQpars <<- Data$NQparPass

  # # Read in the tag file
  TagSpecs <<- ReadTagFile(TagFile,PropnFile,GeneralSpecs,TheData)
  Data <<- append(Data,TagSpecs)

  # Read in the Control file
  ControlSpecs <<- ReadControlFile(ControlFile,GeneralSpecs,TheData)
  Data <<- append(Data,ControlSpecs)

  # Read in the Reproduction file
  ReproSpecs <<- ReadReprodFile(ReprodFile,GeneralSpecs,TheData)
  Data <<- append(Data,ReproSpecs)

    # Read in the projections file
  ProjectSpecs <<- ReadProjFile(ProjFile,GeneralSpecs,Data$Phi1,Data$Catch)
  Data$Catch <<- ProjectSpecs$Catch    # overlay projection-year catch (ProjType==1)
  ProjectSpecs$Catch <- NULL           # avoid a duplicate 'Catch' entry on append
  Data <<- append(Data,ProjectSpecs)

  # Read in the Selectivity file
  SelexSpecs <<- ReadSelexFile(SelexFile,GeneralSpecs)
  Data <<- append(Data,SelexSpecs)

  # Read in the Retention file
  RetenSpecs <<-ReadRetenFile(RetenFile,GeneralSpecs)
  Data <<- append(Data,RetenSpecs)

  # Read in the Recruiment file
  RecruitSpecs <<- ReadRecruitFile(RecruitFile,GeneralSpecs)
  Data <<- append(Data,RecruitSpecs)

  # Read in the Growth file
  GrowthSpecs <<-ReadGrowthFile(GrowthFile,GeneralSpecs)
  Data <<- append(Data,GrowthSpecs)

  # Read in the Movement file
  MoveSpecs <<- ReadMoveFile(MoveFile,GeneralSpecs)
  Data <<- append(Data,MoveSpecs)

  # preset some things
  FullOutput <<- FALSE
  Data$DoProject <<- 0

  outtmp <- isnafunc2()
  if(!is.null(outtmp[[2]]))   { message("\nSome data objects are empty (which can be OK): ", paste(outtmp[[2]], collapse = ', '), '\n') }
}

