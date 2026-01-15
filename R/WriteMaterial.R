
#' Write Comprehensive Model Output Files
#'
#' Master output function to write all IMuLT model results to standardized output
#' files. Creates Output.RL (main results), SDReport.RL (parameter estimates and
#' standard errors), and Model.ppp (parameter values) for use in diagnostics,
#' plotting, and subsequent analyses.
#'
#' @param Report List of model results from objective function evaluation, containing
#'   all predicted values, likelihoods, derived quantities, and population dynamics
#' @param SDrep Data frame of standard error report from RTMB/TMB sdreport(), containing
#'   parameter estimates and standard errors for reported variables
#' @param fullrep Complete standard error report including all parameters and derived
#'   quantities with estimates, standard errors, and names
#' @param pin Named list of initial parameter values (input to optimization)
#' @param pout Named list of final parameter values (output from optimization)
#' @param GeneralSpecs List from ReadGeneralFile() containing model dimensions
#' @param ControlSpecs List from ReadControlFile() containing control specifications
#' @param TheData List from ReadDataFile() containing all model data
#' @param CurrPhase Integer indicating estimation phase (default = 2)
#' @param best Numeric vector of best parameter values from optimization (default = -1)
#' @param grad Numeric vector of final gradients for all parameters (default = -1)
#'
#' @return Invisibly returns NULL. Side effects include creating output files:
#' \itemize{
#'   \item Output.RL - Main results file with all model outputs
#'   \item SDReport.RL - Standard error report (copy of fullrep)
#'   \item Model.ppp - Parameter file with estimated and fixed values
#' }
#'
#' @details
#' This is the comprehensive output writer for IMuLT. It organizes all model results
#' into a standardized format for subsequent analysis, plotting, and reporting. The
#' function automatically handles missing standard errors (when sdreport not run).
#'
#' **Output.RL Structure** (main results file):
#'
#' *Likelihood Components*:
#' \itemize{
#'   \item Total objective function value
#'   \item Individual likelihood components (catch, CPUE, numbers, length, larval)
#'   \item Weighted and unweighted likelihoods by data type
#'   \item Penalty terms (initial N, recruitment, recruitment smoothing)
#'   \item Likelihood breakdown by fleet/group
#' }
#'
#' *Parameter Estimates*:
#' \itemize{
#'   \item Complete parameter table with names, indices, estimates, SEs, gradients, bounds
#'   \item Distinction between estimated (positive phase) and fixed parameters
#'   \item Total count of estimated parameters
#' }
#'
#' *Population Dynamics*:
#' \itemize{
#'   \item Egg production (total and by area)
#'   \item Recruitment (total and by area, with SEs when available)
#'   \item Recruitment size distributions
#'   \item Legal biomass by area (various definitions: >76mm, reference selectivity)
#'   \item Mature biomass by sex and area
#'   \item Virgin/unfished biomass
#'   \item Extended projections when applicable
#' }
#'
#' *Mortality and Exploitation*:
#' \itemize{
#'   \item Harvest rates by management zone
#'   \item Fishing mortality (F) by year, time step, and fleet
#'   \item Natural mortality (M) by area, age, and year (including density dependence)
#' }
#'
#' *Model Fits to Data*:
#' \itemize{
#'   \item Catch data: observed, predicted, and fishing mortality by fleet
#'   \item Index data: observations, predictions, residuals with CVs
#'   \item Catch-in-numbers: fits and residuals
#'   \item Length compositions: observed vs predicted with effective sample sizes
#'   \item Larval/puerulus data: observations vs predictions
#'   \item Fishing efficiency trends (technological creep)
#' }
#'
#' *Biological Processes*:
#' \itemize{
#'   \item Selectivity patterns (actualized for each combination of factors)
#'   \item Retention patterns (legal size selectivity)
#'   \item Growth curves by year, area, sex, and age
#'   \item Recruitment distributions (spatial and size-based)
#'   \item Movement patterns by size
#' }
#'
#' *Population Structure*:
#' \itemize{
#'   \item Initial N-matrix: starting population by area/sex/age/length
#'   \item Simplified N-matrix: population summed across ages by year
#'   \item Full N-matrix: complete population array (area × year × step × sex × age × length)
#'   \item F-matrix: fishing mortality array (year × step × fleet)
#' }
#'
#' *Data Diagnostics*:
#' \itemize{
#'   \item Length composition tuning (Francis multipliers for reweighting)
#'   \item Effective sample sizes vs input sample sizes
#'   \item Sigma estimates for each data type
#'   \item Lambda (overdispersion) parameters
#'   \item Q (catchability) estimates by index
#' }
#'
#' **Model.ppp File**: Contains parameter values in simple format for reading by
#' external programs. Lists both estimated (with final values) and fixed (with
#' initial values) parameters with comments indicating parameter names and indices.
#'
#' **Handling Missing Standard Errors**: When standard errors are unavailable
#' (sdreport not run), the function outputs estimates only, adding 'nose' or 'se=NA'
#' flags as appropriate. This allows output generation even when Hessian is singular.
#'
#' **Year Indexing**: The function carefully converts between model year indices
#' (which may start at 1) and calendar years (Year1 through Year2), including
#' burn-in period adjustments where necessary. Larval data years account for both
#' burn-in and settlement delay (Larval_Offset).
#'
#' **Output Directory Management**: If an 'Output' subdirectory exists, the function
#' changes to it before writing files, then returns to the parent directory afterward.
#'
#' **Data Tuning Statistics**: Calculates Francis (1011) weighting multipliers for
#' length composition data by comparing variance of residuals in mean length. Also
#' computes effective sample sizes using McAllister-Ianelli method for comparison
#' with input sample sizes.
#'
#' The function relies on tidyr, dplyr, and magrittr for data manipulation when
#' handling standard error output from RTMB/TMB.
#'
#' @note This function should be called after successful model optimization to
#' preserve results. The output files are read by MakeOutPut() for report generation
#' and by other analysis functions. File formats are specific to IMuLT and must
#' maintain their structure for compatibility with downstream processing.
#'
#' @references
#' Francis, R.I.C.C. (2011). Data weighting in statistical fisheries stock assessment
#' models. Canadian Journal of Fisheries and Aquatic Sciences, 68(6), 1124-1138.
#'
#' @examples
#' \dontrun{
#' # After optimization
#' obj <- MakeADFun(data, parameters, map = final_map)
#' opt <- nlminb(start, obj$fn, obj$gr, lower = lower, upper = upper)
#'
#' # Get standard errors
#' sdrep <- sdreport(obj)
#'
#' # Write outputs
#' Report <- obj$report()
#' WriteOutput(
#'   Report = Report,
#'   SDrep = summary(sdrep, "report"),
#'   fullrep = summary(sdrep, "all"),
#'   pin = parameters,
#'   pout = opt$par,
#'   GeneralSpecs = GeneralSpecs,
#'   ControlSpecs = ControlSpecs,
#'   TheData = Data,
#'   best = opt$par,
#'   grad = obj$gr(opt$par)
#' )
#' }
#'
#' @seealso
#' \code{\link{MakeOutPut}} for generating HTML reports from output files
#' \code{\link{LoadOutputData}} for reading output files
#'
#' @keywords internal
WriteOutput <- function(Report,SDrep,fullrep,pin,pout,GeneralSpecs,ControlSpecs,TheData,CurrPhase=2,best=rep(-1,1000),grad=rep(-1,1000)){

  library(tidyr)
  library(dplyr)
  library(magrittr)
  files <- list.files()
  if(max(files=='Output')==1) { setwd(paste(getwd(), "/Output",sep=""))  }

  write.table(fullrep,'SDReport.RL',sep=' ', quote=F)
  stdrep <- fullrep
  OutputFile <- "Output.RL"
  ParFileName <- "Model.ppp"

  write("# Likelihood summary",OutputFile)
  write("Catch penalty",OutputFile,append=T)
  write(paste("Total objective function",Report$neglogL),OutputFile,append=T)
  write(paste("Catch likelihood",Report$CatchLike),OutputFile,append=T)
  write(paste("Cpue likelihood",Report$CpueLike, " (Weighted ", Report$Weighted_CpueLike, ")"),OutputFile,append=T)
  write(paste("Numbers likelihood",Report$NumbersLike, " (Weighted ", Report$Weighted_NumbersLike, ")"),OutputFile,append=T)
  write(paste("Length likelihood",Report$LengthLike, " (Weighted ", Report$Weighted_LengthLike, ")"),OutputFile,append=T)
  write(paste("Larval likelihood",Report$LarvalLike, " (Weighted ", Report$Weighted_LarvalLike, ")"),OutputFile,append=T)
 # write(paste("Tag1 Likelihood",sum(Report$TagLike1)," (Weighted ", Report$Weighted_TagLike1,")"),OutputFile,append=T)
 # write(paste("Tag2 Likelihood",sum(Report$TagLike2)," (Weighted ", Report$Weighted_TagLike2,")"),OutputFile,append=T)
  write(paste("Initial N penalty",Report$Initial_pen),OutputFile,append=T)
  write(paste("Recruitment penalty",Report$Rec_Penal),OutputFile,append=T)
  write(paste("Recruitment smooth penalty",Report$Rec_Penal_Smooth),OutputFile,append=T)

  write("\n# Likelihood by fleet",OutputFile,append=T)
  write(paste("Cpue likelihood",paste(Report$CpueLikeComp[])),OutputFile,append=T)
  write(paste("Numbers likelihood",Report$NumbersLikeComp),OutputFile,append=T)
  write(paste("Length likeliood",Report$LengthLikeComps),OutputFile,append=T)
  write(paste("Larval likeliood",Report$LarvalLikeComps),OutputFile,append=T)

  write("\n# Likelihood by group",OutputFile,append=T)
 # write(paste("Tag Likelihood 1",Report$TagLike1),OutputFile,append=T)
 # write(paste("Tag Likelihood 2",Report$TagLike2),OutputFile,append=T)


  #print("AEP IS STILL WORKING ON THIS")
  write("\n# Labels",OutputFile,append=T)
  write("#Simplified N-matrix",OutputFile,append=T)
  write("#F-matrix",OutputFile,append=T)
  write("#Recruitment patterns by timestep, area and sex",OutputFile,append=T)
  write("#Recruitment patterns by sex and size",OutputFile,append=T)


  write("\n#Initiation Option",OutputFile,append=T)
  write(Data$InitOpt,OutputFile,append=T)


  write("\n# parameter table",OutputFile,append=T)
  write("# Parameter Par_cnt Estpar_cnt Estimate SD Gradient lwrBound uprBound",OutputFile,append=T)
  ParName <- names(pin)
  Ipnt <- 0; Iqnt <- 0
  write("# parameters",ParFileName)
  write("# dummy",ParFileName,append=T)
  write("0",ParFileName,append=T)
  for (ParName in names(pin))
   {
    if (ParName != "dummy")
     {
      print(ParName)
      ThePar <- InitialVars[[ParName]]
      for (Ipar in 1:length(ThePar$Initial))
       {
        Iqnt <- Iqnt + 1
        write(paste("#",ParName,"_",Ipar," ",Iqnt," ",sep=""),ParFileName,append=T)
        ## Add to the par out file
        if (ThePar$Phase[Ipar] > 0)
         {
          Ipnt <- Ipnt + 1;
          xx <- paste(ParName,"_",Ipar," ",Iqnt," " ,Ipnt," ", stdrep[Ipnt,1]," ",stdrep[Ipnt,2]," ",as.vector(grad)[Iqnt]," ",ThePar$Bnd[Ipar,1]," ",ThePar$Bnd[Ipar,2],sep="")
          write(best[Ipnt],ParFileName,append=T)
         }
        else
         {
          xx <- paste(ParName,"_",Ipar," ",Iqnt," NA ", ThePar$Initial[Ipar],sep="")
          write(ThePar$Initial[Ipar],ParFileName,append=T)
         }
        ## Now add o the Ouptput file
        write(xx,OutputFile,append=T)
       }
     } # Parameters within ParName
    } # ParNames
   write(paste("#Total estimated parameters:",Ipnt),OutputFile,append=T)

  IvarPnt <- Ipnt + 1

  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1

  #### =====================================================================================
  write("\n#Egg Production Total",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+GeneralSpecs$MaxProjYr
  if('MatBio'%in%rownames(stdrep)){
    tmp <- as.data.frame(stdrep[rownames(stdrep)=='MatBio',]) %>% mutate(year=(GeneralSpecs$Year1-GeneralSpecs$BurnIn) ) %>% filter(Estimate>0)
    tmp$year <- tmp$year+(0:(nrow(tmp)-1))
    tmp %<>% dplyr::select(year, est=Estimate,se='Std. Error')
    write("#Year est se",OutputFile,append=T)
    write(t(tmp),ncol=3,OutputFile,append=T)
    } else {
      tmp <- data.frame(year=(GeneralSpecs$Year1-GeneralSpecs$BurnIn), est=Report$MatBio, se=NA) %>% filter(est>0)
      tmp$year <- tmp$year+(0:(nrow(tmp)-1))
      write("#Year est nose",OutputFile,append=T)
      write(t(tmp),ncol=3,OutputFile,append=T)
    }

  write("\n#Egg Production by Area",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+GeneralSpecs$MaxProjYr
  if('MatBioArea'%in%rownames(stdrep)){
    tmp <- as.data.frame(stdrep[rownames(stdrep)=='MatBioArea',]) %>% mutate(year=(GeneralSpecs$Year1-GeneralSpecs$BurnIn),area=NA) %>% filter(Estimate>0)
    yr1 <- GeneralSpecs$Year1-GeneralSpecs$BurnIn
    yr2 <- yr1 + nrow(tmp)/GeneralSpecs$Narea - 1
    tmp$year <- rep(yr1:yr2, each=GeneralSpecs$Narea)
    tmp$area <- 1:GeneralSpecs$Narea
    tmp %<>% dplyr::select(year, area, est=Estimate,se='Std. Error')
    write("#Year area est se",OutputFile,append=T)
    write(t(tmp),ncol=4,OutputFile,append=T)
  } else {
    tmp <- as.data.frame(t(Report$MatBioArea))
    tmp <- tmp[tmp[,1]>0, ]
    tmp$year <- (GeneralSpecs$Year1-GeneralSpecs$BurnIn)
    tmp$year <- tmp$year+(0:(nrow(tmp)-1))
    colnames(tmp) <- c(1:GeneralSpecs$Narea, 'year')
    tmp %<>% pivot_longer(!year,names_to = 'area')
    write("#Year area est",OutputFile,append=T)
    write(t(tmp),ncol=3,OutputFile,append=T)
  }

  write("\n#Total Recruitment",OutputFile,append=T)
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  write("#Year est",OutputFile,append=T)
  write(t(cbind(Years,Report$Recruits)),ncol=2,OutputFile,append=T)

  write("\n#Recruitment by area (1+SD?)",OutputFile,append=T)
  Year <- GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)
  if('RecruitmentByArea'%in%rownames(stdrep)){
    tmp <- as.data.frame(stdrep[rownames(stdrep)=='RecruitmentByArea',]) %>% mutate(year=Year) %>% filter(Estimate>0)
    tmp$area <- 1:GeneralSpecs$Narea
    tmp$year <- rep(Year:(Year+nrow(tmp)/GeneralSpecs$Narea-1), each=GeneralSpecs$Narea)
    tmp %<>% dplyr::select(year, area, est=Estimate,se='Std. Error')
    write("#Year area est se",OutputFile,append=T)
    write(t(tmp),ncol=4,OutputFile,append=T)
  } else {
    tmp <- as.data.frame(t(Report$RecruitmentByArea))
    tmp <- tmp[tmp[,1]>0, ]
    tmp$year <- (GeneralSpecs$Year1-GeneralSpecs$BurnIn)
    tmp$year <- tmp$year+(0:(nrow(tmp)-1))
    colnames(tmp) <- c(1:GeneralSpecs$Narea, 'year')
    tmp %<>% pivot_longer(!year,names_to = 'area')
    write("#Year area est",OutputFile,append=T)
    write(t(tmp),ncol=3,OutputFile,append=T)
  }

  write("\n#Recruitment Fractions",OutputFile,append=T)
  lens <- Data$LowLenBin[1,1:dim(Report$RecruitFrac)[2]]
  write(t(cbind(lens,t(Report$RecruitFrac))),ncol=(dim(Report$RecruitFrac)[1]+1),OutputFile,append=T)

  write("\n#Legal Biomass by area (se?)",OutputFile,append=T) ## LegalBioAll runs first yr to last year
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  if('LegalBioAll'%in%rownames(stdrep)){
    tmp <- as.data.frame(stdrep[rownames(stdrep)=='LegalBioAll',]) %>% mutate(year=(GeneralSpecs$Year1-GeneralSpecs$BurnIn) ) %>% filter(Estimate>0)
    tmp$year <- tmp$year+(0:((nrow(tmp)/GeneralSpecs$Narea)-1))
    tmp$area <- rep(1:GeneralSpecs$Narea, each=nrow(tmp)/GeneralSpecs$Narea)
    tmp %<>% dplyr::select(year, area, est=Estimate,se='Std. Error')
    write("#Year area est se",OutputFile,append=T)
    write(t(tmp),ncol=4,OutputFile,append=T)
  } else {
    tmp <- as.data.frame((Report$LegalBioAll))
    tmp <- tmp[tmp[,1]>0, ]
    tmp$year <- (GeneralSpecs$Year1-GeneralSpecs$BurnIn)
    tmp$year <- tmp$year+(0:(nrow(tmp)-1))
    colnames(tmp) <- c(1:GeneralSpecs$Narea, 'year')
    tmp %<>% pivot_longer(!year,names_to = 'area')
    write("#Year area est",OutputFile,append=T)
    write(t(tmp),ncol=3,OutputFile,append=T)
  }

  write("\n#Legal Biomass sex by area (in predefined timestep weight of all lobster using 'Reference selectivity pattern')",OutputFile,append=T)  ## LegalBioAllbySex runs Burn-in to lastyear
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  tmp <- as.data.frame((Report$LegalBioAllbySex))
  tmp$year <- GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)
  tmp <- tmp[tmp[,1]>0, ]
  tmp$year <- tmp$year+(0:(nrow(tmp)-1))
  Sexes <- ifelse(GeneralSpecs$Nsex==1,'U',c('F','M'))
  colnames(tmp) <- c(paste(rep(Sexes,each=GeneralSpecs$Narea),rep(1:GeneralSpecs$Narea,GeneralSpecs$Nsex)), 'year')
  tmp %<>% pivot_longer(!year,names_to = 'sex area') %>% mutate(sex=substr(`sex area`,1,1), area=as.numeric(do.call('rbind',strsplit(`sex area`,' '))[,2]  )) %>% dplyr::select(year,sex,area,est=value)
  write("#Year sex area est",OutputFile,append=T)
  write(t(tmp),ncol=4,OutputFile,append=T)

  write("\n#Mature Biomass sex by area (in predefined timestep weight of all lobster of mature age)",OutputFile,append=T)
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1 # MatureBioAllbySex runs Burn-in to lastyear
  tmp <- as.data.frame((Report$MatureBioAllbySex))
  tmp$year <- GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)
  tmp <- tmp[tmp[,1]>0, ]
  tmp$year <- tmp$year+(0:(nrow(tmp)-1))
  Sexes <- ifelse(GeneralSpecs$Nsex==1,'U',c('F','M'))
  colnames(tmp) <- c(paste(rep(Sexes,each=GeneralSpecs$Narea),rep(1:GeneralSpecs$Narea,GeneralSpecs$Nsex)), 'year')
  tmp %<>% pivot_longer(!year,names_to = 'sex area') %>% mutate(sex=substr(`sex area`,1,1), area=as.numeric(do.call('rbind',strsplit(`sex area`,' '))[,2]  )) %>% dplyr::select(year,sex,area,est=value)
  write("#Year sex area est",OutputFile,append=T)
  write(t(tmp),ncol=4,OutputFile,append=T)


  write("\n#Biomass >76 by area (se?)",OutputFile,append=T)
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  if('sLegalBio76'%in%rownames(stdrep)){
    tmp <- as.data.frame(stdrep[rownames(stdrep)=='sLegalBio76',]) %>% mutate(year=GeneralSpecs$Year1) %>% filter(Estimate>0)
    tmp$year <- tmp$year+(0:((nrow(tmp)/GeneralSpecs$Narea)-1))
    tmp$area <- rep(1:GeneralSpecs$Narea, each=nrow(tmp)/GeneralSpecs$Narea)
    tmp %<>% dplyr::select(year, area, est=Estimate,se='Std. Error')
    write("#Year area est se",OutputFile,append=T)
    write(t(tmp),ncol=4,OutputFile,append=T)
  } else {
    tmp <- as.data.frame((Report$sLegalBio76))
    tmp$year <- (GeneralSpecs$Year1)
    tmp <- tmp[tmp[,1]>0,]
    tmp$year <- tmp$year+(0:(nrow(tmp)-1))
    colnames(tmp) <- c(1:GeneralSpecs$Narea, 'year')
    tmp %<>% pivot_longer(!year,names_to = 'area')
    write("#Year area est",OutputFile,append=T)
    write(t(tmp),ncol=3,OutputFile,append=T)
  }

  write("\n#Extended Legal Biomass by area (1+SD?)",OutputFile,append=T)
  for (Iarea in 1:GeneralSpecs$Narea)
  {
    Years <- (GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)):(GeneralSpecs$Year2+GeneralSpecs$MaxProjYr)
    Ncol <- 3
    LBOut <- cbind(rep(Iarea,length(Years)),Years,Report$LegalBioAll[,Iarea])
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  }

  write("\n#HarvestRate by zone",OutputFile,append=T)
  Years <- 1:(Nyears+1)+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1
  if('HarvestRate'%in%rownames(stdrep)){
    tmp <- as.data.frame(stdrep[rownames(stdrep)=='HarvestRate',]) %>% mutate(year=(GeneralSpecs$Year1)) %>% filter(Estimate>0)
    tmp$year <- tmp$year+(0:((nrow(tmp)/ControlSpecs$Nzone)-1))
    tmp$zone <- rep(LETTERS[1:ControlSpecs$Nzone], each=nrow(tmp)/ControlSpecs$Nzone)
    tmp %<>% dplyr::select(year, zone, est=Estimate,se='Std. Error')
    write("#Year zone est se",OutputFile,append=T)
    write(t(tmp),ncol=4,OutputFile,append=T)
  } else {
    tmp <- as.data.frame((Report$HarvestRate))
    tmp <- as.data.frame(tmp[tmp[,1]>0, ])
    tmp$year <- (GeneralSpecs$Year1)
    tmp$year <- tmp$year+(0:(nrow(tmp)-1))
    colnames(tmp) <- c(LETTERS[1:ControlSpecs$Nzone], 'year')
    tmp %<>% pivot_longer(!year,names_to = 'zone')
    write("#Year zone est",OutputFile,append=T)
    write(t(tmp),ncol=3,OutputFile,append=T)
  }

  write("\n#Natural Mortality by Area, Age and Year",OutputFile,append=T)
  for (Iage in 1:GeneralSpecs$Nage){
    for (Iarea in 1:GeneralSpecs$Narea){
    Years <- (GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)):(GeneralSpecs$Year2+GeneralSpecs$MaxProjYr)
    ## BurnIn+Nyear+MaxProjYr+1
    Ncol <- 4
    LBOut <- cbind(rep(Iarea,length(Years)),rep(Iage,length(Years)),Years,rep(Report$M[Iarea,Iage], length(Years)))
    write(t(LBOut),ncol=Ncol,OutputFile,append=T)
  } }

  write("\n#Virgin Biomass by area",OutputFile,append=T)
  write(Report$VirginBio,ncol=1,OutputFile,append=T)

  write("\n#Virgin Legal Biomass by area",OutputFile,append=T)
  write(Report$VirginLegalBio,ncol=1,OutputFile,append=T)

  write("\n#Harvest rate by zone (1+SD?)",OutputFile,append=T)
  for (Iarea in 1:ControlSpecs$Nzone)
   {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    HROut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$HarvestRate[,Iarea])
    for (Ivar in 1:ControlSpecs$NvarTypes)
     if (ControlSpecs$VarTypes[Ivar]==3)
      {
       HRVar <- rep(0,GeneralSpecs$Nyear)
       Ncol <- Ncol + 1
       for (IvarYr in 1:GeneralSpecs$Nyear)
        {
         IvarPnt <- IvarPnt + 1
         HRVar[IvarYr] <- fullrep[IvarPnt,2]
        }
       HROut <- cbind(HROut,HRVar)
      }
    write(t(HROut),ncol=Ncol,OutputFile,append=T)
   }

  write("\n#Harvest 76 rate by zone (1+SD?)",OutputFile,append=T)
  for (Iarea in 1:ControlSpecs$Nzone)
  {
    Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
    Ncol <- 3
    HROut <- cbind(rep(Iarea,GeneralSpecs$Nyear),Years,Report$HarvestRate76[,Iarea])
    write(t(HROut),ncol=Ncol,OutputFile,append=T)
  }


  #### =====================================================================================
  write("\n#Catches",OutputFile,append=T)
  write("#Fleet Area Year Step Observed Predicted Fishing_mortality",OutputFile,append=T)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
  {
    TotalCatch <- sum(TheData$Catch[,,Ifleet])
    if (TotalCatch > 0)
    {
      Iarea = ControlSpecs$Fleet_area[Ifleet]
      for (Iyear in 1:GeneralSpecs$Nyear)
        for (Istep in 1:GeneralSpecs$Nstep)
        {
          Summ <- c(Ifleet,Iarea+1,Iyear+GeneralSpecs$Year1-1,Istep,TheData$Catch[Iyear,Istep,Ifleet],Report$CatchCheck[Iyear,Istep,Ifleet],Report$Hrate[Iyear+max(GeneralSpecs$BurnIn),Istep,Ifleet])
          write(Summ,ncol=7,OutputFile,append=T)
        }
    }
  }

  #### =====================================================================================

  write("\n#Fishing Efficiency Commercial",OutputFile,append=T)
  write("#id Year Predicted",OutputFile,append=T)
  fcreep <- Report$CpueEcreep
  for (Iid in 1:ncol(fcreep)) {
    tmp <- cbind(rep(Iid, nrow(fcreep)), Data$Year1:(Data$Year1+nrow(fcreep)-1), fcreep[,Iid])
    write(t(tmp),ncol=3,OutputFile,append=T)
    #for(wte in 1:nrow(fcreep)) write(tmp[wte,],ncol=3,OutputFile,append=T)
      }

  #### =====================================================================================
  write("\n#Index data",OutputFile,append=T)
  write("#Data_set Fleet Sex Year Time_step Observed Relative_CV Predicted Residual",OutputFile,append=T)
  ResCpue <- matrix(0,nrow=TheData$Ncpue,ncol=9)
  for (Ipnt in 1:TheData$Ncpue)
   {
    ResCpue[Ipnt,] <- c(TheData$IndexI[Ipnt,]+c(1,1,0,GeneralSpecs$Year1,1),TheData$IndexR[Ipnt,],Report$PredCpue[Ipnt,])
   }
  write(t(ResCpue),OutputFile,ncol=9,append=T)
  write("\n#Data_set N Sigma Likelihood",OutputFile,append=T)
  for (IdataSet in (0:(Data$NcpueDataSeries-1)))
   {
    IndexPoints <- which(TheData$FixSigmaCpue==IdataSet)-1
    if (length(IndexPoints) >0)
     {
      Use <- TheData$IndexI[,1] %in% IndexPoints
      Summ <- paste(IdataSet+1,sum(Use),Report$SigmaCpue[IdataSet+1],Report$CpueLikeComp[IdataSet+1],Report$LambdaCpue2[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }
   }
  write("\n#Data_set N Q",OutputFile,append=T)
  for (IdataSet in (0:(Data$NcpueDataSeries-1)))
  {
    IndexPoints <- which(TheData$TreatQcpue==IdataSet)-1
    if (length(IndexPoints) >0)
    {
      Use <- TheData$IndexI[,1] %in% IndexPoints
      Summ <- paste(IdataSet+1,sum(Use),Report$CpueQ[IdataSet+1])
      write(Summ,OutputFile,append=T)
    }
  }

  #### --------------------------------------------------------------------------------------
  write("\n#Numbers data",OutputFile,append=T)
  write("#Index Fleet Year Time_step Observed Relative_CV Predicted Residual",OutputFile,append=T)
  ResNumbers <- matrix(0,nrow=TheData$Nnumbers,ncol=8)
  if (TheData$Nnumbers > 0)
   {
    for (Ipnt in 1:TheData$Nnumbers)
     {
      ResNumbers[Ipnt,] <- c(TheData$NumbersI[Ipnt,]+c(1,1,GeneralSpecs$Year1-1,1),TheData$NumbersR[Ipnt,],Report$PredNumbers[Ipnt,])
     }
    write(t(ResNumbers),OutputFile,ncol=8,append=T)
    write("\n#Data_set N Sigma Likelihood",OutputFile,append=T)
    for (IdataSet in (0:(Data$NcatchDataSeries-1)))
     {
      Use <- TheData$NumbersI[,1] == IdataSet
      Summ <- paste(IdataSet+1,sum(Use),Report$SigmaNumbers[IdataSet+1],Report$NumbersLikeComp[IdataSet+1],Report$LambdaNumbers2[IdataSet+1])
      write(Summ,OutputFile,append=T)
     }
  }

  write("\n#Growth Curves",OutputFile,append=T)   ## (Nyear,Narea,Nsex,Nage,MaxLen)
  write(c("#year area sex age lengthbins"),OutputFile,append=T)
  Years <- 1:(GeneralSpecs$Nyear)+GeneralSpecs$Year1-1
  tempout <- expand.grid(year=Years, area=1:GeneralSpecs$Narea, sex=1:GeneralSpecs$Nsex, age=1:GeneralSpecs$Nage)
  Mat <- matrix(NA, nrow=nrow(tempout), ncol=GeneralSpecs$Nlen[1])
  for (i in 1:nrow(tempout))   {
    GROut <- c(as.numeric(tempout[i,]),  Report$GrowthOut[which(Years==tempout$year[i]),tempout$area[i],tempout$sex[i],tempout$age[i],]  )
    write(GROut,ncol=length(GROut), OutputFile,append=T) }

  #### =====================================================================================
#   write("\n#Tagging data - Tag numbers by ",OutputFile,append=T)
#   write("#sex Relgrp area year step repSplit Obs Est",OutputFile,append=T)
#   for (Isex in 1:GeneralSpecs$Nsex)
#    for (Igrp in 1:TagSpecs$NtagGroups)
#     for (Iarea in 1:GeneralSpecs$Narea)
#      for (Iyear in 1:TagSpecs$NyearTags)
#       for (Istep in 1:GeneralSpecs$Nstep)
#        for (IrepSplit in 1:TagSpecs$NrepSplit)
#         if (TagSpecs$RecapObs[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep] > 0)
#          {
#           VAL1<- TagSpecs$RecapObs[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep]*TagSpecs$NrelTotal[Isex,Igrp]
#           VAL2 <- Report$RecapNum[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep]*TagSpecs$NrelTotal[Isex,Igrp]
#           Summ <- paste(Isex,Igrp,Iarea,Iyear+TagSpecs$Year1Tag[1]-1,Istep,IrepSplit,VAL1,VAL2)
#           write(Summ,OutputFile,append=T)
#          }
#
#   write("\n#Tagging length data - Tag numbers by ",OutputFile,append=T)
#   write("#sex Relgrp area repSplit size ObsTot ObsProp EstProp",OutputFile,append=T)
#   for (Isex in 1:GeneralSpecs$Nsex)
#    for (Igrp in 1:TagSpecs$NtagGroups)
#     for (Iarea in 1:GeneralSpecs$Narea)
#      for (IrepSplit in 1:TagSpecs$NrepSplit)
#       if (TagSpecs$FitTagSizes[IrepSplit] == 1)
#        {
#         ObsSS = 0
#         for (Iyear in 1:TagSpecs$NyearTags)
#          for (Istep in 1:GeneralSpecs$Nstep)
#           ObsSS = ObsSS + TagSpecs$TagRec[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep,1]
#         if (ObsSS > 0)
#          {
#           Vec <- rep(0,GeneralSpecs$Nlen[Isex])
#           for (Isize in 1:GeneralSpecs$Nlen[Isex])
#            {
#             ObsEE <- 0
#             for (Iyear in 1:TagSpecs$NyearTags)
# 	     for (Istep in 1:GeneralSpecs$Nstep)
# 	      ObsEE = ObsEE + TagSpecs$TagRec[Isex,Igrp,Iarea,IrepSplit,Iyear,Istep,1+Isize]
# 	     Vec[Isize] <- ObsEE
# 	     PredEE <- Report$PredTagSize[Isex,Igrp,Iarea,Isize]
#              Summ <- paste(Isex,Igrp,Iarea,IrepSplit,Isize,ObsSS,Vec[Isize],PredEE)
#              write(Summ,OutputFile,append=T)
#             }
#
#          }
#        }



  #### =====================================================================================
  write("\n#Length data tuning",OutputFile,append=T)
  #TagSpecs

  # Calculate Francis weights and McAllister-Ianelli mode
  write("#Fleet Sex Npnts Francis_Multiplier",OutputFile,append=T)
  EffN <- rep(0,TheData$NlenComp)
  for (Ifleet in 1:GeneralSpecs$Nfleet)
   for (Isex in 1:GeneralSpecs$Nsex)
     {
      Use <-  TheData$LenCompI[,1]+1 == Ifleet & TheData$LenCompI[,2]+1 == Isex
      Indexes <- c(1:TheData$NlenComp)[Use]
      Residuals <- NULL
      for (Index2 in 1:length(Indexes))
       {
        Index <- Indexes[Index2]
        Top <- 0; Bot <- 0
        Nlens <- GeneralSpecs$Nlen[Isex]
        Length <- GeneralSpecs$MidLenBin[Isex,1:Nlens]
        ObsProp <- TheData$LenCompR[Index,1:Nlens]/sum(TheData$LenCompR[Index,1:Nlens])
        PredProp <- Report$PredLengthComp[Index,1:Nlens]
        for (Ilen in 1:Nlens)
         {
          Bot <- Bot + (ObsProp[Ilen] - PredProp[Ilen])^2
          Top <- Top + PredProp[Ilen]*(1.0-PredProp[Ilen])
         }
        MeanOL <- sum(ObsProp*Length)
        MeanPL <- sum(PredProp*Length)
        SD <- sqrt((sum(PredProp*Length^2)-MeanPL^2)/TheData$Stage1W[Index])
        Residual <- (MeanOL-MeanPL)/SD
        Residuals <- c(Residuals,Residual)
        EffN[Index] <- Top/Bot
       }
      if (!is.na(Residuals[1]))
       {
        LenWghtMultipliers <- 1.0/var(Residuals)
        if(is.na(LenWghtMultipliers)) LenWghtMultipliers <- 1
        Summ <- c(Ifleet,Isex,length(Residuals),LenWghtMultipliers)
        write(Summ,OutputFile,append=T,ncol=5)
       }
     }
  write("\n#Obs/Pred Fleet Sex Year Step Nsamp EffN proportions",OutputFile,append=T)
  for (Ipnt in 1:TheData$NlenComp)
   {
    Ifleet <-  TheData$LenCompI[Ipnt,1]+1
    Isex <- TheData$LenCompI[Ipnt,2]+1
    Iyear <-  TheData$LenCompI[Ipnt,3]+GeneralSpecs$Year1
    Istep <- TheData$LenCompI[Ipnt,4]+1
    Obs <- c(round(as.vector(TheData$LenCompR[Ipnt,1:GeneralSpecs$Nlen[Isex]]),5))
    Summ <- paste(c("O",Ifleet,Isex,Iyear,Istep,TheData$Stage1W[Ipnt],round(EffN[Ipnt],3),paste(Obs)))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$Nlen[Isex])
    Pred <- c(round(as.vector(Report$PredLengthComp[Ipnt,1:GeneralSpecs$Nlen[Isex]]),5))
    Summ <- paste(c("P",Ifleet,Isex,Iyear,Istep,TheData$Stage1W[Ipnt],round(EffN[Ipnt],3),paste(Pred)))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$Nlen[Isex])
   }

  #### --------------------------------------------------------------------------------------
  write("\n#Larval data",OutputFile,append=T)
  write("#Area Year Observed SD Predicted Residual",OutputFile,append=T)
  ResLarval <- matrix(0,nrow=TheData$NLarvalData,ncol=6)
  if (TheData$NLarvalData > 0)
   {
    for (Ipnt in 1:TheData$NLarvalData)
     {
      ResLarval[Ipnt,] <- c(TheData$Lar_dataI[Ipnt,]+c(1,GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)),TheData$Lar_dataR[Ipnt,],Report$PredLarval[Ipnt,])
     }
   # print(ResLarval)
    write(t(ResLarval),OutputFile,ncol=6,append=T)
   }


  #### =====================================================================================
  write("\n#Full Selectivity",OutputFile,append=T)
  NselexPatterns <-length(Report$ActSelex[,1])
  for (Ipattern in 1:NselexPatterns)
   {
    Summ <- c(Ipattern,c(as.vector(round(Report$ActSelex[Ipattern,],5))))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$MaxLen)
  }

  write("\n#Legal Selectivity by sex age fleet year tstep",OutputFile,append=T)
  out <- expand.grid(sex=1:GeneralSpecs$Nsex, age=1:GeneralSpecs$Nage, fleet=1:GeneralSpecs$Nfleet, year=1:GeneralSpecs$Nyear, tstep=1:GeneralSpecs$Nstep)
  #dim(Report$ActLegal)
  #dim(Data$LegalFleetPnt)
  out2 <- matrix(NA, nrow=nrow(out),ncol=ncol(Report$ActLegal))
  for (I in 1:nrow(out))  {
    out2 <-  Report$ActLegal[ (Data$LegalFleetPnt[out$sex[I],out$age[I],out$fleet[I],out$year[I],out$tstep[I]]+1), ]
    out3 <- c(out$sex[I],out$age[I],out$fleet[I],(GeneralSpecs$Year1+out$year[I]-1),out$tstep[I],out2)
    write(out3,OutputFile,append=T,ncol=length(out3))
  }

  write("\n#Retention",OutputFile,append=T)
  NretenPatterns <-length(Report$ActReten[,1])
  for (Ipattern in 1:NretenPatterns)
   {
    Summ <- c(Ipattern,c(as.vector(round(Report$ActReten[Ipattern,],5))))
    write(Summ,OutputFile,append=T,ncol=7+GeneralSpecs$MaxLen)
   }

  #### =====================================================================================
  write("\n#Recruitment",OutputFile,append=T)

  write("#Recruitment patterns by timestep, area and sex ",OutputFile,append=T)
  write("#year, tstep, sex area1, area2, area3, area4, area5, area6, area7, area8",OutputFile,append=T)
  for (Iyear in 1:dim(Report$ActRecruitAreaSexDist)[1]){
    for (Istep in 1:dim(Report$ActRecruitAreaSexDist)[2]){
        for (Isex in 1:dim(Report$ActRecruitAreaSexDist)[4]){
     write(paste(Iyear,Istep,Isex,Report$ActRecruitAreaSexDist[Iyear,Istep, ,Isex]),OutputFile,append=T)
    }}}

  write("#Recruitment patterns by sex and size ",OutputFile,append=T)
  write("#Pattern, sex, size ",OutputFile,append=T)
  for (Ipattern in 1:dim(Report$ActRecruitLenDist)[1]){
   for (Isex in 1:dim(Report$ActRecruitLenDist)[2]){
     write(paste(Ipattern,Iarea,Isex,Report$ActRecruitLenDist[Ipattern,Isex,]),OutputFile,append=T)
    }}

  #### =====================================================================================
  write("\n#Movement",OutputFile,append=T)

  write("#Movement patterns by size ",OutputFile,append=T)
  write("#Pattern, size",OutputFile,append=T)
  NmovePatterns <- length(Report$ActMove[,1])
  for (Ipattern in 1:NmovePatterns)
   write(paste(Ipattern,Report$ActMove[Ipattern,]),OutputFile,append=T)

  #### =====================================================================================
  write("\n#Initial N-matrix\n#Area\tSex\tAge",OutputFile,append=T)
  Nout <-rep(0,3+GeneralSpecs$MaxLen)

  for (Iarea in 1:GeneralSpecs$Narea)
   for (Isex in 1:GeneralSpecs$Nsex)
    for (Iage in 1:(GeneralSpecs$Nage))
     {
      Nout[1:3]  <- c(Iarea,Isex,Iage)
      Nout[4:(3+GeneralSpecs$MaxLen)] <- round(Report$Ninit[Iarea,Isex,Iage,],2)
      write(t(Nout),OutputFile,append=T,ncol=3+GeneralSpecs$MaxLen)
     }

  #### =====================================================================================
  write("\n#Simplified N-matrix",OutputFile,append=T)
  ncol <- GeneralSpecs$MaxLen+4
  for (Iarea in 1:GeneralSpecs$Narea)
   for (Isex in 1:GeneralSpecs$Nsex)
    {
     write("#Area Sex Year(raw) Year(act) Lengths",OutputFile,append=T,ncol=100)
     Nrow = (Nyears)
     Nout <- matrix(0,Nrow,4+GeneralSpecs$MaxLen)
     Ipnt <- 0
     for (Iyear in 1:Nyears)
      {
       Ipnt <- Ipnt + 1
       Nout[Ipnt,1:4]  <- c(Iarea,Isex,Iyear-max(GeneralSpecs$BurnIn)+1,Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1)
       for (Iage in 1:(GeneralSpecs$Nage))
        Nout[Ipnt,(5:(4+GeneralSpecs$MaxLen))] <- Nout[Ipnt,(5:(4+GeneralSpecs$MaxLen))]+ round(Report$N[Iarea,Iyear,1,Isex,Iage,],2)
      }
     write(t(Nout),OutputFile,append=T,ncol=ncol)
   }


  #### =====================================================================================
  write("\n#F-matrix",OutputFile,append=T)
  write("#Year Time-step F_by_fleet",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)
  Nrow = GeneralSpecs$Nstep*Nyears
  Fout <- matrix(0,Nrow,2+GeneralSpecs$Nfleet)
  Ipnt <- 0
  for (Iyear in 1:Nyears)
    for (Istep in 1:GeneralSpecs$Nstep)
    {
     Ipnt <- Ipnt + 1
     Fout[Ipnt,1:2]  <- c(Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1,Istep-1)
     Fout[Ipnt,(3:(2+GeneralSpecs$Nfleet))] <- round(Report$Hrate[Iyear,Istep,],5)
    }
  write(t(Fout),OutputFile,append=T,ncol=GeneralSpecs$Nfleet+2)

  write("\n#Full N-matrix",OutputFile,append=T)
  write("#Area Sex Age Year Time-step Lengths",OutputFile,append=T)
  Nyears <- GeneralSpecs$Nyear+max(GeneralSpecs$BurnIn)+1
  Nrow = (GeneralSpecs$Narea*GeneralSpecs$Nsex*GeneralSpecs$Nage*Nyears*GeneralSpecs$Nstep)
  Nout <- matrix(NA,Nrow,5+GeneralSpecs$MaxLen)
  ncol <- GeneralSpecs$MaxLen+5
  Ipnt <- 0
  for (Iarea in 1:GeneralSpecs$Narea)
    for (Isex in 1:GeneralSpecs$Nsex)
      for (Iage in 1:GeneralSpecs$Nage)
       for (Iyear in 1:Nyears)
        for (Istep in 1:GeneralSpecs$Nstep)
         {
          Ipnt <- Ipnt + 1
          Nout[Ipnt,1:5]  <- c(Iarea,Isex,Iage,Iyear+GeneralSpecs$Year1-max(GeneralSpecs$BurnIn)-1,Istep)
          Nout[Ipnt,(6:(5+GeneralSpecs$MaxLen))] <- 0
          Nout[Ipnt,(6:(5+GeneralSpecs$MaxLen))] <- Report$N[Iarea,Iyear,Istep,Isex,Iage,]
         }
  write(t(round(Nout,1)),OutputFile,append=T,ncol=ncol)
  setwd(dirname(getwd()))

}

