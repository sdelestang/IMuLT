#' Configure Parameter Initialization and Estimation Maps for Current Phase
#'
#' Internal function that prepares parameter structures for RTMB/TMB optimization
#' by setting initial values, bounds, and creating a "map" that controls which
#' parameters are estimated versus fixed for the current estimation phase.
#'
#' @param ParOld Numeric vector of parameter values from previous phase(s). Used
#'   to initialize parameters that were estimated in earlier phases at their
#'   converged values
#' @param parameters Named list of parameter arrays/vectors matching the RTMB/TMB
#'   model structure. Must include all parameter types defined in InitialVars
#' @param InitialVars Nested list from ReadInitialValues() containing Initial,
#'   Bnd, and Phase components for each parameter type
#' @param CurrPhase Integer indicating the current estimation phase (1, 2, 3, etc.)
#'
#' @return List containing components needed for RTMB/TMB optimization:
#' \itemize{
#'   \item map - Named list of factor vectors controlling parameter estimation.
#'     NA values indicate fixed parameters; integer values group parameters
#'     (same integer = estimated as single shared parameter)
#'   \item parameters - Parameter list with values initialized appropriately
#'     for the current phase
#'   \item EstVec - Numeric vector of parameter values to be estimated in
#'     current phase (initial values for optimization)
#'   \item lowBnd - Lower bounds for parameters being estimated
#'   \item uppBnd - Upper bounds for parameters being estimated
#' }
#'
#' @details
#' This function implements phased parameter estimation for RTMB/TMB models.
#' Parameters are estimated in sequential phases, with earlier phases converging
#' before later parameters are activated. This approach:
#' \itemize{
#'   \item Improves optimization stability for complex models
#'   \item Reduces parameter correlations
#'   \item Allows core parameters to stabilize first
#'   \item Makes troubleshooting easier
#' }
#'
#' **Phase logic for each parameter**:
#' \itemize{
#'   \item Phase < 0 or 0: Always fixed at initial value
#'   \item Phase < CurrPhase: Previously estimated, fixed at converged value (from ParOld)
#'   \item Phase = CurrPhase: Currently being estimated, included in map and EstVec
#'   \item Phase > CurrPhase: Not yet estimated, fixed at initial value
#' }
#'
#' The "map" is RTMB/TMB's mechanism for controlling which parameters are estimated.
#' It's a named list parallel to the parameter structure, containing factor vectors
#' where NA indicates a fixed parameter and integers identify parameters to estimate.
#'
#' **Special handling of dummy parameter**: If no parameters are estimated in the
#' current phase, a dummy parameter is activated to prevent optimization failures.
#' This ensures the model can run even when all biological parameters are fixed.
#'
#' The function uses a nested SinglePhase() helper that processes each parameter
#' type individually, building up the complete map, estimation vector, and bounds.
#'
#' @examples
#' \dontrun{
#' # Phase 1: Estimate only core parameters
#' phase1 <- SetInitialAndPhases(ParOld = NULL, parameters, InitialVars, CurrPhase = 1)
#' obj <- RTMB::MakeADFun(data, phase1$parameters, map = phase1$map)
#' opt1 <- nlminb(phase1$EstVec, obj$fn, obj$gr,
#'                lower = phase1$lowBnd, upper = phase1$uppBnd)
#'
#' # Phase 2: Add selectivity parameters, holding phase 1 at converged values
#' phase2 <- SetInitialAndPhases(opt1$par, parameters, InitialVars, CurrPhase = 2)
#' obj <- RTMB::MakeADFun(data, phase2$parameters, map = phase2$map)
#' opt2 <- nlminb(phase2$EstVec, obj$fn, obj$gr,
#'                lower = phase2$lowBnd, upper = phase2$uppBnd)
#' }
#'
#' @keywords internal
SetInitialAndPhases <- function(ParOld,parameters,InitialVars,CurrPhase)
{

  # Pointer to old parameters
  Ipnt <- 0
  SinglePhase<-function(Est,Bnd,Phase,CurrPhase)
  {
    map <- NULL
    estVec <- NULL; lowBnd <- NULL; uppBnd <- NULL
    Npar <- length(Est)

    for (Ipar in 1:Npar)
    {
      if (Phase[Ipar] > 0 & Phase[Ipar] <= (CurrPhase-1))
      {
        Ipnt <<- Ipnt + 1
        Est[Ipar] <- ParOld[Ipnt]
      }
      if (Phase[Ipar] > 0 & Phase[Ipar] <= CurrPhase)
      {
        map <- c(map,Ipar)
        estVec <- c(estVec,Est[Ipar])
        lowBnd <- c(lowBnd,Bnd[Ipar,1])
        uppBnd <- c(uppBnd,Bnd[Ipar,2])
      }
      else
        map <- c(map,as.factor(NA))
    }
    ReturnObj <- NULL
    ReturnObj$map <- as.factor(map)
    ReturnObj$estVec <- estVec
    ReturnObj$lowBnd <- lowBnd
    ReturnObj$uppBnd <- uppBnd
    return(ReturnObj)
  }  # SinglePhase

  map <- list()
  estvec <- NULL; lowBnd <- NULL; uppBnd <- NULL
  for (ParName in names(parameters))
  {
    if (ParName != "dummy")
    {
      if (length(InitialVars[[ParName]]$Initial) >0)
      {
        ThePar <- InitialVars[[ParName]]
        PhaseOut <- SinglePhase(ThePar$Initial,ThePar$Bnd,ThePar$Phase,CurrPhase)
        map <- append(map,list(ParName=PhaseOut$map))
        estvec <- c(estvec,PhaseOut$estVec)
        lowBnd <- c(lowBnd,PhaseOut$lowBnd)
        uppBnd <- c(uppBnd,PhaseOut$uppBnd)
      }
      else
      {
        estvec <- c(estvec,0)
        map <- append(map,list(ParName=factor(NA)))
      }
    }
    else
    {
      if(max(as.vector(sapply(map, function(x) max(as.numeric(!is.na(x))))))==0) {
        map <- append(map,list(dummy=factor(1)))
        estvec <- c(estvec,0)
        lowBnd <- c(lowBnd,-1)
        uppBnd <- c(uppBnd,1)
      } else { map <- append(map,list(dummy=factor(NA)))}
    }
  }
  names(map) <- names(parameters)
  # dummy
  ReturnObj <- NULL
  ReturnObj$map = map
  ReturnObj$parameters <- parameters
  ReturnObj$EstVec <- estvec
  ReturnObj$lowBnd <- lowBnd
  ReturnObj$uppBnd <- uppBnd
  return(ReturnObj)
}

#' Get Maximum Estimation Phase Number
#'
#' Helper function to determine the highest estimation phase across all parameter
#' types in the model specification.
#'
#' @param InitialVars Nested list from ReadInitialValues() containing parameter
#'   specifications with Phase components
#'
#' @return Integer indicating the maximum phase number found across all parameters
#'
#' @details
#' This function scans through all parameter types and finds the highest phase
#' number specified. This is useful for:
#' \itemize{
#'   \item Determining how many optimization phases are needed
#'   \item Loop control in sequential phase estimation
#'   \item Validation that all phases are sequential
#' }
#'
#' The function returns a minimum value of 1 even if all parameters have negative
#' phases (all fixed), ensuring at least one phase exists for model evaluation.
#'
#' @examples
#' \dontrun{
#' InitVals <- ReadInitialValues(...)
#' max_phase <- getPhase(InitVals)
#' # Run optimization through all phases
#' for (phase in 1:max_phase) {
#'   # ... estimation code ...
#' }
#' }
#'
#' @seealso \code{\link{SetInitialAndPhases}} for phase-based parameter configuration
#'
#' @keywords internal
getPhase <- function(InitialVars){
  mX <- 1
  num <- length(InitialVars)
  for(n in 1:num){
    tmp <- InitialVars[[n]]
    mX <- max(c(mX, tmp$Phase))
  }
  return(mX)
}

csemod <- function(x){
  mod <-  as.numeric(dlg_input(c('Choose a model:',paste(1:length(x), x, sep=(" : ") )), 1)$res)
  if (!length(mod)) {# The user clicked the 'cancel' button
    cat(paste("OK, the default model is",x[1],"\n"))
  } else {
    cat(paste("Model", x[mod], "has been chosen"), "\n")
  }
  return(x[mod])}

#' Interactively Select IMuLT Model Run Directory
#'
#' Opens an interactive dialog to select from available model run directories
#' and sets the working directory to the chosen run. If no selection is made,
#' the first matching directory is used by default.
#'
#' @param pattern Character string pattern for matching run directories.
#'   Default is 'Run' to match directories like "8Area8AgeRun92_23".
#'
#' @return Invisibly returns the name of the selected directory as a character string.
#'   The function also changes the working directory to the selected run folder as a side effect.
#'
#' @details
#' This function is typically used at the start of an analysis workflow to select
#' which model run to work with. It:
#' \itemize{
#'   \item Searches for directories matching the pattern in the current working directory
#'   \item Displays an interactive dialog listing all matching directories
#'   \item Changes the working directory to the selected folder
#'   \item Returns the directory name for potential further use
#' }
#'
#' @examples
#' \dontrun{
#' # Select from available "Run" directories interactively
#' chosen_run <- choose_model()
#'
#' # Select from directories matching a custom pattern
#' chosen_run <- choose_model(pattern = "AgeRun")
#'
#' # The working directory is now set to the chosen run
#' getwd()
#' }
#'
#' @seealso \code{\link{BuildInputFiles}} for creating new model run directories
#'
#' @export
choose_model <- function(pattern = 'Run') {
  x <- list.files(pattern = pattern)

  mod <- as.numeric(dlg_input(c('Choose a model:',
                                paste(1:length(x), x, sep = " : ")),
                              1)$res)

  if (!length(mod)) {
    cat(paste("OK, the default model is", x[1], "\n"))
    selected <- x[1]
  } else {
    cat(paste("Model", x[mod], "has been chosen"), "\n")
    selected <- x[mod]
  }

  setwd(file.path(getwd(), selected))
  invisible(selected)  # Returns but doesn't auto-print
}

#' Update Model Input Files with Estimated Parameters
#'
#' Updates IMuLT model input files (.DAT files) with parameter estimates from
#' a completed model run. This allows you to use converged parameter values
#' as starting values for subsequent runs or modified scenarios.
#'
#' @param todo Character string indicating whether to update parameters.
#'   Options are 'Yes' or 'No'. If not specified (default ' '), an interactive
#'   dialog will prompt the user to choose.
#'
#' @return NULL (invisibly). The function modifies input files as a side effect.
#'
#' @details
#' The function reads final parameter estimates from 'Output/model final.par' and
#' updates the following input files with these values:
#' \itemize{
#'   \item CONTROL.DAT - Main parameters, Q parameters, efficiency parameters, recruitment deviations, and spatial recruitment deviations
#'   \item RECRUITSPEC.DAT - Recruitment parameters and puerulus power parameters
#'   \item MOVESPEC.DAT - Movement/migration parameters
#'   \item SELEXSPEC.DAT - Selectivity parameters
#' }
#'
#' This is useful when:
#' \itemize{
#'   \item Starting a new run from converged values
#'   \item Running projections with estimated parameters
#'   \item Testing model sensitivity with slightly modified starting values
#' }
#'
#' @note This function must be run from within a model run directory that contains
#'   the input .DAT files and an 'Output/model final.par' file from a completed run.
#'
#' @examples
#' \dontrun{
#' # Set working directory to a completed model run
#' choose_model()
#'
#' # Update parameters interactively (prompts for Yes/No)
#' UpdatePars()
#'
#' # Update parameters without prompting
#' UpdatePars("Yes")
#'
#' # Skip updating (useful in scripts)
#' UpdatePars("No")
#' }
#'
#' @seealso \code{\link{choose_model}} for selecting a model run directory
#'
#' @export
UpdatePars <- function(todo=' '){
  if(todo==' ') {todo <- dlg_list(c('Yes','No','                '),  title=c('Update Parameters?                    '))$res  }
  ## Get estimated parameters  KeyWord <- locs$id[i]
  if(todo=='Yes'){
find <- function(KeyWord, DataFile, Offset){
    KeyWord <- unlist(strsplit(as.character(KeyWord),' '))
    if(length(KeyWord)==1) pos1 <- which(grepl(KeyWord,DataFile[,2]))+Offset
    if(length(KeyWord)==2) pos1 <- which(2==(grepl(KeyWord[1],DataFile[,2])+grepl(KeyWord[2],DataFile[,3])))+Offset
    if(length(KeyWord)==3) pos1 <- which(3==(grepl(KeyWord[1],DataFile[,2])+grepl(KeyWord[2],DataFile[,3])+grepl(KeyWord[3],DataFile[,4])))+Offset
    return(pos1)}

  pout <- read.delim('Output/model final.par',sep='\t',stringsAsFactors =F)
  unique(pout$name)
  locs <-data.frame(par=c('MainPars', 'RecruitPars','PuerPowPars','RecDevs','Qpars','efpars','RecSpatDevs','MovePars','SelPars'), file=c('CONTROL.DAT', 'RECRUITSPEC.DAT', 'RECRUITSPEC.DAT','CONTROL.DAT','CONTROL.DAT','CONTROL.DAT','CONTROL.DAT','MOVESPEC.DAT','SELEXSPEC.DAT'), id=c('Basic parameters','Recuitment1 parameters','Puerulus Power for','Prespecify_rec_devs','Q parameters','Efficiency parameters','Prespecify_spatial_rec_devs','Movement parameters','Selectivity Parameters'), off=c(1,2,3,2,1,2,2,2,2), col=c(3,3,3,1,3,3,1,3,3))

  for (i in 1:nrow(locs)){
    DataFile <- read.table(locs$file[i],comment.char = "?",fill=T, blank.lines.skip=F,stringsAsFactors=F,col.names=1:200)
    ptmp <- pout$est[grepl(locs$par[i], pout$name)]
    roff <- locs$off[i]
    coff <- locs$col[i]
    pos <- find(c(unlist(strsplit(locs$id[i],' '))), DataFile, roff)
    DataFile[pos:(pos+length(ptmp)-1),coff] <- ptmp
    write.table(DataFile, locs$file[i], na=" ", sep=" ", row.names = F, col.names = F, quote=F)
    print(paste("Parameters upated: ", locs$par[i]))
  }
}}


#' Fit IMuLT Stock Assessment Model Using TMB
#'
#' Main function for estimating model parameters using Template Model Builder (TMB)
#' and sequential phased optimization. Minimizes negative log-likelihood using
#' nlminb optimizer with bounds constraints and optional report generation.
#'
#' @param phit Integer. Maximum number of function evaluations for intermediate
#'   estimation phases. Default is 500.
#' @param lphit Integer. Maximum number of function evaluations for the final
#'   estimation phase. Default is 1000 (allows more iterations for convergence).
#' @param mxph Integer. Maximum phase number to run. Default is MaxPhase (set by
#'   LoadPars()). Use lower values to run partial estimation sequences.
#' @param PrintLag Integer. Progress is printed every PrintLag function calls.
#'   Default is 50. Lower values give more frequent updates.
#' @param report Logical. If TRUE, generates full diagnostic report including
#'   SD report and saves BigSave.lda and Output.RL files. Default is FALSE
#'   (faster, for intermediate runs). Set TRUE for final model run.
#'
#' @return NULL. Creates output files as side effects:
#' \itemize{
#'   \item Output/model [phase].par - Parameter values after each phase
#'   \item Output/model final.par - Final converged parameters
#'   \item Output/BigSave.lda - Complete model object (if report=TRUE)
#'   \item Output/Output.RL - Formatted results for diagnostics (if report=TRUE)
#' }
#'
#' @details
#' The function implements sequential phased estimation:
#' \enumerate{
#'   \item For each phase (1 to MaxPhase):
#'   \item Sets active parameters based on phase specification
#'   \item Initializes TMB model object
#'   \item Runs nlminb optimization with parameter bounds
#'   \item Saves parameter estimates
#'   \item Uses estimates as starting values for next phase
#' }
#'
#' Progress monitoring displays:
#' \itemize{
#'   \item Current phase and iteration number
#'   \item Negative log-likelihood value
#'   \item Percent improvement from previous best
#'   \item Number of active parameters
#' }
#'
#' Convergence is indicated by:
#' \itemize{
#'   \item Convergence code = 0 (successful)
#'   \item Maximum gradient < 0.001 (well converged)
#'   \item Maximum gradient < 0.01 (acceptable)
#' }
#'
#' @note
#' \itemize{
#'   \item Must run LoadPars() before calling this function
#'   \item Set report=TRUE only for final production runs (much slower)
#'   \item Monitor convergence - may need to increase lphit if not converging
#'   \item Large models may take hours to run with report=TRUE
#' }
#'
#' @examples
#' \dontrun{
#' # Standard workflow
#' choose_model()
#' LoadPars()
#'
#' # Quick test run (no report)
#' SolveModelNew(phit = 100, lphit = 200, report = FALSE)
#'
#' # Full production run with reports
#' SolveModelNew(phit = 500, lphit = 1000, report = TRUE)
#'
#' # Run only first 2 phases for testing
#' SolveModelNew(mxph = 2, report = FALSE)
#'
#' # More frequent progress updates
#' SolveModelNew(PrintLag = 10)
#' }
#'
#' @seealso
#' \code{\link{LoadPars}} for loading parameters before estimation,
#' \code{\link{AdjustPhase}} for modifying estimation phases,
#' \code{\link{MakeDiagReport}} for generating diagnostic outputs,
#' \code{\link{choose_model}} for selecting model directory
#'
#' @export
SolveModelNew <- function(phit=500,lphit=1000, mxph=MaxPhase, PrintLag = 50, report=F){
  MaxPhase=ifelse(mxph==0,1,mxph)
  for (CurrPhase in 1:MaxPhase) {
    MaXeVaL <- ifelse(CurrPhase<MaxPhase, phit, lphit)
    parameters <- list(MainPars=InitialVars$MainPars$Initial,RecruitPars=InitialVars$RecruitPars$Initial,PuerPowPars=InitialVars$PuerPowPars$Initial,SelPars=InitialVars$SelPars$Initial,RetPars=InitialVars$RetPars$Initial,RecDevs=InitialVars$RecDevs$Initial,Qpars=InitialVars$Qpars$Initial,efpars=InitialVars$efpars$Initial,InitPars=InitialVars$InitPars$Initial,RecSpatDevs=InitialVars$RecSpatDevs$Initial,MovePars=InitialVars$MovePars$Initial,GrowthPars=InitialVars$GrowthPars$Initial,dummy=0)
    RunSpecs <- SetInitialAndPhases(ParOld,parameters,InitialVars,CurrPhase=CurrPhase)  # Set parameters and mapping
    ## Make model
    cat("Making model object that will solve for",sum(!is.na(unlist(RunSpecs$map))) ,"parameters.","Phase =",CurrPhase,"\n")
    pnames <- names(unlist(RunSpecs$map)[!is.na(unlist(RunSpecs$map))]);
    nam <- stringr::str_extract(pnames, "[\\p{Letter}]+")
    num <- stringr::str_extract(pnames, "\\d+$")
    unnam <- nam[!duplicated(nam)]
    for(iii in 1:length(unnam))  { print(paste(unnam[iii], length(num[nam==unnam[iii]]),'parameters'))  }
    model <- MakeADFun(Data, parameters, map=RunSpecs$map, DLL="IMuLT",silent=T)
    BestFn <- model$fn()
    initBestFn <- BestFn
    FnCallNo <<- 0;
    model$fn_Orig <- model$fn
    yy <- 1e+10
    model$fn <- function(x)  {
      tyy <- model$fn_Orig(x)
      yy <<- ifelse(is.na(tyy),yy,tyy)
      FnCallNo <<- FnCallNo + 1
      if(BestFn>yy){
        if ((FnCallNo %% PrintLag)==0) {
          cat("Phase ", CurrPhase," ",FnCallNo," -LogLike / Delta: ",yy," / ",round(100*(1-(yy/BestFn)),6),"%; npar = ",length(x),"\n",sep="")
          BestFn <<- yy }
      }
      return(tyy);  }
    model$par <- RunSpecs$EstVec  # Assign new parameters associated with the correct phase
    # Run model
    BestFn <- model$fn(model$par)
    initBestFn <- BestFn
    mout<-nlminb(model$par,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL, eval.max=MaXeVaL))
    pars <- mout$par; names(pars) <- pnames; ParOld <- mout$par;
    Grad <- model$gr(mout$par)
    grad <- round(max(abs(Grad)),6)
    cat("Likelihood: ",round(initBestFn,6),' to ' ,round(mout$objective,6),"Convergence:",ifelse(mout$convergence==0,'Yes','No')," (",mout$convergence,") ","Max Gradient:",round(max(abs(model$gr(mout$par))),6)," Interations:",mout$iterations,"Evalutions:",mout$evaluations,"\n")
    # Store and save parameters
    pout <- unlist(parameters); pout[names(pout)%in%names(pars)] <- pars; suffix <- ifelse(CurrPhase==MaxPhase," final", CurrPhase); write.table(pout, paste("Output/model",suffix,".par",sep=""), sep='\t', col.names = c('name\test'), quote=F)
  if(report==T & CurrPhase==MaxPhase){       cat("Making report object.\n")
      print("Loading report")
      Report <- model$report()
      best <- mout$par
      print("Loading SD report (can take quite a long time)")
      SDrep <- sdreport(model)
      fullrep <- summary(SDrep)
      BigSave <-NULL
      BigSave$Report <- Report
      BigSave$SDrep <- SDrep
      BigSave$map <- RunSpecs$map
      BigSave$Data <- Data
      BigSave$fullrep <- fullrep
      BigSave$parameters <- parameters
      BigSave$pin  <- pout
      BigSave$best <- best
      BigSave$Gradient <- Grad
      #BigSave$lowlike <- model$fn(best)
      if(max(list.files()=='Output')==1) { setwd(paste(getwd(), "/Output",sep=""))  }
      save(BigSave,file="BigSave.lda")
      print("making Output.RL")
      WriteOutput(Report,SDrep,fullrep,parameters,pout,GeneralSpecs,ControlSpecs,TheData,CurrPhase=0,best=best,grad=Grad)
      }
    }
  }



#' Load Initial Parameter Values for Model Estimation
#'
#' Reads initial parameter values, bounds, and estimation phases from all
#' model input files and prepares them for the optimization routine. This
#' function must be called before running the model estimation.
#'
#' @param aask Character string for testing mode. Use 'test' to display a
#'   console summary of parameter inputs. Default is '' (no testing output).
#'
#' @return NULL. Creates global objects InitialVars and parameters in the
#'   parent environment. Sets global variables MaxPhase, ParOld, and CurrPhase.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reads parameter specifications from all .DAT input files
#'   \item Checks data integrity with isnafunc()
#'   \item Initializes parameter list with proper structure for TMB
#'   \item Handles special cases (e.g., single area models)
#'   \item Records which parameters are active for output tracking
#' }
#'
#' Parameter groups loaded include:
#' MainPars, RecruitPars, PuerPowPars, SelPars, RetPars, RecDevs, Qpars,
#' efpars, InitPars, RecSpatDevs, MovePars, and GrowthPars.
#'
#' @note This function modifies global variables and must be run before
#'   model estimation. Requires all input .DAT files to be present in the
#'   working directory.
#'
#' @examples
#' \dontrun{
#' # Standard workflow
#' choose_model()
#' LoadPars()
#'
#' # Test mode with console output
#' LoadPars(aask = 'test')
#' }
#'
#' @seealso
#' \code{\link{AdjustPhase}} for modifying parameter estimation phases,
#' \code{\link{choose_model}} for selecting model directory
#'
#' @export
LoadPars <- function(aask=''){
  #for(i in 1:length(Data)){    isnafunc(Data[[i]],i)   }
  outtmp <- isnafunc2()
  if(!is.null(outtmp[[1]]))   { warning("\nThere are some NA's in your data: ", paste(outtmp[[1]], collapse = ', '), '\n', call. = FALSE) }
  InitialVars <<- ReadInitialValues(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
  if(aask=='test')dlg_message("Check Console for summary of parameter inputs")
  if(Data$Narea==1){## need to trick SetInitialAndPhases because only one area
    InitialVars$RecSpatDevs$Initial <<- 0
    InitialVars$RecSpatDevs$Bnd <<- c(-15,15)
    InitialVars$RecSpatDevs$Phase <<- -1
  }
  Parssolved(InitialVars) ## Records which parameters were used to solve the model for output file
  parameters <- list(MainPars=NULL,RecruitPars=NULL,PuerPowPars=NULL,SelPars=NULL,RetPars=NULL,RecDevs=NULL,Qpars=NULL,efpars=NULL,InitPars=NULL,RecSpatDevs=NULL,MovePars=NULL,GrowthPars=NULL,dummy=0)
  (MaxPhase <<- getPhase(InitialVars));ParOld <<- NULL;CurrPhase <<- 1
}

#' Generate Model Diagnostics Report
#'
#' Creates a comprehensive diagnostics report from a completed model run,
#' including plots, tables, and summary statistics for model assessment.
#'
#' @param is95 Logical. If TRUE (default), uses 95% confidence intervals
#'   in output plots and tables. If FALSE, may use alternative interval levels.
#'
#' @return NULL. Creates output files in the Output/ directory as a side effect.
#'
#' @details
#' The function automatically navigates to the Output directory (or creates it)
#' and generates diagnostic materials including:
#' \itemize{
#'   \item Model fit statistics and residual plots
#'   \item Parameter estimates with confidence intervals
#'   \item Time series of estimated quantities (biomass, recruitment, etc.)
#'   \item Length frequency fit diagnostics
#'   \item CPUE index fits
#' }
#'
#' @note Must be run from a model directory containing completed model output.
#'   The function calls MakeOutPut() which generates the actual report files.
#'
#' @examples
#' \dontrun{
#' # After running model estimation
#' choose_model()
#' # ... run model ...
#'
#' # Generate diagnostics with 95% CI
#' MakeDiagReport()
#'
#' # Generate diagnostics with alternative CI
#' MakeDiagReport(is95 = FALSE)
#' }
#'
#' @seealso \code{\link{choose_model}} for selecting model directory
#'
#' @export
MakeDiagReport <- function(is95=T) {
  ## Run and output diagnostics file
  print("Making Diagnostics report")
  if(max(list.files()%in%'Output')==1) {  setwd(makehtml::filenametopath(getwd(),'Output'))}
  #source('../../R files/MakeOutPut.R')
  MakeOutPut(is95)
}

#' Interactively Adjust Parameter Estimation Phases
#'
#' Provides an interactive interface to modify which parameters are estimated
#' in which phases of the optimization. Phases control the sequential estimation
#' of parameter groups, with lower phases estimated before higher phases.
#'
#' @param dum Character string. Use 'dummy' to set all parameters to phase -1
#'   (fixed) except growth parameters. Default '' provides interactive dialog.
#'
#' @return NULL. Modifies global objects InitialVars and MaxPhase.
#'
#' @details
#' Estimation phases allow sequential parameter estimation:
#' \itemize{
#'   \item Phase -1: Parameter is fixed at initial value (not estimated)
#'   \item Phase 1: Estimated in first optimization phase
#'   \item Phase 2+: Estimated after lower phases converge
#' }
#'
#' The function offers two interactive options for each parameter group:
#' \itemize{
#'   \item 'All': Set all parameters in the group to the same phase
#'   \item 'Individual': Set each parameter's phase separately
#'   \item 'Skip': Leave current phases unchanged
#' }
#'
#' Using 'dummy' mode is useful for testing model structure without full
#' estimation (e.g., checking growth specifications only).
#'
#' @note This function modifies global variables. Run after LoadPars() and
#'   before model estimation to customize the estimation sequence.
#'
#' @examples
#' \dontrun{
#' # Interactive phase adjustment
#' choose_model()
#' LoadPars()
#' AdjustPhase()  # Opens interactive dialogs
#'
#' # Dummy mode - test growth only
#' AdjustPhase(dum = 'dummy')
#' }
#'
#' @seealso
#' \code{\link{LoadPars}} for loading initial parameters,
#' \code{\link{choose_model}} for selecting model directory
#'
#' @export
AdjustPhase <- function(dum=' '){
  MaxPhase <- 0
  for(i in 1:length(names(InitialVars))) {
    if(dum!='dummy') {
      todo <- dlg_list(c('All','Individual','Skip','                '),  title=paste('Phase for',names(InitialVars)[i],'                             ') )$res
      ## Get estimated parameters  KeyWord <- locs$id[i]
      if(todo=='All'){
        nphase <- dlg_list(c('-1','1','2','3','4','                '),  title=paste('Phase for all',names(InitialVars)[i],'                             ') )$res
        InitialVars[[i]]$Phase <- rep(as.numeric(nphase), length(InitialVars[[i]]$Phase))
      }
      if(todo=='Individual'){
        for(ip in 1:length(InitialVars[[i]]$Initial)){
          nphase <- dlg_list(c('-1','1','2','3','4','                '),  title=paste('Phase',names(InitialVars)[i],'par #',ip,'                             ') )$res
          InitialVars[[i]]$Phase[ip] <- as.numeric(nphase)
        }
      }
      if(max(InitialVars[[i]]$Phase)>MaxPhase) MaxPhase <- max(InitialVars[[i]]$Phase)
    }
    if(dum=='dummy') {
      InitialVars[[i]]$Phase <- rep(-1, length(InitialVars[[i]]$Phase))
      tmp <- InitialVars$Growth; tmp$Phase <- 1
      InitialVars$dummy <- tmp
    }
  }
  InitialVars <<- InitialVars
  MaxPhase <<- MaxPhase
}

#' Update Length Frequency Data Weights Based on Tuning Results
#'
#' Updates the weighting of length frequency data in CONTROL.DAT based on
#' model tuning results. This adjusts relative weights to improve model fit
#' to length composition data after iterative tuning analyses.
#'
#' @param todo Character string indicating whether to update weights.
#'   Options are 'Yes' or 'No' (default). Use 'Yes' to apply tuning results.
#'
#' @return NULL (invisibly). The function modifies CONTROL.DAT as a side effect.
#'
#' @details
#' The function reads tuning scale factors from 'Output/Summary/result/Tuning.csv'
#' and applies them to length frequency weights in CONTROL.DAT. Each existing
#' weight is multiplied by its corresponding tuning scale factor from the
#' iterative weighting analysis.
#'
#' The tuning process typically involves:
#' \enumerate{
#'   \item Running the model with initial weights
#'   \item Analyzing length frequency fit diagnostics
#'   \item Calculating appropriate scaling factors
#'   \item Applying these factors using this function
#'   \item Re-running the model with updated weights
#' }
#'
#' The function will stop with an error message if the number of predetermined
#' weights doesn't match the number of length composition datasets.
#'
#' @note
#' \itemize{
#'   \item Must be run from within a model run directory containing CONTROL.DAT
#'   \item Requires 'Output/Summary/result/Tuning.csv' from a completed tuning run
#'   \item The Tuning.csv file must contain columns: Fleet, Sex, and Multiscale
#' }
#'
#' @examples
#' \dontrun{
#' # Navigate to model run directory
#' choose_model()
#'
#' # After running tuning analysis, update weights
#' UpdateLFWeights("Yes")
#'
#' # Skip updating (default behavior)
#' UpdateLFWeights()
#' UpdateLFWeights("No")
#' }
#'
#' @seealso
#' \code{\link{UpdatePars}} for updating model parameters,
#' \code{\link{choose_model}} for selecting a model run directory
#'
#' @export
UpdateLFWeights <- function(todo='No'){
  library(dplyr); library(magrittr)
  if(todo=='Yes'){
    tmp  <- read.csv(paste0(getwd(),"/Output/Summary/result/Tuning.csv"))
    # open Control file and find weightings
    DataFile <- read.table('CONTROL.DAT',comment.char = "?",fill=T, blank.lines.skip=F,stringsAsFactors=F,col.names=1:100)
    ## Find the length freqs weights
    pos1 <- which(1==(grepl('length',DataFile[,8])))
    if(length(pos1)!=nrow(tmp)) {
      print("Predetermined weights do not match length compoitions")
      break  }
    for(i in 1:length(pos1)){
      ttmp <- DataFile[pos1[i],1:10]
      Scale <- tmp$Multiscale[(tmp$Fleet-1)==as.numeric(ttmp[,2]) & (tmp$Sex-1)==as.numeric(ttmp[,4])]
      DataFile[pos1[i],5] <-  as.numeric(DataFile[pos1[i],5]) * as.numeric(Scale)
    }
    write.table(DataFile, 'CONTROL.DAT',na=" ", sep=" ", row.names = F, col.names = F, quote=F)
  }
  }

