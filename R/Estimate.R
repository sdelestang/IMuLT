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
    mout<-nlminb(model$par,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL, eval.max=MaXeVaL, rel.tol=1e-12))
    initBestFn <- BestFn
    pars <- mout$par; names(pars) <- pnames; ParOld <- mout$par;
    Grad <- abs(model$gr(mout$par))
    badpar <- paste0("[",pnames[Grad==max(Grad)],"]")
    cat("Likelihood: ",round(initBestFn,6),' to ' ,round(mout$objective,6),"| Convergence:",ifelse(mout$convergence==0,'Yes','No')," (",mout$convergence,") ","| Max Gradient [Par]:",round(max(abs(model$gr(mout$par))),6), badpar, "| Interations:",mout$iterations,"| Evalutions:",mout$evaluations,"\n")

    if(CurrPhase==MaxPhase) {
      cat("Re-run last phase to further reduce the gradient")
      tmppars <- model$env$last.par.best
      mout<-nlminb(start=tmppars,model$fn,model$gr,lower=RunSpecs$lowBnd,upper=RunSpecs$uppBnd,control = list(iter.max = MaXeVaL/2, eval.max=MaXeVaL, rel.tol=1e-12))
      initBestFn <- BestFn
      pars <- mout$par; names(pars) <- pnames; ParOld <- mout$par;
      Grad <- abs(model$gr(mout$par))
      badpar <- paste0("[",pnames[Grad==max(Grad)],"]")
      cat("Likelihood: ",round(initBestFn,6),' to ' ,round(mout$objective,6),"| Convergence:",ifelse(mout$convergence==0,'Yes','No')," (",mout$convergence,") ","| Max Gradient [Par]:",round(max(abs(model$gr(mout$par))),6), badpar, "| Interations:",mout$iterations,"| Evalutions:",mout$evaluations,"\n")
      }
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
#' @param folder_name The name for the folder to contain the outputs (inside Summary).
#' If omitted or left blank it will revert to the default behaviour of storing outputs
#' in 'summary/result/'.
#'
#'
#' @param openfile Whether to open the html file on completion. Default TRUE
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
MakeDiagReport <- function(is95=T,folder_name = '',openfile=TRUE) {
  ## Run and output diagnostics file
  print("Making Diagnostics report")
  current_wd = getwd()
  if(max(list.files()%in%'Output')==1) {  setwd(makehtml::filenametopath(getwd(),'Output'))}
  #source('../../R files/MakeOutPut.R')
  MakeOutPut(is95,folder_name,openfile=openfile)

    #return user to wd
  setwd(current_wd)
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
  Innames <- c("MainPars","RecruitPars","PuerPowPars","SelPars","RecDevs","RecSpatDevs","efpars","MovePars")
  for(i in 1:length(Innames)) {

    iv <- which(names(InitialVars)==Innames[i]) ## which Initivals vars matches the name
    if(dum!='dummy') {
      #todo <- dlg_list(c('All','Individual','Skip this par','End all resets','                '),  title=paste('Phase for',names(InitialVars)[i],'                             ') )$res

      todo <- tk_choice(
        choices = c('All', 'Individual', 'Skip this par', 'End all resets'),
        title   = paste('Phase for', Innames[i]) )

      if(todo=='End all resets') {
        InitialVars <<- InitialVars
        MaxPhase <<- MaxPhase
        return(invisible(NULL))
      }

      if(todo=='All'){
        #nphase <- dlg_list(c('-1','1','2','3','4','                '),  title=paste('Phase for all',names(InitialVars)[i],'                             ') )$res
        nphase <- tk_choice(
          choices = c('-1','1','2','3','4','5','6'),
          title   = paste('Phase for all', Innames[i]) )

        InitialVars[[iv]]$Phase <- rep(as.numeric(nphase), length(InitialVars[[iv]]$Phase))
      }
      if(todo=='Individual'){
        for(ip in 1:length(InitialVars[[iv]]$Initial)){
          nphase <- tk_choice(
            choices = c('-1','1','2','3','4','5','6','End this par'),
            title   = paste('Phase',Innames[i],'par #',ip) )
          if(nphase=='End this par') break
          InitialVars[[iv]]$Phase[ip] <- as.numeric(nphase)
        }
      }
      if(max(InitialVars[[iv]]$Phase)>MaxPhase) MaxPhase <- max(InitialVars[[i]]$Phase)
    }
    if(dum=='dummy') {
      InitialVars[[i]]$Phase <- rep(-1, length(InitialVars[[i]]$Phase))
      tmp <- InitialVars$Growth; tmp$Phase <- 1
      InitialVars$dummy <- tmp
    }
  }
  (MaxPhase <<- getPhase(InitialVars));ParOld <<- NULL;CurrPhase <<- 1
  InitialVars <<- InitialVars
  MaxPhase <<- MaxPhase
  Parssolved(InitialVars) ## Records which parameters were used to solve the model for output file
}

#' Display a styled Tk list selection dialog
#'
#' A styled replacement for \code{svDialogs::dlg_list} using a light blue
#' Tk window with a listbox and OK button.
#'
#' @param choices Character vector of options to display in the listbox.
#' @param title Character string for the window title. Default is \code{"Select"}.
#'
#' @return A single character string corresponding to the selected item.
#'
#' @import tcltk
#'
#' @export
tk_choice <- function(choices, title = "Select") {
  tt <- tktoplevel(background = "#d6eaf8")        # light blue background
  tkwm.title(tt, title)
  tkwm.geometry(tt, "350x300")
  tkwm.resizable(tt, FALSE, FALSE)
  result <- tclVar("")
  # Title label
  lbl <- tklabel(tt, text = title, font = "Arial 11 bold",
                 background = "#d6eaf8", foreground = "#1a5276")
  tkpack(lbl, pady = c(15, 5))
  # Listbox with blue tones
  lb <- tklistbox(tt, height = length(choices), width = 40,
                  selectmode = "single", font = "Arial 12",
                  background = "#eaf4fb",           # very light blue listbox
                  foreground = "#1a5276",            # dark blue text
                  selectbackground = "#2e86c1",      # mid blue selection
                  selectforeground = "white",
                  borderwidth = 0, relief = "flat",
                  highlightthickness = 1,
                  highlightbackground = "#aed6f1")
  for (item in choices) tkinsert(lb, "end", item)
  tkselection.set(lb, 0)
  tkpack(lb, padx = 20, pady = 10)
  onOK <- function() {
    idx <- as.integer(tkcurselection(lb))
    tclvalue(result) <- choices[idx + 1]
    tkdestroy(tt)
  }
  btn <- tkbutton(tt, text = "OK", width = 15, font = "Arial 11 bold",
                  background = "#2e86c1",            # blue button
                  foreground = "white",
                  activebackground = "#1a5276",       # darker on hover
                  activeforeground = "white",
                  relief = "flat", borderwidth = 0,
                  command = onOK)
  tkpack(btn, pady = 15)
  # --- keyboard bindings ---
  tkbind(tt, "<Return>", onOK)
  tkbind(lb, "<Return>", onOK)
  for (k in seq_along(choices)) {
    local({
      idx <- k - 1
      tkbind(tt, as.character(k), function() {
        tkselection.clear(lb, 0, "end")
        tkselection.set(lb, idx)
      })
    })
  }

  tkfocus(tt)
  tkwait.window(tt)
  return(tclvalue(result))
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
    pos1 <- which(grepl('Weights',DataFile[,2]) & grepl('by',DataFile[,3]))
    pos2 <- which(grepl('Basic',DataFile[,2]) & grepl('parameters',DataFile[,3]))
    pos3 <- which(grepl('3',DataFile[,1]))
    posall <- pos3[pos3>pos1 & pos3<pos2]
    if(length(posall)!=nrow(tmp)) {
      stop("Predetermined weights do not match length compositions")  # Changed!
    }
    for(i in 1:length(posall)){
      ttmp <- DataFile[posall[i],1:10]
      Scale <- tmp$Multiscale[(tmp$Fleet-1)==as.numeric(ttmp[,2]) & (tmp$Sex-1)==as.numeric(ttmp[,4])]
      DataFile[posall[i],5] <-  as.numeric(DataFile[posall[i],5]) * as.numeric(Scale)
    }
    write.table(DataFile, 'CONTROL.DAT',na=" ", sep=" ", row.names = F, col.names = F, quote=F)
  }
}



# ── Helper: build the live-trace plot ─────────────────────────────────────────

#' Build Live NLL Trace Plot
#'
#' Constructs a multi-panel base-graphics figure showing optimisation progress.
#' The top panel shows the full history of the negative log-likelihood across
#' all phases and restarts. Below it, any optimisation stage that recorded two
#' or more print points gets its own sub-panel, arranged side by side in
#' chronological order. Uses \code{\link[graphics]{layout}} for the panel
#' structure.
#'
#' @param df Data frame with columns \code{eval} (cumulative function
#'   evaluations), \code{nll} (negative log-likelihood), and \code{stage}
#'   (character label for the optimisation stage).
#' @param current_nll Numeric. Current best negative log-likelihood, used in
#'   the global panel title.
#'
#' @return Called for its side effect (drawing to the active graphics device).
#'   Returns \code{invisible(NULL)}.
#'
#' @keywords internal
.plot_trace <- function(df, current_nll) {
  # Which stages have >= 2 points?
  stage_order  <- unique(df$stage)
  stage_counts <- table(df$stage)
  plot_stages  <- stage_order[stage_counts[stage_order] >= 2]
  n_sub        <- length(plot_stages)

  n_panels <- 1 + n_sub
  # Grid with equal-sized cells: fill row-wise
  ncol <- min(n_panels, 2)
  nrow <- ceiling(n_panels / ncol)
  par(mfrow = c(nrow, ncol), mar = c(4, 5, 2, 1))

  # Global panel
  plot(df$eval, df$nll, type = "l", lwd = 2, col = "steelblue",
       xlab = "Function evaluations", ylab = "-log L",
       main = paste0("Global  |  -logL = ", round(current_nll, 2)))
  points(tail(df$eval, 1), tail(df$nll, 1), pch = 19, col = "red", cex = 1.5)

  # Stage panels
  for (s in plot_stages) {
    sub <- df[df$stage == s, ]
    plot(sub$eval, sub$nll, type = "l", lwd = 2, col = "steelblue",
         xlab = "Function evaluations", ylab = "-log L", main = s)
    points(tail(sub$eval, 1), tail(sub$nll, 1), pch = 19, col = "red", cex = 1.5)
  }
}


# ── Main function ─────────────────────────────────────────────────────────────

#' Fit the IMuLT Stock Assessment Model
#'
#' Fits the IMuLT TMB model using a phased optimisation approach. Parameters are
#' progressively introduced across phases, with the final phase employing
#' sandwich restarts (alternating L-BFGS-B and nlminb) to escape local minima.
#' Optional Newton polishing steps can further refine the solution.
#'
#' @param phit Integer. Maximum number of function evaluations per phase for
#'   all phases except the last. Default 500.
#' @param lphit Integer. Maximum number of function evaluations for the last
#'   phase (per optimizer call). Default 1000.
#' @param mxph Integer. Maximum phase number. If 0, treated as 1. Defaults to
#'   the global \code{MaxPhase}.
#' @param PrintLag Integer. Print and plot progress every \code{PrintLag}
#'   function evaluations. Default 50.
#' @param report Logical. If \code{TRUE}, produce SD report and save outputs
#'   after the final phase. Default \code{FALSE}.
#' @param nRestarts Integer. Number of sandwich restarts (L-BFGS-B then nlminb)
#'   in the final phase. Default 3.
#' @param newtonSteps Integer. Number of Newton polishing steps after sandwich
#'   restarts. Default 0 (disabled).
#' @param PrintNll Logical. If \code{TRUE}, display a live multi-panel trace
#'   plot of the negative log-likelihood during fitting. Default \code{TRUE}.
#'
#' @details
#' \strong{Phased estimation:} Parameters are activated across phases via
#' \code{SetInitialAndPhases()}. Early phases use \code{phit} evaluations;
#' the final phase uses \code{lphit}.
#'
#' \strong{Sandwich restarts:} In the final phase, after the initial nlminb run,
#' the optimizer alternates between L-BFGS-B (which uses a different Hessian
#' approximation) and nlminb. This helps escape ridges and saddle points that
#' trap a single algorithm. Each restart resets progress tracking for clean
#' reporting.
#'
#' \strong{Newton polishing:} If \code{newtonSteps > 0}, exact Newton steps
#' using the full Hessian are attempted after the sandwich restarts, followed
#' by a final nlminb call.
#'
#' \strong{Live trace plot:} When \code{PrintNll = TRUE}, a multi-panel
#' base-graphics figure is updated every \code{PrintLag} function evaluations
#' using \code{\link[grDevices]{dev.hold}}/\code{\link[grDevices]{dev.flush}}
#' for flicker-free rendering. The top panel shows the full optimisation
#' history across all phases and restarts. Below it, any stage (phase or
#' restart) that accumulates two or more print points receives its own
#' sub-panel, arranged side by side in chronological order. A red dot marks
#' the current best value.
#'
#' \strong{Prior penalties:} The model supports normal (type 1) and gamma
#' (type 2) priors on main parameters via \code{MparsPrior}. Gamma priors
#' are recommended for strictly positive parameters such as natural mortality.
#'
#' \strong{Global side effects:} The function writes to the global environment
#' via \code{<<-}. The trace is stored in \code{TraceDF} (a data frame with
#' columns \code{eval}, \code{nll}, \code{stage}), cumulative evaluations in
#' \code{TotalEval}, and the current stage label in \code{CurrentStage}.
#'
#' @return Invisibly returns \code{NULL}. Side effects include writing
#'   parameter files to \code{Output/}, and if \code{report = TRUE}, saving
#'   \code{BigSave.lda} and calling \code{WriteOutput()}.
#'
#' @examples
#' \dontrun{
#' # Quick fit with defaults
#' FitModel(500, 1000)
#'
#' # Full fit with reporting, extra restarts, and Newton polishing
#' FitModel(500, 3000, report = TRUE, nRestarts = 5, newtonSteps = 3)
#'
#' # Fit without live plotting
#' FitModel(500, 1000, PrintNll = FALSE)
#' }
#'
#' @export
FitModel <- function(phit = 500, lphit = 1000, mxph = MaxPhase,
                     PrintLag = 50, report = FALSE,
                     nRestarts = 0, newtonSteps = 0, PrintNll = TRUE) {

  MaxPhase <- ifelse(mxph == 0, 1, mxph)

  # Initialise trace bookkeeping
  TraceDF      <<- data.frame(eval = numeric(0), nll = numeric(0),
                              stage = character(0),
                              stringsAsFactors = FALSE)
  TotalEval    <<- 0
  CurrentStage <<- ""

  # ---- Local helper to append to the trace ----
  .trace_append <- function(ev, nll) {
    TraceDF <<- rbind(TraceDF,
                      data.frame(eval = ev, nll = nll,
                                 stage = CurrentStage,
                                 stringsAsFactors = FALSE))
  }

  for (CurrPhase in 1:MaxPhase) {

    MaXeVaL <- ifelse(CurrPhase < MaxPhase, phit, lphit)

    parameters <- list(
      MainPars    = InitialVars$MainPars$Initial,
      RecruitPars = InitialVars$RecruitPars$Initial,
      PuerPowPars = InitialVars$PuerPowPars$Initial,
      SelPars     = InitialVars$SelPars$Initial,
      RetPars     = InitialVars$RetPars$Initial,
      RecDevs     = InitialVars$RecDevs$Initial,
      Qpars       = InitialVars$Qpars$Initial,
      efpars      = InitialVars$efpars$Initial,
      InitPars    = InitialVars$InitPars$Initial,
      RecSpatDevs = InitialVars$RecSpatDevs$Initial,
      MovePars    = InitialVars$MovePars$Initial,
      GrowthPars  = InitialVars$GrowthPars$Initial,
      dummy       = 0
    )

    # Set parameters and mapping
    RunSpecs <- SetInitialAndPhases(ParOld, parameters, InitialVars,
                                    CurrPhase = CurrPhase)

    ## Identify active parameters
    pnames <- names(unlist(RunSpecs$map)[!is.na(unlist(RunSpecs$map))])
    nam    <- stringr::str_extract(pnames, "[\\p{Letter}]+")
    num    <- stringr::str_extract(pnames, "\\d+$")
    unnam  <- nam[!duplicated(nam)]

    cat("Making model object that will solve for",
        sum(!is.na(unlist(RunSpecs$map))),
        "parameters.", "Phase =", CurrPhase, "\n")
    for (iii in seq_along(unnam)) {
      print(paste(unnam[iii], length(num[nam == unnam[iii]]), "parameters"))
    }

    ## Build AD model
    model <- MakeADFun(Data, parameters, map = RunSpecs$map,
                       DLL = "IMuLT", silent = TRUE)

    BestFn      <- model$fn()
    initBestFn  <- BestFn
    LastPrintFn <<- BestFn
    FnCallNo    <<- 0
    model$fn_Orig <- model$fn

    # Record starting point
    CurrentStage <<- paste0("Phase ", CurrPhase, " \u2013 Initial")
    .trace_append(TotalEval, BestFn)

    yy <- 1e+10

    # ── Wrapper to track progress and update live plot ────────
    model$fn <- function(x) {
      tyy <- model$fn_Orig(x)
      yy <<- ifelse(is.na(tyy), yy, tyy)
      FnCallNo <<- FnCallNo + 1

      if (BestFn > yy) {
        BestFn <<- yy
        if ((FnCallNo %% PrintLag) == 0) {
          delta <- 100 * (1 - (yy / LastPrintFn))
          cat("Phase ", CurrPhase, " ", FnCallNo, " -LogLike: ",
              round(yy, 3), " | Delta: ", round(delta, 6), "%\n", sep = "")
          LastPrintFn <<- yy

          # Append to trace
          .trace_append(TotalEval + FnCallNo, yy)

          # Live plot
          if (PrintNll) {
            dev.hold()
            .plot_trace(TraceDF, yy)
            dev.flush()
          }
        }
      }
      return(tyy)
    }

    model$par <- RunSpecs$EstVec

    # ---- Control list used throughout ----
    ctrl <- list(iter.max = MaXeVaL, eval.max = MaXeVaL,
                 rel.tol = 1e-12, x.tol = 1e-12, abs.tol = 0)

    has_bounds <- !is.null(RunSpecs$lowBnd) && !is.null(RunSpecs$uppBnd)

    # ===========================================================
    # Initial nlminb run
    # ===========================================================
    CurrentStage <<- paste0("Phase ", CurrPhase, " \u2013 nlminb")
    BestFn       <- model$fn(model$par)
    initBestFn   <- BestFn
    LastPrintFn  <<- BestFn
    FnCallNo     <<- 0

    .trace_append(TotalEval, BestFn)

    mout <- nlminb(model$par, model$fn, model$gr,
                   lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                   control = ctrl)

    .report_fit(mout, model, pnames, initBestFn, label = "  Initial nlminb")
    TotalEval <<- TotalEval + FnCallNo

    # ===========================================================
    # Sandwich restarts (final phase only)
    # ===========================================================
    if (CurrPhase == MaxPhase) {

      for (restart in seq_len(nRestarts)) {

        cat("\n--- Sandwich restart", restart, "of", nRestarts, "---\n")

        # ---- Step A: L-BFGS-B ────────────────────────────────
        CurrentStage <<- paste0("Restart ", restart, " \u2013 L-BFGS-B")
        cat("  Step A: L-BFGS-B\n")
        FnCallNo     <<- 0
        BestFn       <- model$fn(model$env$last.par.best)
        LastPrintFn  <<- BestFn
        bfgs_start   <- model$env$last.par.best

        .trace_append(TotalEval, BestFn)

        if (has_bounds) {
          fit_bfgs <- optim(bfgs_start, model$fn, model$gr,
                            method  = "L-BFGS-B",
                            lower   = RunSpecs$lowBnd,
                            upper   = RunSpecs$uppBnd,
                            control = list(maxit = MaXeVaL, factr = 1e-15))
        } else {
          fit_bfgs <- optim(bfgs_start, model$fn, model$gr,
                            method  = "BFGS",
                            control = list(maxit = MaXeVaL, reltol = 1e-12))
        }

        bfgs_grad <- max(abs(model$gr(fit_bfgs$par)))
        cat("  L-BFGS-B complete: obj =", round(fit_bfgs$value, 6),
            "| max|grad| =", round(bfgs_grad, 6),
            "| convergence:", fit_bfgs$convergence, "\n")
        TotalEval <<- TotalEval + FnCallNo

        # ---- Step B: nlminb from BFGS solution ───────────────
        CurrentStage <<- paste0("Restart ", restart, " \u2013 nlminb")
        cat("  Step B: nlminb\n")
        FnCallNo     <<- 0
        BestFn       <- fit_bfgs$value
        initBestFn   <- BestFn
        LastPrintFn  <<- BestFn

        .trace_append(TotalEval, BestFn)

        mout <- nlminb(fit_bfgs$par, model$fn, model$gr,
                       lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                       control = ctrl)

        .report_fit(mout, model, pnames, initBestFn,
                    label = paste("  nlminb restart", restart))
        TotalEval <<- TotalEval + FnCallNo

        # ---- Early exit if gradient is small enough ----
        cur_grad <- max(abs(model$gr(mout$par)))
        cat("  Restart", restart, "complete: max|grad| =",
            round(cur_grad, 6), "\n")
        if (cur_grad < 1e-3) {
          cat("  Gradient < 1e-3 \u2014 exiting restart loop early.\n")
          break
        }
      }

      # ===========================================================
      # Newton polishing steps
      # ===========================================================
      if (newtonSteps > 0) {
        cat("\n--- Newton polishing steps ---\n")
        CurrentStage <<- "Newton polish"
        newton_par   <- model$env$last.par.best

        for (ns in seq_len(newtonSteps)) {
          tryCatch({
            H    <- optimHess(newton_par, model$fn, model$gr)
            g    <- as.vector(model$gr(newton_par))
            step <- solve(H, g)
            newton_par <- newton_par - step

            # Respect bounds if they exist
            if (has_bounds) {
              newton_par <- pmax(newton_par, RunSpecs$lowBnd)
              newton_par <- pmin(newton_par, RunSpecs$uppBnd)
            }

            ng <- max(abs(model$gr(newton_par)))
            cat("  Newton step", ns, "- obj:",
                round(model$fn(newton_par), 6),
                "| max|grad|:", round(ng, 8), "\n")

            if (ng < 1e-3) {
              cat("  Gradient < 1e-3 after Newton \u2014 stopping.\n")
              break
            }
          }, error = function(e) {
            cat("  Newton step", ns,
                "failed (Hessian singular?):", conditionMessage(e), "\n")
          })
        }

        # Final nlminb from Newton-polished parameters
        CurrentStage <<- "Post-Newton nlminb"
        cat("  Final nlminb (post-Newton)\n")
        FnCallNo     <<- 0
        BestFn       <- model$fn(newton_par)
        initBestFn   <- BestFn
        LastPrintFn  <<- BestFn

        .trace_append(TotalEval, BestFn)

        mout <- nlminb(newton_par, model$fn, model$gr,
                       lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                       control = ctrl)
        .report_fit(mout, model, pnames, initBestFn,
                    label = "  Final nlminb (post-Newton)")
        TotalEval <<- TotalEval + FnCallNo
      }

      # ===========================================================
      # Bound diagnostics
      # ===========================================================
      .check_bounds(mout$par, RunSpecs$lowBnd, RunSpecs$uppBnd, pnames, model)
    }

    # ===========================================================
    # Store and save parameters
    # ===========================================================
    pars <- mout$par
    names(pars) <- pnames
    ParOld <- mout$par

    pout <- unlist(parameters)
    pout[names(pout) %in% names(pars)] <- pars
    suffix <- ifelse(CurrPhase == MaxPhase, " final", CurrPhase)
    write.table(pout, paste0("Output/model", suffix, ".par"),
                sep = "\t", col.names = c("name\test"), quote = FALSE)

    # ===========================================================
    # Report (final phase only)
    # ===========================================================
    if (report && CurrPhase == MaxPhase) {
      cat("Making report object.\n")
      print("Loading report")
      Report <- model$report()
      best   <- mout$par

      print("Loading SD report (can take quite a long time)")
      SDrep   <- sdreport(model)
      fullrep <- summary(SDrep)

      BigSave <- list(
        Report     = Report,
        SDrep      = SDrep,
        map        = RunSpecs$map,
        Data       = Data,
        fullrep    = fullrep,
        parameters = parameters,
        pin        = pout,
        best       = best,
        Gradient   = abs(model$gr(best))
      )

      if (max(list.files() == "Output") == 1) {
        setwd(paste0(getwd(), "/Output"))
      }
      save(BigSave, file = "BigSave.lda")

      print("making Output.RL")
      WriteOutput(Report, SDrep, fullrep, parameters, pout,
                  GeneralSpecs, ControlSpecs, TheData,
                  CurrPhase = 0, best = best,
                  grad = abs(model$gr(best)))
    }
  }
}


#' Report Optimiser Result
#'
#' Prints a summary line for an optimiser run including starting and ending
#' likelihood, convergence status, maximum gradient and the associated
#' parameter name, iteration and evaluation counts.
#'
#' @param mout Output list from \code{\link[stats]{nlminb}}.
#' @param model TMB model object created by \code{\link[TMB]{MakeADFun}}.
#' @param pnames Character vector of active parameter names.
#' @param initBestFn Numeric. Objective value at the start of this optimiser
#'   call.
#' @param label Character. Label prefix for the printed line.
#'
#' @return Called for its side effect (printing to the console). Returns
#'   \code{invisible(NULL)}.
#'
#' @keywords internal
.report_fit <- function(mout, model, pnames, initBestFn, label = "") {
  Grad   <- abs(model$gr(mout$par))
  badpar <- paste0("[", pnames[Grad == max(Grad)], "]")
  cat(label, "- Likelihood:", round(initBestFn, 6), "to", round(mout$objective, 6),
      "| Convergence:", ifelse(mout$convergence == 0, "Yes", "No"),
      "(", mout$convergence, ")",
      "| Max|grad| [par]:", round(max(Grad), 6), badpar,
      "| Iter:", mout$iterations,
      "| Eval:", mout$evaluations, "\n")
}


#' Check for Parameters at Bounds
#'
#' Prints a warning if any estimated parameters are sitting at or near their
#' lower or upper bounds, along with the associated gradient. Parameters at
#' bounds with large gradients indicate the bound is constraining the solution.
#'
#' @param par Numeric vector. Estimated parameter values.
#' @param lower Numeric vector or \code{NULL}. Lower bounds.
#' @param upper Numeric vector or \code{NULL}. Upper bounds.
#' @param pnames Character vector. Parameter names.
#' @param model TMB model object created by \code{\link[TMB]{MakeADFun}}.
#' @param tol Numeric. Tolerance for detecting bound proximity. Default 1e-4.
#'
#' @return Called for its side effect (printing to the console). Returns
#'   \code{invisible(NULL)}.
#'
#' @keywords internal
.check_bounds <- function(par, lower, upper, pnames, model, tol = 1e-4) {
  if (is.null(lower) || is.null(upper)) return(invisible(NULL))

  at_lower <- which(abs(par - lower) < tol)
  at_upper <- which(abs(par - upper) < tol)
  at_bound <- c(at_lower, at_upper)

  if (length(at_bound) > 0) {
    Grad <- abs(model$gr(par))
    cat("\n*** WARNING: Parameters at or near bounds ***\n")
    for (idx in at_bound) {
      side <- ifelse(idx %in% at_lower, "LOWER", "UPPER")
      cat("  ", pnames[idx], "=", round(par[idx], 6),
          " [", side, " bound:", ifelse(side == "LOWER", lower[idx], upper[idx]), "]",
          " |grad| =", round(Grad[idx], 6), "\n")
    }
    cat("  If these have large gradients, widen the bound or fix via map.\n\n")
  } else {
    cat("\n  No parameters at bounds.\n")
  }
}
