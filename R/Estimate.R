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
    ReturnObj$Est <- Est
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
        parameters[[ParName]] <- PhaseOut$Est
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
  selected <- dlg_list(x, title = "Choose a model:")$res
  if (!length(selected)) {
    cat(paste("OK, the default model is", x[1], "\n"))
    selected <- x[1]
  } else {
    cat(paste("Model", selected, "has been chosen\n"))
  }
  setwd(file.path(getwd(), selected))
  invisible(selected)
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
    ## Ensure model final.par exists first
    par_file <- "Output/model final.par"
    if (!file.exists(par_file)) {
      stop("'", par_file, "' not found. The model must be solved first, e.g.:\n",
           "  FitModel(1000, 2000)\n",
           "before parameters can be updated.", call. = FALSE)
    }

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


#' Load Initial Parameter Values for Model Estimation
#'
#' Reads initial parameter values, bounds, and estimation phases from all
#' model input files and prepares them for the optimisation routine. This
#' function must be called before running the model estimation. Prints a
#' parameter summary table to the console showing counts of estimated,
#' linked, and offset parameters per group, followed by a detailed listing
#' of any linked parameters to aid in diagnosing mis-specified links.
#'
#' @param aask Character string for testing mode. Use 'test' to display a
#'   dialog box prompting the user to check the console. Default is ''
#'   (no dialog).
#'
#' @return NULL (invisible). Creates global objects \code{InitialVars} and
#'   \code{parameters} in the parent environment. Sets global variables
#'   \code{MaxPhase}, \code{ParOld}, and \code{CurrPhase}.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Reads parameter specifications from all .DAT input files
#'   \item Checks data integrity with \code{isnafunc2()}
#'   \item Initialises parameter list with proper structure for TMB
#'   \item Handles special cases (e.g., single area models)
#'   \item Records which parameters are active for output tracking
#'   \item Prints a summary table of parameter groups with columns for
#'     Total, Estimated, Phases, Linked (positive links = copy), and
#'     Offset (negative links = additive offset)
#'   \item Lists each linked parameter individually, showing source and
#'     target within the group, to help catch mis-specified link indices
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
#' # Test mode with dialog prompt
#' LoadPars(aask = 'test')
#' }
#'
#' @seealso
#' \code{\link{AdjustPhase}} for modifying parameter estimation phases,
#' \code{\link{choose_model}} for selecting model directory
#'
#' @export
LoadPars <- function(aask=''){
  outtmp <- isnafunc2()
  if(!is.null(outtmp[[1]]))   { warning("\nThere are some NA's in your data: ", paste(outtmp[[1]], collapse = ', '), '\n', call. = FALSE) }
  InitialVars <<- ReadInitialValues(ControlFile,SelexFile,RetainFile,RecruitFile,GrowthFile,MoveFile,GeneralSpecs,ControlSpecs,SelexSpecs,RetenSpecs,GrowthSpecs,MoveSpecs)
  if(aask=='test')dlg_message("Check Console for summary of parameter inputs")
  if(Data$Narea==1){
    InitialVars$RecSpatDevs$Initial <<- 0
    InitialVars$RecSpatDevs$Bnd <<- c(-15,15)
    InitialVars$RecSpatDevs$Phase <<- -1
  }
  Parssolved(InitialVars)

  # ── Bounds-encompass-initial check ──────────────────────────────────────────
  # For every estimated parameter (phase > 0), verify that the initial value
  # sits within [lower, upper]. Catches CTL file errors such as offset parameters
  # with absolute-value bounds (e.g. lower=30 on an offset that starts at 0).
  .check_bounds_initial <- function(grp_name, initial, lower, upper, phase) {
    if (is.null(initial) || is.null(lower) || is.null(upper) || is.null(phase))
      return(NULL)
    # Bounds may be stored as a 2-element vector (global) or per-parameter matrix
    if (length(lower) == 1) lower <- rep(lower, length(initial))
    if (length(upper) == 1) upper <- rep(upper, length(initial))
    # Only check estimated parameters
    est_idx <- which(phase > 0)
    if (length(est_idx) == 0) return(NULL)
    bad <- est_idx[initial[est_idx] < lower[est_idx] |
                     initial[est_idx] > upper[est_idx]]
    if (length(bad) == 0) return(NULL)
    data.frame(
      Group   = grp_name,
      Index   = bad,
      Initial = initial[bad],
      Lower   = lower[bad],
      Upper   = upper[bad],
      stringsAsFactors = FALSE
    )
  }

  bnd_groups <- list(
    MainPars    = InitialVars$MainPars,
    RecruitPars = InitialVars$RecruitPars,
    PuerPowPars = InitialVars$PuerPowPars,
    SelPars     = InitialVars$SelPars,
    RecDevs     = InitialVars$RecDevs,
    efpars      = InitialVars$efpars,
    MovePars    = InitialVars$MovePars
  )

  bnd_issues <- do.call(rbind, lapply(names(bnd_groups), function(grp) {
    g <- bnd_groups[[grp]]
    # Bounds stored as Bnd (2-col matrix or 2-element vector) or Lower/Upper
    if (!is.null(g$Bnd)) {
      if (is.matrix(g$Bnd)) {
        lower <- g$Bnd[, 1]; upper <- g$Bnd[, 2]
      } else {
        lower <- g$Bnd[1];   upper <- g$Bnd[2]
      }
    } else if (!is.null(g$Lower) && !is.null(g$Upper)) {
      lower <- g$Lower;      upper <- g$Upper
    } else {
      return(NULL)
    }
    .check_bounds_initial(grp, g$Initial, lower, upper, g$Phase)
  }))

  if (!is.null(bnd_issues) && nrow(bnd_issues) > 0) {
    cat("\n*** WARNING: Initial values outside bounds ***\n")
    cat("  These parameters will be clamped or cause immediate bound-hitting.\n")
    cat("  Check CTL file — offset/linked parameters likely have wrong bounds.\n\n")
    print(bnd_issues, row.names = FALSE)
    cat("\n")
  } else {
    cat("  Bounds check: all initial values within bounds.\n")
  }
  # ────────────────────────────────────────────────────────────────────────────

  ## Parameter link summary
  par_info <- list(
    MainPars    = list(phase = InitialVars$MainPars$Phase,      link = Data$MparsLink),
    RecruitPars = list(phase = InitialVars$RecruitPars$Phase,   link = Data$RecparsLink),
    PuerPowPars = list(phase = InitialVars$PuerPowPars$Phase,   link = NULL),
    SelPars     = list(phase = InitialVars$SelPars$Phase,       link = Data$SelparsLink),
    #RetPars     = list(phase = InitialVars$RetPars$Phase,       link = NULL),
    RecDevs     = list(phase = InitialVars$RecDevs$Phase,       link = NULL),
    EffPars     = list(phase = InitialVars$efpars$Phase,        link = Data$EffparsLink),
    RecSpatDevs = list(phase = InitialVars$RecSpatDevs$Phase,   link = NULL),
    MovePars    = list(phase = InitialVars$MovePars$Phase,      link = Data$MoveparsLink)#,
    #GrowthPars  = list(phase = InitialVars$GrowthPars$Phase,    link = NULL)
  )
  par_summary <- do.call(rbind, lapply(names(par_info), function(grp) {
    ph <- par_info[[grp]]$phase
    lk <- par_info[[grp]]$link
    if(is.null(ph)) return(NULL)
    if(is.null(lk)) lk <- rep(0, length(ph))
    n_est    <- sum(ph > 0)
    ph_used  <- paste(sort(unique(ph[ph > 0])), collapse = ",")
    if(ph_used == "") ph_used <- "-"
    n_linked <- sum(lk > 0)
    n_offset <- sum(lk < 0)
    data.frame(Group = grp, Total = length(ph), Estimated = n_est,
               Phases = ph_used,
               Linked = ifelse(n_linked == 0, "-", n_linked),
               Offset = ifelse(n_offset == 0, "-", n_offset),
               stringsAsFactors = FALSE)
  }))
  cat("\n--- Parameter Summary ---\n")
  print(par_summary, row.names = FALSE, right = FALSE)
  for(grp in names(par_info)) {
    lk <- par_info[[grp]]$link
    if(is.null(lk)) next
    idx <- which(lk != 0)
    if(length(idx) > 0) {
      for(j in idx) {
        ltype <- ifelse(lk[j] > 0, "link", "link+offset")
        cat(sprintf("  %s_%d -> %s_%d (%s)\n", grp, j, grp, abs(lk[j]), ltype))
      }
    }
  }
  cat("------------------------\n\n")
  parameters <- list(MainPars=NULL,RecruitPars=NULL,PuerPowPars=NULL,SelPars=NULL,RetPars=NULL,RecDevs=NULL,Qpars=NULL,efpars=NULL,InitPars=NULL,RecSpatDevs=NULL,MovePars=NULL,dummy=0)
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

  ## Check pars and data loaded first.
  if (!exists("Data", envir = .GlobalEnv)) {
    stop("Data not loaded. Run LoadData() first.", call. = FALSE)
  }
  if (!exists("InitialVars", envir = .GlobalEnv)) {
    stop("Parameters not loaded. Run LoadData() then LoadPars() first.", call. = FALSE)
  }

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
    ## Check Output.RL exists
    rl_file <- "Output/Output.RL"
    if (!file.exists(rl_file)) {
      message("Output.RL not found. Run FitModel with report = TRUE,..  first, e.g.:\n",
              "  FitModel(1000, 2000, report = TRUE)\n",
              "Skipping length-frequency weight update.")
      return(invisible(NULL))
    }

    ## Function to calculate the slope
    lmslope <- function(x,y)  return(coefficients(lm(y~x-1)))

    output <- read.table(rl_file, comment.char = "?", fill=T, blank.lines.skip=F, stringsAsFactors=F, col.names=1:200)
    tdat <- findNclean(c('Obs/Pred','Fleet'), output, 1, convert=1) %>% filter(`Obs/Pred`=='P')
    tdat <- tdat[,1:7]
    DataFile   <- read.table('CONTROL.DAT', comment.char = "?", fill=T, blank.lines.skip=F, stringsAsFactors=F, col.names=1:200)
    didtune <- findNclean(c('#','Weights','by','fleet'), DataFile, 3)
    tdat$ScaleNsamp <- didtune[match(tdat$Fleet,(didtune$Fleet+1)),5]
    tdat %<>% mutate(ScaleNsamp=ifelse(is.na(ScaleNsamp),1,ScaleNsamp)) %>% group_by(Sex, Fleet) %>% summarise(Multiscale=lmslope(ScaleNsamp,EffN))

    if (identical(tdat, NA) || is.null(tdat) || !is.data.frame(tdat)) {
      message("No Francis tuning data found in Output.RL. Skipping weight update.")
      return(invisible(NULL))
    }

    pos1 <- which(grepl('Weights',DataFile[,2]) & grepl('by',DataFile[,3]))
    pos2 <- which(grepl('Basic',DataFile[,2]) & grepl('parameters',DataFile[,3]))
    pos3 <- which(grepl('3',DataFile[,1]) & nchar(DataFile[,1])==1)
    posall <- pos3[pos3>pos1 & pos3<pos2]
    if(length(posall)!=nrow(tdat)) {
      stop("Predetermined weights do not match length compositions")
    }
    for (i in 1:length(posall)) {
      ttmp <- DataFile[posall[i], 1:10]
      Scale <- as.numeric(tdat$Multiscale[(as.numeric(tdat$Fleet) - 1) == as.numeric(ttmp[, 2]) & (as.numeric(tdat$Sex) - 1) == as.numeric(ttmp[, 4])])
      DataFile[posall[i], 5] <- as.numeric(DataFile[posall[i], 5]) * as.numeric(Scale)
    }
    write.table(DataFile, 'CONTROL.DAT', na=" ", sep=" ", row.names = F, col.names = F, quote=F)
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
  if(n_panels>4) ncol <- min(n_panels, 3)
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
#' sandwich restarts (alternating BFGS on logit-transformed parameters and
#' nlminb) to escape local minima. Optional Newton polishing steps can further
#' refine the solution.
#'
#' Non-final phases include a plateau diagnostic: after each phase completes,
#' the rolling delta history is inspected and a message is printed if a plateau
#' was detected, suggesting a more efficient \code{phit} value for future runs.
#' nlminb is never interrupted — the diagnostic is informational only.
#'
#' @param phit Integer. Maximum number of function evaluations per phase for
#'   all phases except the last. Default 500.
#' @param lphit Integer. Maximum number of function evaluations for the last
#'   phase (per optimiser call). Default 1000.
#' @param mxph Integer. Maximum phase number. If 0, treated as 1. Defaults to
#'   the global \code{MaxPhase}.
#' @param PrintLag Integer. Print and plot progress every \code{PrintLag}
#'   function evaluations. Default 50.
#' @param report Logical. If \code{TRUE}, produce SD report and save outputs
#'   after the final phase. Default \code{FALSE}.
#' @param nRestarts Logical. If \code{TRUE}, sandwich restarts (alternating
#'   BFGS and nlminb) run adaptively in the final phase until either the
#'   gradient falls below \code{newton_grad_thresh} (ready for Newton) or two
#'   successive cycles fail to improve the NLL (stuck). If \code{FALSE}, no
#'   sandwich is performed. Default \code{TRUE}.
#' @param newtonSteps Integer. Maximum number of damped Newton polishing steps
#'   after the sandwich. Only attempted if the gradient is below
#'   \code{newton_grad_thresh} at sandwich exit; otherwise a warning is
#'   printed. Default 5.
#' @param newton_grad_thresh Numeric. Gradient threshold below which Newton
#'   polishing is considered safe. Also used as the sandwich convergence
#'   criterion. Default 0.1.
#' @param PrintNll Logical. If \code{TRUE}, display a live multi-panel trace
#'   plot of the negative log-likelihood during fitting. Default \code{TRUE}.
#' @param plateau_k Integer. Number of consecutive reporting intervals (each
#'   \code{PrintLag} evaluations) used in the rolling window for plateau
#'   detection. Default 4.
#' @param plateau_cv Numeric. Coefficient of variation threshold below which
#'   the delta window is declared a plateau. A diagnostic message is printed
#'   after the phase suggesting a more efficient \code{phit}. Set to 0 to
#'   disable. Default 0.05 (5\%).
#'
#' @details
#' \strong{Phased estimation:} Parameters are activated across phases via
#' \code{SetInitialAndPhases()}. Early phases use \code{phit} evaluations;
#' the final phase uses \code{lphit}.
#'
#' \strong{Plateau diagnostic:} In non-final phases, after the full nlminb call
#' completes, the rolling window of \code{plateau_k} delta values is inspected.
#' If the coefficient of variation falls below \code{plateau_cv}, a message is
#' printed reporting the eval at which the plateau was first detected and
#' suggesting that \code{phit} could be reduced to that value for future runs.
#' nlminb is never interrupted — this is purely informational. Set
#' \code{plateau_cv = 0} to disable. The final phase is never reported on.
#'
#' \strong{Sandwich restarts:} If \code{nRestarts = TRUE}, the final phase
#' runs adaptive sandwich cycles (BFGS then nlminb) until the gradient falls
#' below \code{newton_grad_thresh} or two successive cycles produce no NLL
#' improvement. BFGS operates in logit-transformed unconstrained space so it
#' can never step outside bounds or into non-finite regions. The gradient is
#' adjusted via the chain rule. After BFGS, parameters are back-transformed
#' and passed to nlminb with full bounds. A safety ceiling of 20 cycles
#' prevents runaway loops.
#'
#' \strong{Finite penalty wrapper:} The objective function wrapper returns a
#' large finite penalty (\code{last_good * 1.5 + 1e6}) for any non-finite
#' TMB evaluation, preventing optimiser crashes.
#'
#' \strong{Gradient wrapper:} The gradient wrapper returns a zero vector for
#' any non-finite or error-producing gradient evaluation, preventing crashes
#' in optimisers that require finite gradients.
#'
#' \strong{Newton polishing:} If \code{newtonSteps > 0} and the gradient at
#' sandwich exit is below \code{newton_grad_thresh}, damped Newton steps are
#' attempted using the full Hessian. Each step performs a backtracking line
#' search (halving the step size up to 10 times) to ensure the objective
#' actually improves before accepting the step. A final nlminb call follows
#' to clean up any remaining gradient. If the gradient is above
#' \code{newton_grad_thresh} at sandwich exit, Newton is skipped with a
#' warning — the Hessian is unlikely to be well-conditioned and steps would
#' be unreliable.
#'
#' \strong{Global side effects:} Writes to the global environment via
#' \code{<<-}. The trace is stored in \code{TraceDF}, cumulative evaluations
#' in \code{TotalEval}, and the current stage label in \code{CurrentStage}.
#'
#' @return Invisibly returns \code{NULL}. Side effects include writing
#'   parameter files to \code{Output/}, and if \code{report = TRUE}, saving
#'   \code{BigSave.lda} and calling \code{WriteOutput()}.
#'
#' @examples
#' \dontrun{
#' # Quick fit, no sandwich
#' FitModel(500, 1000, nRestarts = FALSE)
#'
#' # Full fit — adaptive sandwich + Newton polishing
#' FitModel(500, 3000, report = TRUE)
#'
#' # Tighter Newton threshold (only polish if gradient < 0.01)
#' FitModel(500, 3000, report = TRUE, newton_grad_thresh = 0.01)
#'
#' # Fit without live plotting
#' FitModel(500, 1000, PrintNll = FALSE)
#'
#' # Tighter plateau diagnostic
#' FitModel(500, 1000, plateau_k = 5, plateau_cv = 0.02)
#'
#' # Disable plateau diagnostic
#' FitModel(500, 1000, plateau_cv = 0)
#' }
#'
#' @export
FitModel <- function(phit = 500, lphit = 1000, mxph = MaxPhase,
                     PrintLag = 50, report = FALSE,
                     nRestarts = TRUE, newtonSteps = 5,
                     newton_grad_thresh = 0.1,
                     PrintNll = TRUE,
                     plateau_k = 4, plateau_cv = 0.05) {

  MaxPhase <- ifelse(mxph == 0, 1, mxph)

  # ── Enforce negative phase for linked parameters ──────────────────────────
  link_map <- list(MainPars    = Data$MparsLink,
                   RecruitPars = Data$RecparsLink,
                   SelPars     = Data$SelparsLink,
                   efpars      = Data$EffparsLink)
  for (pname in names(link_map)) {
    lv <- link_map[[pname]]
    if (!is.null(lv)) {
      for (i in seq_along(lv)) {
        if (!is.na(lv[i]) && lv[i] > 0 && InitialVars[[pname]]$Phase[i] > 0) {
          warning(paste(pname, "parameter", i,
                        "is linked but has positive phase — forcing negative"))
          InitialVars[[pname]]$Phase[i] <- -abs(InitialVars[[pname]]$Phase[i])
        }
      }
    }
  }

  # ── Initialise trace bookkeeping ──────────────────────────────────────────
  TraceDF      <<- data.frame(eval = numeric(0), nll = numeric(0),
                              stage = character(0),
                              stringsAsFactors = FALSE)
  TotalEval    <<- 0
  CurrentStage <<- ""

  .trace_append <- function(ev, nll) {
    TraceDF <<- rbind(TraceDF,
                      data.frame(eval  = ev,
                                 nll   = nll,
                                 stage = CurrentStage,
                                 stringsAsFactors = FALSE))
  }

  # ── Logit transform helpers (bounded <-> unconstrained) ───────────────────
  # Maps x in [lo, hi] to y in (-Inf, Inf) and back.
  # Keeps x strictly inside bounds to avoid log(0).
  .to_unbounded <- function(x, lo, hi) {
    x_c <- pmax(pmin(x, hi - 1e-8), lo + 1e-8)
    log((x_c - lo) / (hi - x_c))
  }

  .to_bounded <- function(y, lo, hi) {
    lo + (hi - lo) / (1 + exp(-y))
  }

  # Chain-rule correction: dL/dy = dL/dx * dx/dy
  # dx/dy for logit = (hi - lo) * sigmoid(y) * (1 - sigmoid(y))
  .chain_grad <- function(g, y, lo, hi) {
    sig  <- 1 / (1 + exp(-y))
    dxdy <- (hi - lo) * sig * (1 - sig)
    g * dxdy
  }

  # ── Phase loop ────────────────────────────────────────────────────────────
  for (CurrPhase in 1:MaxPhase) {

    MaXeVaL    <- ifelse(CurrPhase < MaxPhase, phit, lphit)
    is_final   <- CurrPhase == MaxPhase

    parameters <- list(
      MainPars    = InitialVars$MainPars$Initial,
      RecruitPars = InitialVars$RecruitPars$Initial,
      PuerPowPars = InitialVars$PuerPowPars$Initial,
      SelPars     = InitialVars$SelPars$Initial,
      RetPars     = InitialVars$RetPars$Initial,
      RecDevs     = InitialVars$RecDevs$Initial,
      Qpars       = InitialVars$Qpars$Initial,
      efpars      = InitialVars$efpars$Initial,
      #InitPars   = InitialVars$InitPars$Initial,
      RecSpatDevs = InitialVars$RecSpatDevs$Initial,
      MovePars    = InitialVars$MovePars$Initial,
      GrowthPars  = InitialVars$GrowthPars$Initial,
      dummy       = 0
    )

    RunSpecs   <- SetInitialAndPhases(ParOld, parameters, InitialVars,
                                      CurrPhase = CurrPhase)
    parameters <- RunSpecs$parameters

    ## Identify active parameters
    pnames <- names(unlist(RunSpecs$map)[!is.na(unlist(RunSpecs$map))])
    nam    <- stringr::str_extract(pnames, "[\\p{Letter}]+")
    num    <- stringr::str_extract(pnames, "\\d+$")
    unnam  <- nam[!duplicated(nam)]

    cat("Making model object that will solve for",
        sum(!is.na(unlist(RunSpecs$map))),
        "parameters. Phase =", CurrPhase, "\n")
    for (iii in seq_along(unnam))
      print(paste(unnam[iii], length(num[nam == unnam[iii]]), "parameters"))

    ## Build AD model
    model <- MakeADFun(Data, parameters, map = RunSpecs$map,
                       DLL = "IMuLT", silent = TRUE)

    model$par  <- RunSpecs$EstVec
    BestFn     <- model$fn()
    if (is.na(BestFn)) BestFn <- Inf
    initBestFn <- BestFn
    LastPrintFn <<- BestFn
    FnCallNo   <<- 0

    CurrentStage <<- paste0("Phase ", CurrPhase, " \u2013 Initial")
    .trace_append(TotalEval, BestFn)

    last_good <- BestFn   # tracks last finite fn value for penalty fallback

    # ── Plateau diagnostic state ──────────────────────────────────────────
    # Tracks delta history across all phases. In non-final phases, after the
    # full nlminb call completes, the history is inspected post-hoc and a
    # suggestion is printed if a plateau was detected. nlminb is never
    # interrupted — this is purely informational.
    delta_history  <- numeric(0)
    plateau_eval   <- NA_integer_   # eval at which plateau was first detected
    min_evals_diag <- PrintLag * (plateau_k + 2)  # warm-up guard

    # ── Print suppression flag ────────────────────────────────────────────
    # Set TRUE during sandwich and Newton to silence per-eval NLL prints.
    # Summary lines from .report_fit are unaffected (outside model$fn).
    suppress_print <- FALSE

    # ── Store originals then wrap both fn and gr ──────────────────────────
    model$fn_Orig <- model$fn
    model$gr_Orig <- model$gr

    # Objective wrapper: finite-penalty + progress tracking + delta accumulation
    model$fn <- function(x) {
      tyy <- model$fn_Orig(x)

      # Return large finite penalty instead of NA/Inf — keeps optimisers alive
      if (is.na(tyy) || !is.finite(tyy)) return(last_good * 1.5 + 1e6)

      last_good <<- tyy
      FnCallNo  <<- FnCallNo + 1

      if (!is.na(BestFn) && BestFn > tyy) {
        BestFn <<- tyy
        if ((FnCallNo %% PrintLag) == 0 && !suppress_print) {
          delta <- 100 * (1 - (tyy / LastPrintFn))
          cat("Phase ", CurrPhase, " ", FnCallNo, " -LogLike: ",
              round(tyy, 3), " | Delta: ", round(abs(delta), 6), "%\n", sep = "")
          LastPrintFn <<- tyy
          .trace_append(TotalEval + FnCallNo, tyy)
          if (PrintNll) {
            dev.hold()
            .plot_trace(TraceDF, tyy)
            dev.flush()
          }

          # ── Accumulate delta for post-hoc plateau diagnostic ───────────
          # Only in non-final phases, past the warm-up guard, and when the
          # feature is enabled. Records the first eval at which the rolling
          # window CV drops below the threshold.
          if (!is_final && plateau_cv > 0 && FnCallNo >= min_evals_diag) {
            delta_history <<- c(delta_history, abs(delta))
            if (length(delta_history) >= plateau_k && is.na(plateau_eval)) {
              recent <- tail(delta_history, plateau_k)
              cv     <- sd(recent) / abs(mean(recent))
              if (cv < plateau_cv) plateau_eval <<- FnCallNo
            }
          }
          # ──────────────────────────────────────────────────────────────
        }
      }
      return(tyy)
    }

    # Gradient wrapper: return zero gradient for non-finite regions instead
    # of crashing optimisers that require finite gradients
    model$gr <- function(x) {
      g <- tryCatch(model$gr_Orig(x), error = function(e) NULL)
      if (is.null(g) || any(!is.finite(g))) return(rep(0, length(x)))
      return(g)
    }
    # ─────────────────────────────────────────────────────────────────────

    has_bounds <- !is.null(RunSpecs$lowBnd) && !is.null(RunSpecs$uppBnd)

    # ── Initial nlminb ────────────────────────────────────────────────────
    CurrentStage <<- paste0("Phase ", CurrPhase, " \u2013 nlminb")
    BestFn       <- model$fn(model$par)
    initBestFn   <- BestFn
    LastPrintFn  <<- BestFn
    FnCallNo     <<- 0
    .trace_append(TotalEval, BestFn)

    ctrl <- list(iter.max = MaXeVaL, eval.max = MaXeVaL,
                 rel.tol = 1e-12, x.tol = 1e-12, abs.tol = 0)
    mout <- nlminb(model$par, model$fn, model$gr,
                   lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                   control = ctrl)

    # ── Post-hoc plateau diagnostic (non-final phases only) ───────────────
    if (!is_final && plateau_cv > 0 && !is.na(plateau_eval)) {
      recent <- tail(delta_history, plateau_k)
      cat("  [Plateau diagnostic] Phase ", CurrPhase,
          " plateaued at eval ~", plateau_eval,
          " (CV = ", round(sd(recent) / abs(mean(recent)), 4),
          ", mean delta = ", round(mean(recent), 4), "%)",
          " — consider phit = ", plateau_eval, " for future runs\n", sep = "")
    }

    .report_fit(mout, model, pnames, initBestFn, label = "  Initial nlminb")
    TotalEval <<- TotalEval + FnCallNo

    # ── Sandwich restarts (final phase only) ──────────────────────────────
    if (is_final && isTRUE(nRestarts)) {

      suppress_print <<- TRUE        # silence per-eval prints during sandwich
      sandwich_done <- FALSE
      no_improve    <- 0L                    # consecutive cycles with no grad improvement
      max_sandwich  <- 20L                   # safety ceiling
      prev_grad     <- max(abs(model$gr_Orig(mout$par)))  # gradient at sandwich entry

      for (restart in seq_len(max_sandwich)) {

        cat("\n--- Sandwich restart", restart, "---\n")

        # Tighter tolerances from restart 2 onwards — by then we're in the
        # basin and want fine-grained polishing rather than broad exploration.
        bfgs_reltol  <- ifelse(restart == 1, 1e-12, 1e-15)
        nlminb_rtol  <- ifelse(restart == 1, 1e-12, 1e-15)
        nlminb_xtol  <- ifelse(restart == 1, 1e-12, 1e-15)

        # ---- Step A: BFGS in logit-transformed unconstrained space ─────
        CurrentStage <<- paste0("Restart ", restart, " \u2013 BFGS")
        cat("  Step A: BFGS (logit-transformed)\n")
        FnCallNo <<- 0

        bfgs_start     <- model$env$last.par.best
        fn_check       <- model$fn(bfgs_start)
        pre_bfgs_grad  <- max(abs(model$gr_Orig(bfgs_start)))
        if (!is.finite(fn_check)) {
          cat("  WARNING: last.par.best non-finite, falling back to mout$par\n")
          bfgs_start    <- mout$par
          pre_bfgs_grad <- max(abs(model$gr_Orig(bfgs_start)))
        }

        if (has_bounds) {
          lo <- RunSpecs$lowBnd
          hi <- RunSpecs$uppBnd

          bfgs_start_u <- .to_unbounded(bfgs_start, lo, hi)

          fn_u <- function(y) model$fn(.to_bounded(y, lo, hi))
          gr_u <- function(y) {
            x <- .to_bounded(y, lo, hi)
            .chain_grad(model$gr(x), y, lo, hi)
          }

          fit_bfgs <- optim(bfgs_start_u, fn_u, gr_u,
                            method  = "BFGS",
                            control = list(maxit  = MaXeVaL,
                                           reltol = bfgs_reltol))

          fit_bfgs$par <- .to_bounded(fit_bfgs$par, lo, hi)

        } else {
          fit_bfgs <- optim(bfgs_start, model$fn, model$gr,
                            method  = "BFGS",
                            control = list(maxit  = MaXeVaL %/% 2,
                                           reltol = bfgs_reltol))
        }

        bfgs_grad <- max(abs(model$gr_Orig(fit_bfgs$par)))
        cat("  BFGS complete: obj =", round(fit_bfgs$value, 6),
            "| max|grad| =", round(bfgs_grad, 6),
            "| convergence:", fit_bfgs$convergence, "\n")
        TotalEval <<- TotalEval + FnCallNo

        # ---- BFGS gradient revert guard ─────────────────────────────────
        # If BFGS degraded the gradient (common in tight basins on restart
        # >= 2 where the cold-start Hessian approximation overshoots), revert
        # to the pre-BFGS parameters so nlminb starts from the better point.
        if (bfgs_grad > pre_bfgs_grad) {
          cat("  BFGS degraded gradient (", round(pre_bfgs_grad, 4), "->",
              round(bfgs_grad, 4), ") — reverting to pre-BFGS parameters\n",
              sep = "")
          fit_bfgs$par   <- bfgs_start
          fit_bfgs$value <- model$fn(bfgs_start)
        }

        # ---- Step B: nlminb from BFGS solution ──────────────────────────
        CurrentStage <<- paste0("Restart ", restart, " \u2013 nlminb")
        cat("  Step B: nlminb\n")
        FnCallNo    <<- 0
        BestFn      <- fit_bfgs$value
        initBestFn  <- BestFn
        LastPrintFn <<- BestFn
        .trace_append(TotalEval, BestFn)

        ctrl <- list(iter.max = MaXeVaL, eval.max = MaXeVaL,
                     rel.tol  = nlminb_rtol,
                     x.tol    = nlminb_xtol,
                     abs.tol  = 0)
        mout <- nlminb(fit_bfgs$par, model$fn, model$gr,
                       lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                       control = ctrl)

        .report_fit(mout, model, pnames, initBestFn,
                    label = paste("  nlminb restart", restart))
        TotalEval <<- TotalEval + FnCallNo

        cur_grad <- max(abs(model$gr_Orig(mout$par)))
        cat("  Restart", restart, "complete: max|grad| =",
            round(cur_grad, 6), "\n")

        # ---- Exit criteria ───────────────────────────────────────────────
        # 1. Gradient low enough for Newton — clean exit
        if (cur_grad < newton_grad_thresh) {
          cat("  Gradient <", newton_grad_thresh,
              "\u2014 sandwich converged, proceeding to Newton.\n")
          sandwich_done <- TRUE
          break
        }

        # 2. Gradient stall — two consecutive cycles with no meaningful
        #    gradient reduction. NLL may be flat (basin floor) while the
        #    gradient is still improving, so we track gradient not NLL.
        grad_improvement <- prev_grad - cur_grad
        if (grad_improvement < newton_grad_thresh * 0.01) {
          no_improve <- no_improve + 1L
          cat("  Gradient not improving (delta =",
              round(grad_improvement, 6), ") — strike", no_improve, "of 2\n")
          if (no_improve >= 2L) {
            cat("  Sandwich stalled on gradient — exiting after",
                restart, "cycles.\n")
            break
          }
        } else {
          no_improve <- 0L    # reset strike counter on meaningful improvement
        }
        prev_grad <- cur_grad
      }

      # ── Newton polishing ───────────────────────────────────────────────
      suppress_print <<- FALSE       # restore printing for Newton and beyond
      if (newtonSteps > 0) {

        cur_grad <- max(abs(model$gr_Orig(mout$par)))

        if (cur_grad >= newton_grad_thresh) {
          cat("\n  WARNING: gradient =", round(cur_grad, 4),
              ">= newton_grad_thresh (", newton_grad_thresh, ")",
              "— skipping Newton (Hessian likely ill-conditioned).\n")
        } else {
          cat("\n--- Newton polishing steps ---\n")
          CurrentStage <<- "Newton polish"
          newton_par   <- model$env$last.par.best

          for (ns in seq_len(newtonSteps)) {
            improved <- FALSE
            tryCatch({
              H    <- optimHess(newton_par, model$fn, model$gr)
              g    <- as.vector(model$gr_Orig(newton_par))
              step <- solve(H, g)

              # Backtracking line search — halve step until objective improves
              base_obj   <- model$fn(newton_par)
              step_size  <- 1.0
              for (backstep in seq_len(10)) {
                candidate <- newton_par - step_size * step
                if (has_bounds) {
                  candidate <- pmax(candidate, RunSpecs$lowBnd)
                  candidate <- pmin(candidate, RunSpecs$uppBnd)
                }
                if (model$fn(candidate) < base_obj) {
                  newton_par <- candidate
                  improved   <- TRUE
                  break
                }
                step_size <- step_size / 2
              }

              if (!improved) {
                cat("  Newton step", ns,
                    "— line search failed, stopping Newton.\n")
              } else {
                ng <- max(abs(model$gr_Orig(newton_par)))
                cat("  Newton step", ns,
                    "- obj:", round(model$fn(newton_par), 6),
                    "| max|grad|:", round(ng, 8),
                    "| step size:", round(step_size, 6), "\n")
                if (ng < newton_grad_thresh / 100) {
                  cat("  Gradient converged after Newton \u2014 stopping.\n")
                  break
                }
              }
            }, error = function(e) {
              cat("  Newton step", ns, "failed (Hessian singular?):",
                  conditionMessage(e), "\n")
            })
            if (!improved) break
          }

          # Final nlminb from Newton-polished parameters
          CurrentStage <<- "Post-Newton nlminb"
          cat("  Final nlminb (post-Newton)\n")
          FnCallNo    <<- 0
          BestFn      <- model$fn(newton_par)
          initBestFn  <- BestFn
          LastPrintFn <<- BestFn
          .trace_append(TotalEval, BestFn)

          ctrl <- list(iter.max = MaXeVaL, eval.max = MaXeVaL,
                       rel.tol = 1e-12, x.tol = 1e-12, abs.tol = 0)
          mout <- nlminb(newton_par, model$fn, model$gr,
                         lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                         control = ctrl)
          .report_fit(mout, model, pnames, initBestFn,
                      label = "  Final nlminb (post-Newton)")
          TotalEval <<- TotalEval + FnCallNo
        }
      }

      # ── Bound diagnostics ───────────────────────────────────────────────
      .check_bounds(mout$par, RunSpecs$lowBnd, RunSpecs$uppBnd, pnames, model)
    }

    # ── Store and save parameters ─────────────────────────────────────────
    pars        <- mout$par
    names(pars) <- pnames
    ParOld      <- mout$par

    pout <- unlist(parameters)
    pout[names(pout) %in% names(pars)] <- pars
    suffix <- ifelse(is_final, " final", CurrPhase)
    write.table(pout, paste0("Output/model", suffix, ".par"),
                sep = "\t", col.names = c("name\test"), quote = FALSE)

    # ── Report (final phase only) ─────────────────────────────────────────
    if (report && is_final) {
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
        Gradient   = abs(model$gr_Orig(best))
      )
      save(BigSave, file = "Output/BigSave.lda")

      print("making Output.RL")
      WriteOutput(Report, SDrep, fullrep, parameters, pout,
                  GeneralSpecs, ControlSpecs, TheData,
                  CurrPhase = 0, best = best,
                  grad = abs(model$gr_Orig(best)))
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
      "| Max grad [par]:", round(max(Grad), 6), badpar,
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
    cat("  If these have large gradients, use a prior, widen the bound or fix via map.\n\n")
  } else {
    cat("\n  No parameters at bounds :) \n\n")
  }
}

