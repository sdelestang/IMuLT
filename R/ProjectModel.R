#' Project a Fitted IMuLT Model Forward
#'
#' Loads a completed model fit and re-evaluates it with the projection years
#' switched on (\code{Data$DoProject <- 1}), returning a report extended
#' beyond the assessment period according to the specifications parsed from
#' PROJECTIONS.DAT.
#'
#' @param bigsave_file Character. Name/relative path of the saved fit to
#'   project from, located via \code{\link{find_model_file}} (searched for in
#'   the working directory and, if not found there, up to \code{up} parent
#'   directories). Default \code{"Output/BigSave.lda"}, i.e. the file written
#'   by \code{\link{FitModel}(..., report = TRUE)}.
#' @param up Integer. Passed to \code{find_model_file()} — how many parent
#'   directories to search if \code{bigsave_file} isn't in the working
#'   directory. Default 3.
#' @param save Logical. If \code{TRUE} (default), writes the projection
#'   results to \code{"Output/BigSaveProj.lda"} (in the same \code{Output/}
#'   folder as \code{bigsave_file}).
#'
#' @return Invisibly, a list with elements:
#' \itemize{
#'   \item \code{Report} - the model report, extended into the projection years
#'   \item \code{ReportOrig} - the report from the original (unprojected) fit,
#'     for comparison
#'   \item \code{Data} - the \code{Data} list used to build the projection
#'     model (\code{DoProject} set to 1)
#'   \item \code{model} - the \code{MakeADFun} object for the projection run
#' }
#'
#' @details
#' \code{Data$Nproj}, \code{Data$ProjType}, \code{Data$ProjHarvestRate},
#' \code{Data$Catch} (extended into the projection years),
#' \code{Data$SelPntFut}, \code{Data$RetPntFut}, \code{Data$LegalFleetPntFut}
#' and \code{Data$Phi} are already parsed from PROJECTIONS.DAT by
#' \code{ReadProjFile()} inside \code{LoadData()}, and travel with
#' \code{Data} inside \code{BigSave.lda} — so nothing needs to be re-read
#' here. This function only needs to flip \code{Data$DoProject} to 1 and
#' rebuild the AD model so the C++ side runs its projection loop instead of
#' stopping at \code{Nyear}. \code{Data$ProjType} (1 = catch-based,
#' 2 = harvest-rate-based) controls which schedule the compiled model
#' actually uses for the projection years.
#'
#' No re-optimisation happens here: the model is evaluated once at the saved
#' MLE (\code{BigSave$best}), so the objective value should match the
#' original fit (the projection years add no likelihood contribution) and the
#' projected years reflect a straightforward forward simulation, not a
#' re-estimated one.
#'
#' @examples
#' \dontrun{
#' choose_model()
#' proj <- ProjectModel()
#' print(proj$Report$MatBio)
#' print(proj$ReportOrig$MatBio)
#' }
#'
#' @seealso \code{\link{FitModel}}, \code{\link{choose_model}},
#'   \code{\link{find_model_file}}
#'
#' @export
ProjectModel <- function(bigsave_file = "Output/BigSave.lda", up = 3L, save = TRUE) {

  # ── Locate and load the completed fit ─────────────────────────────────────
  bigsave_path <- find_model_file(bigsave_file, up = up)
  cat("Loading", bigsave_path, "\n")

  ee <- new.env()
  load(file = bigsave_path, envir = ee)
  if (!exists("BigSave", envir = ee, inherits = FALSE))
    stop("'", bigsave_path, "' did not contain an object called 'BigSave'.",
         call. = FALSE)
  BigSave <- ee$BigSave

  need <- c("Data", "map", "best", "parameters", "Report")
  miss <- need[!need %in% names(BigSave)]
  if (length(miss))
    stop("ProjectModel: BigSave is missing ", paste(miss, collapse = ", "),
         " \u2014 was it saved by FitModel(..., report = TRUE)?", call. = FALSE)

  Data       <- BigSave$Data
  map        <- BigSave$map
  bestvals   <- BigSave$best
  parameters <- BigSave$parameters
  ReportOrig <- BigSave$Report

  # MatBio (and the other REPORT()'d arrays) are always sized to
  # BurnIn+Nyear+MaxProjYr+1 -- MaxProjYr is the array's maximum capacity, not
  # how many years this run actually filled in. Index the real last computed
  # year explicitly rather than tail(x,1), or you'll read an unfilled zero
  # slot whenever MaxProjYr > the number of years actually run.
  last_hist <- Data$BurnIn + Data$Nyear
  cat("Original fit MatBio, final assessment year:",
      round(ReportOrig$MatBio[last_hist], 3), "\n")

  if (is.null(Data$Nproj) || Data$Nproj == 0)
    warning("Data$Nproj is 0 \u2014 PROJECTIONS.DAT specified no projection ",
            "years, so ProjectModel() will run but nothing beyond the ",
            "assessment period will actually be projected.", call. = FALSE)

  # ── Switch on projections and rebuild the AD model ─────────────────────────
  Data$DoProject <- 1

  cat("Making projection model object (DoProject = 1, Nproj =",
      Data$Nproj, ")\n")
  model <- MakeADFun(Data, parameters, map = map, DLL = "IMuLT", silent = TRUE)

  model$par <- bestvals
  nll <- model$fn(model$par)
  cat("Objective at saved parameter values:", round(nll, 6), "\n")

  cat("Extracting projection report\n")
  Report <- model$report()

  last_proj <- Data$BurnIn + Data$Nyear + Data$Nproj
  cat("Projected MatBio, final projection year:",
      round(Report$MatBio[last_proj], 3), "\n")

  # ── Save ─────────────────────────────────────────────────────────────────
  if (save) {
    outdir <- dirname(bigsave_path)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    ProjSave <- list(Data       = Data,
                     Report     = Report,
                     ReportOrig = ReportOrig,
                     best       = bestvals,
                     parameters = parameters,
                     map        = map)
    outpath <- file.path(outdir, "BigSaveProj.lda")
    save(ProjSave, file = outpath)
    cat("Saved projection to", outpath, "\n")
  }

  invisible(list(Report = Report, ReportOrig = ReportOrig,
                 Data = Data, model = model))
}
