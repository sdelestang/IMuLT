#' Project a Fitted IMuLT Model Forward
#'
#' Loads a completed model fit and re-evaluates it with the projection years
#' switched on (\code{Data$DoProject <- 1}), then (by default) builds an HTML
#' report from the result covering the projection specification, the
#' selectivity/retention/legal-size patterns actually used in the projection
#' years, the input catch or harvest-rate schedule, discards, fishing
#' efficiency, legal biomass, and mature (egg-production proxy) biomass.
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
#' @param report Logical. If \code{TRUE} (default), builds the HTML
#'   projection report described above straight after the projection run,
#'   using the \code{Report}/\code{Data} objects already in memory (no reload
#'   from disk).
#' @param is95 Logical. Use 95\% CIs (\code{TRUE}, default) or 68\% (\code{FALSE})
#'   where relevant in the report.
#' @param folder_name Name for the report folder inside \code{Summary/}. A
#'   \code{"Proj"} suffix is appended if not already present. Blank (default)
#'   prompts interactively, defaulting to \code{"resultProj"}. Ignored if
#'   \code{report = FALSE}.
#' @param openfile Whether to open the HTML report on completion. Default
#'   TRUE. Ignored if \code{report = FALSE}.
#' @param ProjType Integer, \code{1} or \code{2}, or \code{NULL} (default).
#'   If \code{NULL}, uses whatever \code{ProjType} PROJECTIONS.DAT specified
#'   (parsed by \code{ReadProjFile()} at \code{LoadData()} time, travelling
#'   with \code{Data} inside \code{BigSave.lda}). If set, overrides
#'   \code{Data$ProjType} for this run only -- e.g. \code{ProjectModel(ProjType = 2)}
#'   runs a harvest-rate-based projection even if PROJECTIONS.DAT specified
#'   catch-based (1), without needing to edit or regenerate the file. Since
#'   both the catch and harvest-rate schedules are always parsed into
#'   \code{Data} regardless of the file's ProjType flag, overriding here is
#'   safe as long as the schedule you actually want was populated in
#'   PROJECTIONS.DAT.
#' @param sdreport Logical. If \code{TRUE}, calls \code{TMB::sdreport(model)}
#'   after the projection report is extracted, computing delta-method
#'   standard errors for the ADREPORT()'d projection quantities (\code{MatBio},
#'   \code{MatBioArea}, \code{RecruitmentByArea}, \code{LegalBioAll},
#'   \code{HarvestRate}, \code{PredCpue}, \code{CpueEcreep}). Default
#'   \code{FALSE}, since this can be slow (a full Hessian at the projection's
#'   parameter dimension). The result is stored as \code{SDrep} in both the
#'   returned list and \code{BigSaveProj.lda}.
#'
#' @return Invisibly, a list with elements:
#' \itemize{
#'   \item \code{Report} - the model report, extended into the projection years
#'   \item \code{ReportOrig} - the report from the original (unprojected) fit,
#'     for comparison
#'   \item \code{Data} - the \code{Data} list used to build the projection
#'     model (\code{DoProject} set to 1)
#'   \item \code{model} - the \code{MakeADFun} object for the projection run
#'   \item \code{SDrep} - result of \code{TMB::sdreport(model)} if
#'     \code{sdreport = TRUE}, else \code{NULL}
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
#' The report step works directly off the in-memory \code{Report}/\code{Data}
#' objects rather than text output files -- there is no projection equivalent
#' of \code{WriteOutput()} yet. Only fields confirmed to actually extend into
#' the projection years are used: \code{LegalBioAll}, \code{LegalBioAllbySex},
#' \code{MatureBioAllbySex}, \code{MatBio}, \code{MatBioArea},
#' \code{DiscardWt}, \code{DeadDiscardWt}, \code{Hrate}, \code{CpueEcreep},
#' \code{VirginLegalBio}, \code{ActSelex}, \code{ActReten}, \code{ActLegal}.
#' An "Index" (CPUE fit) tab is deliberately omitted, since there is no
#' future CPUE data to fit against.
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
#'   \code{\link{find_model_file}}, \code{\link{MakeOutPut}}
#'
#' @export
ProjectModel <- function(bigsave_file = "Output/BigSave.lda", up = 3L,
                         save = TRUE, report = TRUE, is95 = TRUE,
                         folder_name = '', openfile = TRUE, ProjType = NULL,
                         sdreport = FALSE) {

  # ── Locate and load the completed fit ─────────────────────────────────────
  bigsave_path <- find_model_file(bigsave_file, up = up)
  cat("Loading", bigsave_path, "\n")

  ee <- new.env()
  load(file = bigsave_path, envir = ee)
  if (!exists("BigSave", envir = ee, inherits = FALSE))
    stop("'", bigsave_path, "' did not contain an object called 'BigSave'.",
         call. = FALSE)
  BigSave <- ee$BigSave

  Data       <- BigSave$Data
  map        <- BigSave$map
  bestvals   <- BigSave$best
  parameters <- BigSave$parameters
  ReportOrig <- BigSave$Report

  # Find projections file and reload
  ProjFilePath <- find_model_file("PROJECTIONS.DAT", up = 5L, highest = TRUE)
  ProjFile <- read.table(ProjFilePath, comment.char = "?", fill = TRUE,
                         blank.lines.skip = TRUE, stringsAsFactors = FALSE,col.names = 1:200)
  ProjectSpecs <- ReadProjFile(ProjFile, Data, Data$Phi1, Data$Catch)
  # Overwrite projection info
  Data$Nproj            <- ProjectSpecs$Nproj
  Data$SelPntFut        <- ProjectSpecs$SelPntFut
  Data$RetPntFut        <- ProjectSpecs$RetPntFut
  Data$LegalFleetPntFut <- ProjectSpecs$LegalFleetPntFut
  Data$Phi              <- ProjectSpecs$Phi
  if (!is.null(ProjType)) {
    cat("Overriding PROJECTIONS.DAT's ProjType (", ProjectSpecs$ProjType,
        ") with ProjType =", ProjType, "\n")
    Data$ProjType <- ProjType
  } else {
    Data$ProjType <- ProjectSpecs$ProjType
  }
  Data$Catch           <- ProjectSpecs$Catch
  Data$ProjHarvestRate <- ProjectSpecs$ProjHarvestRate

  need <- c("Data", "map", "best", "parameters", "Report")
  miss <- need[!need %in% names(BigSave)]
  if (length(miss))
    stop("ProjectModel: BigSave is missing ", paste(miss, collapse = ", "),
         " \u2014 was it saved by FitModel(..., report = TRUE)?", call. = FALSE)

  # LegalBioAll (and the other REPORT()'d arrays) are always sized to
  # BurnIn+Nyear+MaxProjYr -- MaxProjYr is the array's maximum capacity, not
  # how many years this run actually filled in. Index the real last computed
  # year explicitly rather than tail(x,1), or you'll read an unfilled zero
  # slot whenever MaxProjYr > the number of years actually run.
  last_hist <- Data$BurnIn + Data$Nyear
  relbio_hist <- sum(ReportOrig$LegalBioAll[last_hist,]) / sum(ReportOrig$VirginLegalBio)
  cat("Original fit relative biomass (B/B0), final assessment year:",
      round(relbio_hist, 3), "\n")

  if (is.null(Data$Nproj) || Data$Nproj == 0)
    warning("Data$Nproj is 0 \u2014 PROJECTIONS.DAT specified no projection ",
            "years, so ProjectModel() will run but nothing beyond the ",
            "assessment period will actually be projected.", call. = FALSE)

  # ── Switch on projections and rebuild the AD model ─────────────────────────
  Data$DoProject <- 1

  cat("Making projection model object (DoProject = 1, Nproj =",
      Data$Nproj, ", ProjType =", Data$ProjType, ")\n")
  model <- MakeADFun(Data, parameters, map = map, DLL = "IMuLT", silent = TRUE)

  model$par <- bestvals
  nll <- model$fn(model$par)
  cat("Objective at saved parameter values:", round(nll, 6), "\n")

  cat("Extracting projection report\n")
  Report <- model$report()

  SDrep <- NULL
  if (sdreport) {
    cat("Running sdreport() for the projection's error structure ",
        "(this can take a while) ...\n", sep = "")
    SDrep <- TMB::sdreport(model)
    cat("sdreport() complete\n")
  }

  last_proj <- Data$BurnIn + Data$Nyear + Data$Nproj
  relbio_proj <- sum(Report$LegalBioAll[last_proj,]) / sum(Report$VirginLegalBio)
  cat("Projected relative biomass (B/B0), final projection year:",
      round(relbio_proj, 3), "\n")

  # ── Save ─────────────────────────────────────────────────────────────────
  outdir <- dirname(bigsave_path)
  if (save) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    file.copy(ProjFilePath, file.path(outdir, "PROJECTIONS.DAT"), overwrite = TRUE)
    ProjSave <- list(Data       = Data,
                     Report     = Report,
                     ReportOrig = ReportOrig,
                     best       = bestvals,
                     parameters = parameters,
                     map        = map,
                     SDrep      = SDrep)
    outpath <- file.path(outdir, "BigSaveProj.lda")
    save(ProjSave, file = outpath)
    cat("Saved projection to", outpath, "\n")
  }

  # ── Report ───────────────────────────────────────────────────────────────
  if (report) {
    .MakeProjectionReport(Data = Data, Report = Report, ReportOrig = ReportOrig,
                          rundir_base = outdir, up = up, is95 = is95,
                          folder_name = folder_name, openfile = openfile,
                          SDrep = SDrep)
  }

  invisible(list(Report = Report, ReportOrig = ReportOrig,
                 Data = Data, model = model, SDrep = SDrep))
}


#' Build the HTML Projection Report
#'
#' Internal helper called from \code{\link{ProjectModel}} when
#' \code{report = TRUE}. Not exported -- always called with the
#' \code{Data}/\code{Report} objects already in memory from the current
#' projection run, never reloaded from disk.
#'
#' @param Data,Report,ReportOrig The projection run's data/report objects.
#' @param rundir_base The model's Output/ directory (where BigSaveProj.lda
#'   lives) -- the report is written to Output/Projections/<folder_name>/.
#'   DATA.DAT, CONTROL.DAT, and ModelStructure.xlsx are located separately
#'   via find_model_file()/load_model_structure(), which search from the
#'   working directory rather than from this path.
#' @param up Passed to \code{find_model_file()} for locating \code{DATA.DAT}
#'   and \code{ModelStructure.xlsx}.
#' @param is95,folder_name,openfile See \code{\link{ProjectModel}}.
#' @param SDrep Result of \code{TMB::sdreport(model)}, or \code{NULL}
#'   (default). When provided, CI ribbons (delta-method SEs) are added to
#'   the Fishing Efficiency, Biomass, and Harvest Rate plots -- the three
#'   the user asked for, out of everything currently \code{ADREPORT()}'d.
#'
#' @return Invisibly NULL. Called for its side effect of writing the report.
#' @keywords internal
.MakeProjectionReport <- function(Data, Report, ReportOrig, rundir_base, up,
                                  is95, folder_name, openfile, SDrep = NULL) {

  suppressPackageStartupMessages({
    library(makehtml); library(hplot); library(dplyr); library(magrittr)
    library(tidyr); library(ggplot2); library(reshape2); library(openxlsx)
  })
  options(dplyr.summarise.inform = FALSE)
  starttime <- as.character(Sys.time())
  SclErr <- ifelse(is95, 1.96, 0.842)

  # Reshaped SEs (same array/matrix shape as the corresponding Report
  # element) for whichever quantities were actually ADREPORT()'d -- NULL
  # entries (e.g. if SDrep is NULL, or a quantity was commented out of
  # ADREPORT in the .cpp) are handled per-plot below by falling back to a
  # plain line with no ribbon.
  SDse <- if (!is.null(SDrep)) as.list(SDrep, what = "Std. Error", report = TRUE) else NULL

  Nyear   <- Data$Nyear
  Nproj   <- Data$Nproj
  BurnIn  <- Data$BurnIn
  Year1   <- Data$Year1
  Nstep   <- Data$Nstep
  Nfleet  <- Data$Nfleet
  Narea   <- Data$Narea
  Nsex    <- Data$Nsex
  ProjType <- Data$ProjType
  projyears <- (Year1+Nyear):(Year1+Nyear+Nproj-1)

  # Row-index helpers -- two conventions are in play: BurnIn-offset arrays
  # (LegalBioAll, MatBio, Hrate, ...) vs Year1-origin arrays (DiscardWt,
  # CpueEcreep, Data$Catch, ...).
  .row_burnin <- function(year) BurnIn + (year - Year1) + 1
  .row_year1  <- function(year) (year - Year1) + 1
  # Years are integers; ggplot's default pretty() breaks can still pick
  # fractional positions (e.g. 2027.5) on a narrow range -- force integers.
  .int_breaks <- function(x) unique(round(pretty(x)))

  # ── Output directory ─────────────────────────────────────────────────────
  indir <- filenametopath(rundir_base, "Projections")
  if (folder_name == '') {
    run_name <- dlg_input("Name for this projection report (blank = 'resultProj'):")$res
    if (run_name == '') run_name <- 'resultProj'
  } else run_name <- folder_name
  if (!grepl("Proj$", run_name)) run_name <- paste0(run_name, "Proj")
  rundir <- filenametopath(indir, run_name)

  if (dir.exists(rundir)) {
    f <- list.files(rundir, include.dirs = FALSE, full.names = TRUE, recursive = TRUE)
    if (length(f) > 0) {
      Archive <- toupper(dlg_input(paste0("'", run_name, "' already exists. Archive it first? (Y or N)"))$res)
      if (Archive == 'Y') {
        ctime <- gsub(' ', '', gsub(".", "", format(Sys.time(), '%Y.%m.%d %H.%M'), fixed = TRUE))
        rundirA <- paste0(indir, '/archive', ctime)
        invisible(dir.create(rundirA))
        invisible(file.copy(rundir, rundirA, recursive = TRUE))
        print("Archived old projection report")
      }
      suppressWarnings(invisible(file.remove(f)))
    }
  }
  dirExists(rundir, verbose = TRUE)
  resfile <- setuphtml(rundir = rundir)

  # ── Reference data: fleet/area names, length bins ───────────────────────
  wb <- load_model_structure()
  fleets <- readWorkbook(wb, sheet = 'Fleetcode', startRow = 2)
  areas  <- readWorkbook(wb, sheet = 'Area', startRow = 2)
  suppressWarnings(fleetarea <- fleets %>%
    mutate(areaname = do.call('rbind', strsplit(description, '_'))[,1], fleettype = group))

  datadat_path <- find_model_file("DATA.DAT", up = up)
  lbin1 <- read.table(datadat_path, comment.char = "?", fill = TRUE,
                      blank.lines.skip = TRUE, stringsAsFactors = FALSE, col.names = 1:200)
  lbin <- findNclean(c('#', 'Lower'), lbin1, 1, convert = 0)
  lbin <- as.numeric(unname(lbin[1,]))
  lbin <- lbin + diff(lbin)[1] / 2
  lbin <- lbin[1:(length(lbin) - 1)]

  ctl_path <- find_model_file("CONTROL.DAT", up = up)
  ctl1 <- read.table(ctl_path, comment.char = "?", fill = TRUE,
                     blank.lines.skip = TRUE, stringsAsFactors = FALSE, col.names = 1:200)
  reflev <- findNclean(c('#', 'Biomass', 'target'), ctl1, 1)
  colnames(reflev) <- c('target', 'threshold', 'limit')

  #### Run summary (native Home page notes, not a separate tab) ####
  homelines <- c(
    paste("Projection years:", Nproj),
    paste("First projection year:", min(projyears)),
    paste("Last projection year:", max(projyears)),
    paste("Basis:", ifelse(ProjType == 1, "Catch-based (ProjType 1)", "Harvest-rate-based (ProjType 2)")),
    paste("Number of fleets:", Nfleet),
    paste("Number of areas:", Narea))

  #### Selectivity / High-grading / Legal / Combined used in projections ####
  print("Making Selectivity/High-grading/Legal-size used in projections")
  .plot_proj_curve <- function(pairs, label, category, caption) {
    # pairs: named list of list(Pnt=<PntFut array>, Act=<pattern table>).
    # Pnt: array (Nsex,Nage,Nfleet,MaxProjYr,Nstep); pulls sex=1,age=1 (age
    # doesn't vary pattern selection in this model) for each fleet/proj-year/
    # step, resolves against Act (0-indexed pointers). When more than one
    # pair is given, the resolved curves are multiplied together elementwise
    # (e.g. selectivity x high-grading x legal-size = overall proportion
    # actually retained).
    rows <- list()
    for (f in 1:Nfleet) {
      for (ys in seq_along(projyears)) {
        for (st in 1:Nstep) {
          curves <- lapply(pairs, function(p) {
            ptr <- p$Pnt[1,1,f,ys,st]
            if (is.na(ptr) || ptr < 0) return(NULL)
            list(ptr = ptr, curve = as.numeric(p$Act[ptr + 1, ])[seq_along(lbin)])
          })
          if (any(vapply(curves, is.null, logical(1)))) next
          combined <- Reduce(`*`, lapply(curves, `[[`, "curve"))
          ptrid <- paste(vapply(curves, function(x) x$ptr, numeric(1)), collapse = ".")
          rows[[length(rows)+1]] <- data.frame(
            fleet = f, year = projyears[ys], step = st, ptr = ptrid,
            lbin = lbin, value = combined)
        }
      }
    }
    if (length(rows) == 0) return(invisible(NULL))
    df <- bind_rows(rows) %>%
      mutate(descrip = fleetarea$description[match(fleet, fleetarea$fleet)],
             curve_id = paste(fleet, ptr, step))
    yr_lab <- df %>% group_by(curve_id, fleet, descrip, ptr) %>%
      summarise(yr1 = min(year), yr2 = max(year), .groups = "drop") %>%
      mutate(lab = ifelse(yr1==yr2, paste0("Pattern ", ptr, " (", yr1, ")"),
                          paste0("Pattern ", ptr, " (", yr1, "-", yr2, ")")))
    df <- left_join(df, yr_lab %>% dplyr::select(curve_id, lab), by = "curve_id")

    fname <- gsub("[^A-Za-z0-9]", "_", label)
    filename <- filenametopath(rundir, paste0(fname, "_projected.png"))
    plotprep(width = 10, height = 10, filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1,1))
    print(ggplot(df, aes(x = lbin, y = value, colour = lab)) +
            geom_line(linewidth = 0.8) +
            facet_wrap(~descrip) +
            ylim(0, 1) +
            labs(x = "Length bin (mm)", y = label, colour = NULL) +
            theme_bw() + theme(legend.position = "bottom"))
    addplot(filen = filename, rundir = rundir, category = category, caption = caption)
  }
  selpair <- list(Pnt = Data$SelPntFut, Act = Report$ActSelex)
  hgpair  <- list(Pnt = Data$RetPntFut, Act = Report$ActReten)
  legpair <- list(Pnt = Data$LegalFleetPntFut, Act = Report$ActLegal)
  .plot_proj_curve(list(sel = selpair), "Selectivity", "Selectivity_Retention",
                   "Selectivity of the gear applied during the projection years, by fleet.")
  .plot_proj_curve(list(hg = hgpair), "High-grading", "Selectivity_Retention",
                   "High-grading pattern applied during the projection years, by fleet.")
  .plot_proj_curve(list(leg = legpair), "Legal size", "Selectivity_Retention",
                   "Legal-size assignment applied during the projection years, by fleet.")
  .plot_proj_curve(list(sel = selpair, hg = hgpair, leg = legpair),
                   "Overall retained", "Selectivity_Retention",
                   paste("Combined proportion of encountered animals actually retained",
                         "(selectivity x high-grading x legal size), by fleet."))

  #### Catches and Harvest Rate schedules -- both always shown regardless of
  #### ProjType (PROJECTIONS.DAT carries both columns in every row), so you
  #### can see what the non-driving one does too ####
  print("Making projected Catch schedule")
  rows <- list()
  for (f in 1:Nfleet) for (yr in projyears) for (st in 1:Nstep) {
    v <- Data$Catch[.row_year1(yr), st, f]
    if (!is.na(v)) rows[[length(rows)+1]] <- data.frame(fleet=f, year=yr, step=st, catch=v)
  }
  catchdf <- bind_rows(rows) %>%
    mutate(areaname = fleetarea$areaname[match(fleet, fleetarea$fleet)],
           descrip  = fleetarea$description[match(fleet, fleetarea$fleet)]) %>%
    filter(catch > 0)
  if (nrow(catchdf) > 0) {
    filename <- filenametopath(rundir, "Projected_Catch.png")
    plotprep(width = 10, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1,1))
    print(ggplot(catchdf %>% group_by(year, areaname) %>% summarise(CatchT = sum(catch)/1000, .groups="drop"),
                 aes(x = year, y = CatchT, fill = areaname)) +
            geom_bar(stat = "identity", position = "stack") +
            viridis::scale_fill_viridis(discrete = TRUE) +
            scale_x_continuous(breaks = .int_breaks) +
            ylab("Catch (t)") + xlab("Year") +
            theme(panel.background = element_rect(fill = "white", colour = NA),
                  panel.border = element_rect(fill = NA, colour = "grey20")))
    catch_caption <- if (ProjType == 1)
      "Input catch schedule for the projection years, by area (ProjType 1 -- this is what drives the projection)."
    else
      paste("Reference catch schedule for the projection years, by area (ProjType 2 -- Harvest",
            "Rate is what actually drives this projection; these numbers are scaled from the last",
            "historical year's catch proportions to the target projected catch level, shown here",
            "for context, not as a realised/fitted value).")
    addplot(filen = filename, rundir = rundir, category = "Catches", caption = catch_caption)

    catchtab <- catchdf %>% group_by(year, descrip) %>% summarise(catch = sum(catch), .groups = "drop") %>% mutate(catch=round(catch/1000,1)) %>%
      pivot_wider(names_from = year, values_from = catch, values_fill = 0)
    catchtab_caption <- if (ProjType == 1)
      "Catch (t) applied by fleet and projection year (ProjType 1 -- this is the input schedule, not a fitted/realised value)."
    else
      "Reference catch (t) by fleet and projection year (ProjType 2 -- context only, not what the imposed harvest rate actually produces)."
    addtable(intable = catchtab, filen = "Projected_Catch.csv", rundir = rundir,
             category = "Catches", caption = catchtab_caption)
  }

  print("Making projected Harvest Rate schedule")
  rows <- list()
  for (f in 1:Nfleet) for (ys in seq_along(projyears)) for (st in 1:Nstep) {
    v <- Data$ProjHarvestRate[ys, st, f]
    if (!is.na(v) && v > 0) rows[[length(rows)+1]] <-
      data.frame(fleet=f, year=projyears[ys], step=st, hrate=v)
  }
  hrdf <- bind_rows(rows) %>%
    mutate(descrip = fleetarea$description[match(fleet, fleetarea$fleet)])
  if (nrow(hrdf) > 0) {
    filename <- filenametopath(rundir, "Projected_HarvestRate_input.png")
    plotprep(width = 10, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1,1))
    print(ggplot(hrdf, aes(x = year, y = hrate, colour = descrip)) +
            geom_line(linewidth = 0.8) + geom_point() +
            scale_x_continuous(breaks = .int_breaks) +
            ylab("Input harvest rate") + xlab("Year") + theme_bw())
    hrate_caption <- if (ProjType == 2)
      "Input harvest-rate schedule for the projection years, by fleet (ProjType 2 -- this is what drives the projection)."
    else
      "Reference harvest-rate schedule for the projection years, by fleet (ProjType 1 -- Catch is what actually drives this projection; shown here for context)."
    addplot(filen = filename, rundir = rundir, category = "Catches", caption = hrate_caption)

    hrtab <- hrdf %>% group_by(year, descrip) %>% summarise(hrate = mean(hrate), .groups = "drop") %>%
      pivot_wider(names_from = year, values_from = hrate, values_fill = 0)
    hrtab_caption <- if (ProjType == 2)
      "Input harvest rate by fleet and projection year (ProjType 2 -- this is what drives the projection)."
    else
      "Reference harvest rate by fleet and projection year (ProjType 1 -- context only, not what actually drives this projection)."
    addtable(intable = hrtab, filen = "Projected_HarvestRate_input.csv", rundir = rundir,
             category = "Catches", caption = hrtab_caption)
  }

  #### Discards ####
  print("Making Discards")
  disc_rows <- list()
  for (f in 1:Nfleet) for (yr in c(projyears[1]-1, projyears)) for (st in 1:Nstep) {
    ri <- .row_year1(yr)
    if (ri < 1 || ri > dim(Report$DiscardWt)[1]) next
    disc_rows[[length(disc_rows)+1]] <- data.frame(
      year = yr, step = st, fleet = f,
      DiscardWt = Report$DiscardWt[ri, st, f],
      DeadDiscardWt = Report$DeadDiscardWt[ri, st, f])
  }
  discdf <- bind_rows(disc_rows) %>%
    mutate(AreaName = fleetarea$areaname[match(fleet, fleetarea$fleet)]) %>%
    group_by(year, AreaName) %>%
    summarise(Discard = sum(DiscardWt)/1000, DeadDiscard = sum(DeadDiscardWt)/1000, .groups = "drop")
  filename <- filenametopath(rundir, "Projected_Discards.png")
  plotprep(width = 10, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  print(ggplot(discdf %>% pivot_longer(c(Discard,DeadDiscard), names_to="Type", values_to="Weight_t"),
               aes(x = year, y = Weight_t, colour = Type)) +
          geom_line() + geom_point() + facet_wrap(~AreaName) +
          scale_color_manual(values = c("red","black"), labels = c("Dead Discards","Total Discards")) +
          scale_x_continuous(breaks = .int_breaks) +
          ylab("Discard Weight (t)") + xlab("Year") +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(angle = 45)))
  caption <- "Projected discards (black) and dead discards (red) by area. Last assessment year included for continuity."
  addplot(filen = filename, rundir = rundir, category = "Discards", caption = caption)

  #### Fishing Efficiency ####
  print("Making Fishing Efficiency")
  eff_rows <- list()
  Nefseries <- ncol(Report$CpueEcreep)
  for (ef in 1:Nefseries) for (yr in c(projyears[1]-1, projyears)) {
    ri <- .row_year1(yr)
    if (ri < 1 || ri > dim(Report$CpueEcreep)[1]) next
    se <- if (!is.null(SDse$CpueEcreep)) SDse$CpueEcreep[ri, ef] else NA_real_
    eff_rows[[length(eff_rows)+1]] <- data.frame(series = ef, year = yr, value = Report$CpueEcreep[ri, ef], se = se)
  }
  effdf <- bind_rows(eff_rows) %>% mutate(lwr = value - se*SclErr, upr = value + se*SclErr)
  filename <- filenametopath(rundir, "Projected_Fishing_Efficiency.png")
  plotprep(width = 8, height = 6, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  p <- ggplot(effdf, aes(x = year, y = value, colour = factor(series)))
  if (!is.null(SDse$CpueEcreep))
    p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr, fill = factor(series)),
                         colour = NA, alpha = 0.2, show.legend = FALSE)
  p <- p + geom_line(linewidth = 0.8) + geom_point() +
    scale_x_continuous(breaks = .int_breaks) +
    labs(x = "Year", y = "Fishing efficiency", colour = "Series") + theme_bw()
  print(p)
  caption <- "Fishing efficiency held at its last estimated value through the projection years."
  addplot(filen = filename, rundir = rundir, category = "Fishing_Efficiency", caption = caption)

  #### Biomass ####
  print("Making Biomass")
  bio_rows <- list()
  for (a in 1:Narea) for (yr in c(projyears[1]-1, projyears)) {
    ri <- .row_burnin(yr)
    se <- if (!is.null(SDse$LegalBioAll)) SDse$LegalBioAll[ri, a] else NA_real_
    bio_rows[[length(bio_rows)+1]] <- data.frame(
      area = a, year = yr, est = Report$LegalBioAll[ri, a],
      virgin = Report$VirginLegalBio[a], se = se)
  }
  biodf <- bind_rows(bio_rows) %>%
    mutate(areaname = areas$Name[match(area, areas$AreaCode)],
           `Legal Biomass (t)` = est/1000, `B/B0` = est/virgin,
           lwr = pmax((est - se*SclErr)/virgin, 0), upr = (est + se*SclErr)/virgin)
  filename <- filenametopath(rundir, "Projected_Legal_Biomass.png")
  plotprep(width = 8, height = 6, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  # Total across areas: SE combined as sqrt(sum(se^2)) (independence approximation,
  # same convention MakeOutPut() already uses for its Egg Production by-BMSA total).
  biosum <- biodf %>% group_by(year) %>%
    summarise(est = sum(est), virgin = sum(virgin), se = sqrt(sum(se^2)), .groups = "drop") %>%
    mutate(`B/B0` = est/virgin, lwr = pmax((est - se*SclErr)/virgin, 0), upr = (est + se*SclErr)/virgin)
  p <- ggplot(biosum, aes(x = year, y = `B/B0`))
  if (!is.null(SDse$LegalBioAll)) p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.5)
  p <- p + geom_line(linewidth = 1) + geom_point() +
    geom_vline(xintercept = projyears[1]-0.5, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = reflev$target, colour = 'green') +
    geom_hline(yintercept = reflev$threshold, colour = 'orange') +
    geom_hline(yintercept = reflev$limit, colour = 'red') +
    scale_x_continuous(breaks = .int_breaks) +
    ylim(0, NA) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20")) +
    ylab("B/B0") + xlab("Year")
  print(p)
  caption <- "Total legal biomass relative to virgin, projected forward (dashed line = last assessment year)."
  addplot(filen = filename, rundir = rundir, category = "Biomass", caption = caption)

  filename <- filenametopath(rundir, "Projected_Legal_Biomass_Area.png")
  plotprep(width = 8, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  p <- ggplot(biodf, aes(x = year, y = `B/B0`))
  if (!is.null(SDse$LegalBioAll)) p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.5)
  p <- p + geom_line() + geom_point(size = 0.75) +
    geom_vline(xintercept = projyears[1]-0.5, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = reflev$target, colour = 'green') +
    geom_hline(yintercept = reflev$threshold, colour = 'orange') +
    geom_hline(yintercept = reflev$limit, colour = 'red') +
    facet_wrap(~areaname) + ylim(0, NA) +
    scale_x_continuous(breaks = .int_breaks) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.text.x = element_text(angle = 45))
  print(p)
  caption <- "Legal biomass relative to virgin by area, projected forward."
  addplot(filen = filename, rundir = rundir, category = "Biomass", caption = caption)

  #### Egg Production (proxy: MatureBioAllbySex) ####
  print("Making Egg Production (Mature Biomass proxy)")
  egg_rows <- list()
  for (a in 1:Narea) for (s in 1:Nsex) for (yr in c(projyears[1]-1, projyears)) {
    ri <- .row_burnin(yr)
    egg_rows[[length(egg_rows)+1]] <- data.frame(
      area = a, sex = s, year = yr, est = Report$MatureBioAllbySex[ri, a, s])
  }
  eggdf <- bind_rows(egg_rows) %>%
    mutate(areaname = areas$Name[match(area, areas$AreaCode)]) %>%
    group_by(year, areaname) %>% summarise(est = sum(est), .groups = "drop") %>%
    group_by(areaname) %>% mutate(rel = est/max(est))
  filename <- filenametopath(rundir, "Projected_Mature_Biomass.png")
  plotprep(width = 8, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  print(ggplot(eggdf, aes(x = year, y = rel)) +
          geom_line() + geom_point(size = 0.75) +
          geom_vline(xintercept = projyears[1]-0.5, linetype = "dashed", colour = "grey50") +
          facet_wrap(~areaname) + ylim(0,1.05) +
          scale_x_continuous(breaks = .int_breaks) +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(angle = 45)) +
          ylab("Relative Mature Biomass") + xlab("Year"))
  caption <- paste("Relative mature biomass (egg production proxy -- MatureBioAllbySex, not the full",
                   "length-based egg-production curve, which does not yet extend into projection years),",
                   "by area, projected forward.")
  addplot(filen = filename, rundir = rundir, category = "Egg_Production", caption = caption)

  #### Harvest Rate ####
  print("Making Harvest Rate")
  hr_rows <- list()
  for (f in 1:Nfleet) for (yr in c(projyears[1]-1, projyears)) for (st in 1:Nstep) {
    ri <- .row_burnin(yr)
    if (ri < 1 || ri > dim(Report$Hrate)[1]) next
    se <- if (!is.null(SDse$Hrate)) SDse$Hrate[ri, st, f] else NA_real_
    hr_rows[[length(hr_rows)+1]] <- data.frame(
      year = yr, step = st, fleet = f, Hrate = Report$Hrate[ri, st, f], se = se)
  }
  hratedf <- bind_rows(hr_rows) %>%
    mutate(descrip = fleetarea$description[match(fleet, fleetarea$fleet)]) %>%
    filter(!is.na(descrip)) %>%
    group_by(fleet) %>% filter(any(Hrate > 0)) %>% ungroup() %>%
    mutate(lwr = pmax(Hrate - se*SclErr, 0), upr = Hrate + se*SclErr)
  filename <- filenametopath(rundir, "Projected_Harvest_Rate.png")
  plotprep(width = 10, height = 10, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  p <- ggplot(hratedf, aes(x = year, y = Hrate))
  if (!is.null(SDse$Hrate)) p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.5)
  p <- p + geom_line() + geom_point(size = 0.75) +
    geom_vline(xintercept = projyears[1]-0.5, linetype = "dashed", colour = "grey50") +
    facet_wrap(~descrip, scales = "free_y") +
    expand_limits(y = 0) +
    scale_x_continuous(breaks = .int_breaks) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          axis.text.x = element_text(angle = 45)) +
    ylab("Harvest rate") + xlab("Year")
  print(p)
  caption <- "Projected harvest rate by fleet (fleets with zero harvest rate throughout are omitted). Last assessment year included for continuity."
  addplot(filen = filename, rundir = rundir, category = "Harvest_Rate", caption = caption)

  #### Finish ####
  endtime <- as.character(Sys.time())
  reportlist <- list(starttime = starttime, endtime = endtime)
  runnotes <- matrix(c(
    homelines,
    "Projection report. Built by Simon de Lestang. Relies on packages developed by Malcolm Haddon."
  ), ncol = 1)

  make_html(replist = reportlist, rundir = rundir, width = 800, openfile = openfile,
           runnotes = runnotes, verbose = FALSE, packagename = "makehtml",
           htmlname = "IMuLTProjection")

  invisible(NULL)
}
