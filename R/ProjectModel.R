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
                         folder_name = '', openfile = TRUE, ProjType = NULL) {

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

  # ── Optional ProjType override ──────────────────────────────────────────
  # Both the catch and harvest-rate schedules are always parsed into Data
  # regardless of what PROJECTIONS.DAT's own ProjType flag says (see
  # ReadProjFile()), so switching which one actually drives this run is safe
  # here without touching the file.
  if (!is.null(ProjType)) {
    if (!ProjType %in% c(1, 2))
      stop("ProjectModel: ProjType must be 1 (catch-based) or 2 ",
           "(harvest-rate-based).", call. = FALSE)
    cat("Overriding PROJECTIONS.DAT's ProjType (", Data$ProjType,
        ") with ProjType =", ProjType, "\n")
    Data$ProjType <- ProjType
  }

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

  last_proj <- Data$BurnIn + Data$Nyear + Data$Nproj
  relbio_proj <- sum(Report$LegalBioAll[last_proj,]) / sum(Report$VirginLegalBio)
  cat("Projected relative biomass (B/B0), final projection year:",
      round(relbio_proj, 3), "\n")

  # ── Save ─────────────────────────────────────────────────────────────────
  outdir <- dirname(bigsave_path)
  if (save) {
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

  # ── Report ───────────────────────────────────────────────────────────────
  if (report) {
    .MakeProjectionReport(Data = Data, Report = Report, ReportOrig = ReportOrig,
                          rundir_base = outdir, up = up, is95 = is95,
                          folder_name = folder_name, openfile = openfile)
  }

  invisible(list(Report = Report, ReportOrig = ReportOrig,
                 Data = Data, model = model))
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
#'
#' @return Invisibly NULL. Called for its side effect of writing the report.
#' @keywords internal
.MakeProjectionReport <- function(Data, Report, ReportOrig, rundir_base, up,
                                  is95, folder_name, openfile) {

  suppressPackageStartupMessages({
    library(makehtml); library(hplot); library(dplyr); library(magrittr)
    library(tidyr); library(ggplot2); library(reshape2); library(openxlsx)
  })
  options(dplyr.summarise.inform = FALSE)
  starttime <- as.character(Sys.time())
  SclErr <- ifelse(is95, 1.96, 0.842)

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

  #### Selectivity / Retention / Legal used in projections ####
  print("Making Selectivity/Retention used in projections")
  .plot_proj_curve <- function(PntFut, ActMat, label, category) {
    # PntFut: array (Nsex,Nage,Nfleet,MaxProjYr,Nstep); pull sex=1,age=1 (age
    # doesn't vary selectivity/retention pattern selection in this model) for
    # each fleet/proj-year/step, resolve against ActMat (0-indexed pointers).
    rows <- list()
    for (f in 1:Nfleet) {
      for (ys in seq_along(projyears)) {
        for (st in 1:Nstep) {
          ptr <- PntFut[1,1,f,ys,st]
          if (is.na(ptr) || ptr < 0) next
          curve <- ActMat[ptr + 1, ]
          rows[[length(rows)+1]] <- data.frame(
            fleet = f, year = projyears[ys], step = st, ptr = ptr,
            lbin = lbin, value = as.numeric(curve)[seq_along(lbin)])
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

    filename <- filenametopath(rundir, paste0(label, "_projected.png"))
    plotprep(width = 10, height = 10, filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1,1))
    print(ggplot(df, aes(x = lbin, y = value, colour = lab)) +
            geom_line(linewidth = 0.8) +
            facet_wrap(~descrip) +
            ylim(0, 1) +
            labs(x = "Length bin (mm)", y = label, colour = NULL) +
            theme_bw() + theme(legend.position = "bottom"))
    caption <- paste(label, "patterns applied during the projection years, by fleet.")
    addplot(filen = filename, rundir = rundir, category = category, caption = caption)
  }
  .plot_proj_curve(Data$SelPntFut, Report$ActSelex, "Selectivity", "Selectivity_Retention")
  .plot_proj_curve(Data$RetPntFut, Report$ActReten, "Retention",   "Selectivity_Retention")
  .plot_proj_curve(Data$LegalFleetPntFut, Report$ActLegal, "Legal size", "Selectivity_Retention")

  #### Catches (ProjType 1) or Harvest Rate schedule (ProjType 2) ####
  if (ProjType == 1) {
    print("Making projected Catch schedule")
    rows <- list()
    for (f in 1:Nfleet) for (yr in projyears) for (st in 1:Nstep) {
      v <- Data$Catch[.row_year1(yr), st, f]
      if (!is.na(v)) rows[[length(rows)+1]] <- data.frame(fleet=f, year=yr, step=st, catch=v)
    }
    catchdf <- bind_rows(rows) %>%
      mutate(areaname = fleetarea$areaname[match(fleet, fleetarea$fleet)]) %>%
      filter(catch > 0)
    if (nrow(catchdf) > 0) {
      filename <- filenametopath(rundir, "Projected_Catch.png")
      plotprep(width = 10, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
      parset(plots = c(1,1))
      print(ggplot(catchdf %>% group_by(year, areaname) %>% summarise(CatchT = sum(catch)/1000, .groups="drop"),
                   aes(x = year, y = CatchT, fill = areaname)) +
              geom_bar(stat = "identity", position = "stack") +
              viridis::scale_fill_viridis(discrete = TRUE) +
              ylab("Catch (t)") + xlab("Year") +
              theme(panel.background = element_rect(fill = "white", colour = NA),
                    panel.border = element_rect(fill = NA, colour = "grey20")))
      caption <- "Input catch schedule for the projection years, by area (ProjType 1)."
      addplot(filen = filename, rundir = rundir, category = "Catches", caption = caption)
    }
  } else {
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
              ylab("Input harvest rate") + xlab("Year") + theme_bw())
      caption <- "Input harvest-rate schedule for the projection years, by fleet (ProjType 2)."
      addplot(filen = filename, rundir = rundir, category = "Catches", caption = caption)
    }
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
    eff_rows[[length(eff_rows)+1]] <- data.frame(series = ef, year = yr, value = Report$CpueEcreep[ri, ef])
  }
  effdf <- bind_rows(eff_rows)
  filename <- filenametopath(rundir, "Projected_Fishing_Efficiency.png")
  plotprep(width = 8, height = 6, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  print(ggplot(effdf, aes(x = year, y = value, colour = factor(series))) +
          geom_line(linewidth = 0.8) + geom_point() +
          labs(x = "Year", y = "Fishing efficiency", colour = "Series") + theme_bw())
  caption <- "Fishing efficiency held at its last estimated value through the projection years."
  addplot(filen = filename, rundir = rundir, category = "Fishing_Efficiency", caption = caption)

  #### Biomass ####
  print("Making Biomass")
  bio_rows <- list()
  for (a in 1:Narea) for (yr in c(projyears[1]-1, projyears)) {
    ri <- .row_burnin(yr)
    bio_rows[[length(bio_rows)+1]] <- data.frame(
      area = a, year = yr, est = Report$LegalBioAll[ri, a],
      virgin = Report$VirginLegalBio[a])
  }
  biodf <- bind_rows(bio_rows) %>%
    mutate(areaname = areas$Name[match(area, areas$AreaCode)],
           `Legal Biomass (t)` = est/1000, `B/B0` = est/virgin)
  filename <- filenametopath(rundir, "Projected_Legal_Biomass.png")
  plotprep(width = 8, height = 6, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  print(ggplot(biodf %>% group_by(year) %>% summarise(est=sum(est), virgin=sum(virgin)) %>%
                 mutate(`B/B0`=est/virgin), aes(x = year, y = `B/B0`)) +
          geom_line(linewidth = 1) + geom_point() +
          geom_vline(xintercept = projyears[1]-0.5, linetype = "dashed", colour = "grey50") +
          geom_hline(yintercept = reflev$target, colour = 'green') +
          geom_hline(yintercept = reflev$threshold, colour = 'orange') +
          geom_hline(yintercept = reflev$limit, colour = 'red') +
          ylim(0, NA) +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20")) +
          ylab("B/B0") + xlab("Year"))
  caption <- "Total legal biomass relative to virgin, projected forward (dashed line = last assessment year)."
  addplot(filen = filename, rundir = rundir, category = "Biomass", caption = caption)

  filename <- filenametopath(rundir, "Projected_Legal_Biomass_Area.png")
  plotprep(width = 8, height = 8, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1,1))
  print(ggplot(biodf, aes(x = year, y = `B/B0`)) +
          geom_line() + geom_point(size = 0.75) +
          geom_vline(xintercept = projyears[1]-0.5, linetype = "dashed", colour = "grey50") +
          geom_hline(yintercept = reflev$target, colour = 'green') +
          geom_hline(yintercept = reflev$threshold, colour = 'orange') +
          geom_hline(yintercept = reflev$limit, colour = 'red') +
          facet_wrap(~areaname) + ylim(0, NA) +
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(angle = 45)))
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
          theme(panel.background = element_rect(fill = "white", colour = NA),
                panel.border = element_rect(fill = NA, colour = "grey20"),
                axis.text.x = element_text(angle = 45)) +
          ylab("Relative Mature Biomass") + xlab("Year"))
  caption <- paste("Relative mature biomass (egg production proxy -- MatureBioAllbySex, not the full",
                   "length-based egg-production curve, which does not yet extend into projection years),",
                   "by area, projected forward.")
  addplot(filen = filename, rundir = rundir, category = "Egg_Production", caption = caption)

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
