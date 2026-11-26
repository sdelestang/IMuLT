#' Parallel likelihood profile for an IMuLT model fit
#'
#' Runs ProfileLT() over several CPU cores at once. The grid of profile values
#' is split into contiguous chunks, one per worker. Each worker is a separate R
#' process working in its OWN copy of the run folder, because FitModel() reads
#' and writes Output/BigSave.lda and Output/model final.par in the working
#' directory and workers sharing one folder would overwrite each other.
#' Within a chunk ProfileLT()'s warm start still applies (each step starts from
#' the previous step's solution); the first step of every chunk starts from the
#' MLE. When all workers finish, their "LPT <label>.txt" (and states) files are
#' merged into one file in `rundir` and plotted with plotLT().
#'
#' Requires in the global environment: `Data`, `InitialVars` (plus anything else
#' named in `export`), and a baseline BigSave.lda from a completed
#' FitModel(report = TRUE) run. The current working directory must be the run
#' folder (the one FitModel() is normally run from).
#'
#' @param group,pos,lower,upper,step,se,values,label As for ProfileLT(). With
#'              few cores, a short explicit `values` vector is usually best.
#' @param ncores  Number of worker processes (capped at the number of chunks).
#'              Each worker holds a full copy of the model in memory.
#' @param chunk   Grid points per worker. Default spreads the grid evenly over
#'              `ncores`. Use chunk = 1 for one point per worker (no warm chain).
#' @param rundir  Directory for the merged output text file and plot.
#' @param baseline Path (relative to the run folder) to the MLE BigSave.lda.
#' @param workroot Where the per-worker folder copies go. Defaults to the
#'              session temp folder, which keeps the heavy file traffic OUT of
#'              OneDrive. Must not be the run folder itself.
#' @param export  Names of global objects copied to every worker.
#' @param setup   Expression evaluated on every worker before fitting, e.g.
#'              quote(library(IMuLT)). Use the INSTALLED package; load_all() on
#'              many workers at once will try to recompile in parallel.
#' @param keep_work If TRUE, keep the worker folders (useful for debugging).
#' @param append,plot,fast,newtonSteps,save_states,warm As for ProfileLT().
#' @param ...     Optimisation-control arguments forwarded to FitModel via
#'              ProfileLT (phit, lphit, mxph, ...). Do NOT pass `report`.
#'
#' @return (invisibly) the merged data.frame, one row per profiled value.
#'
#' @examples
#' \dontrun{
#' # Rbar at 8 points, 4 cores, 2 points per worker:
#' ProfileLTpar("MainPars", pos = 1,
#'              values = 17.86 + c(-0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6),
#'              ncores = 4, label = "Rbar", fast = TRUE, mxph = MaxPhase)
#' }
#' @name ProfileLTpar
#' @export
ProfileLTpar <- function(group, pos,
                         lower = NULL, upper = NULL, step = NULL, se = NULL,
                         values = NULL, label = NULL,
                         ncores = max(1L, parallel::detectCores() - 1L),
                         chunk = NULL, rundir = ".",
                         baseline = "Output/BigSave.lda",
                         workroot = file.path(tempdir(), "LPTwork"),
                         export = c("Data", "InitialVars"),
                         setup = quote(library(IMuLT)),
                         keep_work = FALSE,
                         append = FALSE, plot = TRUE, fast = FALSE, newtonSteps = 0,
                         save_states = TRUE, warm = TRUE, ...) {

  GE <- .GlobalEnv
  .blocks <- c("MainPars","RecruitPars","PuerPowPars","SelPars","RetPars",
               "RecDevs","Qpars","efpars","RecSpatDevs","MovePars","GrowthPars",
               "TemporalMPars")

  ## ---- 0. Checks ---------------------------------------------------------
  for (nm in export)
    if (!exists(nm, envir = GE)) stop("global `", nm, "` not found.")
  if (!file.exists(baseline))
    stop("baseline `", baseline, "` not found - run FitModel(report = TRUE) first.")
  IV <- get("InitialVars", envir = GE)
  if (!group %in% names(IV))
    stop("`group` must be one of: ", paste(names(IV), collapse = ", "))
  if (is.null(label)) label <- paste(group, pos)
  rundir <- normalizePath(rundir, mustWork = FALSE)
  src    <- normalizePath(getwd())

  ## ---- 1. Grid (needs the MLE only for the se-based default) -------------
  if (!is.null(values)) {
    grid <- sort(unique(as.numeric(values)))
  } else {
    if (is.null(lower) || is.null(upper)) {
      if (is.null(se)) stop("Supply values, lower & upper, or `se`.")
      e <- new.env(); load(baseline, envir = e)
      skel <- c(lapply(.blocks, function(b) IV[[b]]$Initial), list(dummy = 0))
      names(skel) <- c(.blocks, "dummy")
      mle <- utils::relist(unname(e$BigSave$pin), skel)[[group]][pos]
      lower <- mle - 3 * se; upper <- mle + 3 * se
    }
    if (is.null(step)) step <- (upper - lower) / 24
    grid <- seq(lower, upper, by = step)
  }
  if (!length(grid)) stop("No profile values.")

  ## ---- 2. Split into contiguous chunks -----------------------------------
  if (is.null(chunk)) chunk <- ceiling(length(grid) / max(1L, ncores))
  chunks <- split(grid, ceiling(seq_along(grid) / chunk))
  ncores <- min(ncores, length(chunks))
  message(sprintf("Profiling %s at %d values: %d chunks on %d cores.",
                  label, length(grid), length(chunks), ncores))

  ## ---- 3. One folder copy per chunk --------------------------------------
  dir.create(workroot, recursive = TRUE, showWarnings = FALSE)
  workroot <- normalizePath(workroot)
  if (identical(workroot, src))
    stop("`workroot` must not be the run folder itself.")
  top <- list.files(src, full.names = TRUE, all.files = FALSE)
  top <- top[!file.info(top)$isdir]                 # input files only, no sub-folders
  extra <- c(baseline, file.path("Output", "model final.par"))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  jobs <- lapply(seq_along(chunks), function(k) {
    wd <- file.path(workroot, sprintf("%s_%s_w%02d", gsub("[^A-Za-z0-9]", "", label), stamp, k))
    dir.create(wd, recursive = TRUE, showWarnings = FALSE)
    file.copy(top, wd, overwrite = TRUE)
    for (x in extra) if (file.exists(x)) {
      dir.create(dirname(file.path(wd, x)), recursive = TRUE, showWarnings = FALSE)
      file.copy(x, file.path(wd, x), overwrite = TRUE)
    }
    list(k = k, wd = wd, values = chunks[[k]])
  })
  if (!keep_work)
    on.exit(unlink(vapply(jobs, `[[`, "", "wd"), recursive = TRUE), add = TRUE)

  ## ---- 4. Start workers ---------------------------------------------------
  logf <- file.path(workroot, paste0("LPT ", label, " ", stamp, " workers.log"))
  cl <- parallel::makeCluster(ncores, outfile = logf)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  message("Worker console output: ", logf)
  parallel::clusterExport(cl, export, envir = GE)
  parallel::clusterCall(cl, function(ex) {
    eval(ex, envir = .GlobalEnv)
    if (requireNamespace("TMB", quietly = TRUE))
      try(TMB::openmp(1), silent = TRUE)            # 1 thread per worker: no oversubscription
    invisible(TRUE)
  }, setup)

  args <- c(list(group = group, pos = pos, label = label, rundir = ".",
                 baseline = baseline, append = FALSE, plot = FALSE, fast = fast,
                 newtonSteps = newtonSteps, save_states = save_states, warm = warm),
            list(...))

  ## ---- 5. Run -------------------------------------------------------------
  t0 <- Sys.time()
  out <- parallel::clusterApplyLB(cl, jobs, function(job, args) {
    setwd(job$wd)
    PL <- get("ProfileLT", envir = .GlobalEnv)      # works for library() or load_all()
    r <- tryCatch({ do.call(PL, c(args, list(values = job$values))); "ok" },
                  error = function(e) conditionMessage(e))
    list(k = job$k, wd = job$wd, status = r)
  }, args = args)
  message(sprintf("All workers finished in %.1f min.",
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  bad <- Filter(function(o) o$status != "ok", out)
  for (o in bad) message(sprintf("  chunk %d failed: %s", o$k, o$status))

  ## ---- 6. Merge -----------------------------------------------------------
  rd <- function(f) if (file.exists(f)) utils::read.delim(f, sep = "\t", check.names = FALSE) else NULL
  prof <- do.call(rbind, lapply(out, function(o) rd(file.path(o$wd, paste0("LPT ", label, ".txt")))))
  if (is.null(prof) || !nrow(prof)) stop("No worker produced results; see ", logf)
  prof <- prof[order(prof$value), ]
  prof$runtag <- stamp

  if (!dir.exists(rundir)) dir.create(rundir, recursive = TRUE)
  fn <- file.path(rundir, paste0("LPT ", label, ".txt"))
  do_app <- isTRUE(append) && file.exists(fn)
  utils::write.table(prof, fn, sep = "\t", row.names = FALSE,
                     col.names = !do_app, append = do_app, quote = FALSE)
  message("Wrote ", fn)

  if (save_states) {
    sx <- do.call(rbind, lapply(out, function(o) rd(file.path(o$wd, paste0("LPT ", label, " states.txt")))))
    if (!is.null(sx)) {
      sx$runtag <- stamp
      fs <- file.path(rundir, paste0("LPT ", label, " states.txt"))
      da <- isTRUE(append) && file.exists(fs)
      utils::write.table(sx, fs, sep = "\t", row.names = FALSE,
                         col.names = !da, append = da, quote = FALSE)
      message("Wrote ", fs)
    }
  }

  ## ---- 7. Plot ------------------------------------------------------------
  if (plot) plotLT(fn, rundir = rundir, label = label)
  invisible(prof)
}
