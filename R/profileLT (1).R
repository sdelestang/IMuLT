#' True likelihood profile for an IMuLT model fit
#'
#' Fixes ONE parameter at a sequence of values and re-estimates every remaining
#' parameter with the full phased sandwich fit at each step (a true profile, not
#' a slice), recording the total negative log-likelihood and each component.
#'
#' This wraps your existing `FitModel()`, which is globals-driven and writes its
#' results to disk. The profiler therefore fixes the target the way FitModel
#' already understands - it sets that element's PHASE negative and its $Initial
#' to the profile value in the global `InitialVars` - then reads the component
#' likelihoods back out of `BigSave$Report`. All globals it touches (`Data`,
#' `InitialVars`) and the baseline `BigSave.lda` file are restored on exit.
#'
#' Requires in the global environment: `Data`, `InitialVars`, `FitModel`, and a
#' baseline `BigSave.lda` from a completed `FitModel(report = TRUE)` run.
#'
#' @param group Character name of the parameter block, e.g. "MainPars". Must be
#'              one of names(InitialVars) that maps to a PARAMETER_VECTOR.
#' @param pos   Integer 1-based position of the parameter within `group`.
#' @param lower,upper,step  Profile range on the parameter's ESTIMATION scale.
#'              Leave lower/upper NULL to default to mle +/- 3*`se`.
#' @param se    Optional SE (estimation scale) for the default range.
#' @param label Display / file name (e.g. the descriptive name "Rbar").
#'              Defaults to paste(group, pos).
#' @param rundir   Directory for the output text file and plot.
#' @param baseline Path to the BigSave.lda from the unprofiled MLE fit.
#' @param append   If TRUE, append to an existing "LPT <label>.txt"; rows carry a
#'              run timestamp and are de-duplicated (latest wins) on read.
#' @param plot   If TRUE, draw and save the component dNLL plot at the end.
#' @param fast   If TRUE, call FitModel(report = FALSE) and read the report from
#'              globals `ProfileReport`/`ProfileGrad` (skips the per-step
#'              sdreport). Requires the 3-line FitModel stash described in the
#'              accompanying notes. Default FALSE uses report = TRUE + BigSave.
#' @param ...   Optimisation-control arguments forwarded to FitModel (phit,
#'              lphit, mxph, newtonSteps, ...). Do NOT pass `report`.
#'
#' @return (invisibly) a data.frame, one row per profiled value.
profileLT <- function(group, pos,
                      lower = NULL, upper = NULL, step = NULL, se = NULL,
                      label = NULL, rundir = ".",
                      baseline = "Output/BigSave.lda",
                      append = FALSE, plot = TRUE, fast = FALSE, ...) {

  GE <- .GlobalEnv
  .blocks <- c("MainPars","RecruitPars","PuerPowPars","SelPars","RetPars",
               "RecDevs","Qpars","efpars","RecSpatDevs","MovePars","GrowthPars")

  ## ---- 0. Pull the globals FitModel relies on ---------------------------
  if (is.null(get0("Data", envir = GE)))        stop("global `Data` not found.")
  if (is.null(get0("InitialVars", envir = GE))) stop("global `InitialVars` not found.")
  if (!is.function(get0("FitModel", envir = GE))) stop("global `FitModel` not found.")
  if (!file.exists(baseline))
    stop("baseline `", baseline, "` not found - run FitModel(report = TRUE) first.")

  Data <- get("Data", envir = GE)
  IV   <- get("InitialVars", envir = GE)

  ## ---- 1. Validate the parameter choice ---------------------------------
  if (!group %in% names(IV))
    stop("`group` must be one of: ", paste(names(IV), collapse = ", "))
  blk_init <- IV[[group]]$Initial
  if (pos < 1L || pos > length(blk_init))
    stop(sprintf("`pos` must be in 1..%d for block '%s'.", length(blk_init), group))

  ## Refuse a MIRRORED parameter (value overwritten by its link inside the .cpp)
  linkmap <- c(MainPars = "MparsLink", RecruitPars = "RecparsLink",
               SelPars  = "SelparsLink", efpars = "EffparsLink",
               MovePars = "MoveparsLink")
  if (group %in% names(linkmap)) {
    lv <- Data[[ linkmap[[group]] ]]
    if (!is.null(lv) && length(lv) >= pos && !is.na(lv[pos])) {
      if (lv[pos] > 0)
        stop(sprintf(paste0("%s[%d] is a mirrored parameter (link = %d): its value ",
                            "is overwritten by the link target inside the model, so ",
                            "profiling it does nothing. Profile the target instead."),
                     group, pos, lv[pos]))
      if (lv[pos] < 0)
        message(sprintf("Note: %s[%d] is OFFSET-linked (link = %d); the profile shifts its offset.",
                        group, pos, lv[pos]))
    }
  }

  if (is.null(label)) label <- paste(group, pos)

  ## ---- 2. Baseline: MLE values, snapshot, file backup -------------------
  load_bigsave <- function(path) {
    e <- new.env(); load(path, envir = e)
    if (is.null(e$BigSave)) stop("no `BigSave` object in ", path)
    e$BigSave
  }
  relist_pin <- function(pin) {
    skel <- c(lapply(.blocks, function(b) IV[[b]]$Initial), list(dummy = 0))
    names(skel) <- c(.blocks, "dummy")
    if (length(pin) != length(unlist(skel)))
      stop("baseline pin length (", length(pin), ") != current parameter ",
           "structure (", length(unlist(skel)), "); InitialVars differs from the ",
           "fit that produced the baseline.")
    utils::relist(unname(pin), skel)
  }

  base    <- load_bigsave(baseline)
  base_pl <- relist_pin(base$pin)
  mle     <- base_pl[[group]][pos]

  base_tmp <- tempfile(fileext = ".lda")
  file.copy(baseline, base_tmp, overwrite = TRUE)
  Data_orig <- Data
  IV_orig   <- IV
  on.exit({
    assign("Data", Data_orig, envir = GE)
    assign("InitialVars", IV_orig, envir = GE)
    if (file.exists(base_tmp)) {            # restore the user's real fit
      file.copy(base_tmp, baseline, overwrite = TRUE); unlink(base_tmp)
    }
  }, add = TRUE)

  ## ---- 3. Build the value grid ------------------------------------------
  if (is.null(lower) || is.null(upper)) {
    if (is.null(se)) stop("Supply lower & upper, or `se` for a mle +/- 3*se range.")
    lower <- mle - 3 * se; upper <- mle + 3 * se
  }
  if (lower >= upper) stop("`lower` must be < `upper`.")
  if (is.null(step))  step <- (upper - lower) / 24
  if (step <= 0)      stop("`step` must be > 0.")
  grid <- seq(lower, upper, by = step)
  if (length(grid) < 3L) stop("Range/step yields fewer than 3 points.")

  ## ---- 4. Prepare the fixed-target globals ------------------------------
  ## Seed every block's $Initial at the MLE (near-optimum warm start for all
  ## grid points), force no projection so the component REPORTs exist, and pin
  ## the target by making its phase negative.
  for (b in .blocks) IV[[b]]$Initial <- base_pl[[b]]
  Data$DoProject <- 0L
  ph <- IV[[group]]$Phase[pos]
  IV[[group]]$Phase[pos] <- if (is.na(ph) || ph == 0) -1L else -abs(ph)

  ## ---- 5. Component extraction ------------------------------------------
  g <- function(rep, nm) if (!is.null(rep[[nm]])) as.numeric(rep[[nm]]) else NA_real_
  collect <- function(rep) c(
    total     = g(rep, "neglogL"),
    Catch     = g(rep, "CatchLike"),
    Cpue      = g(rep, "Weighted_CpueLike"),
    Numbers   = g(rep, "Weighted_NumbersLike"),
    Length    = g(rep, "Weighted_LengthLike"),
    Larval    = g(rep, "Weighted_LarvalLike"),
    Tag1      = g(rep, "Weighted_TagLike1"),
    Tag2      = g(rep, "Weighted_TagLike2"),
    RecPen    = g(rep, "Rec_Penal"),
    RecSmooth = g(rep, "Rec_Penal_Smooth"),
    SumZero   = g(rep, "Rec_Penal_SumZero"),
    InitPen   = g(rep, "Initial_pen"),
    MainPrior = g(rep, "MainParPriorPen"),
    RecPrior  = g(rep, "RecParPriorPen"),
    SelPrior  = g(rep, "SelParPriorPen"),
    EffPrior  = g(rep, "EffParPriorPen"))

  get_step <- function() {
    if (fast && exists("ProfileReport", envir = GE)) {
      list(rep  = get("ProfileReport", envir = GE),
           grad = if (exists("ProfileGrad", envir = GE))
                    max(abs(get("ProfileGrad", envir = GE))) else NA_real_)
    } else {
      bs <- load_bigsave(baseline)          # FitModel(report=TRUE) just overwrote it
      list(rep  = bs$Report,
           grad = if (!is.null(bs$Gradient)) max(abs(bs$Gradient)) else NA_real_)
    }
  }

  ## ---- 6. Profile loop (each step re-fit from the MLE seed) -------------
  res <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    IV[[group]]$Initial[pos] <- grid[i]
    assign("Data", Data, envir = GE)
    assign("InitialVars", IV, envir = GE)
    if (fast)
      suppressWarnings(rm(list = intersect(c("ProfileReport","ProfileGrad"), ls(GE)),
                          envir = GE))
    message(sprintf("Profiling %s = %.5g  (%d/%d)", label, grid[i], i, length(grid)))

    ok <- tryCatch({ FitModel(report = !fast, ...); TRUE },
                   error = function(e) { message("  step failed: ", conditionMessage(e)); FALSE })
    st   <- if (ok) tryCatch(get_step(), error = function(e) NULL) else NULL
    comp <- if (!is.null(st) && !is.null(st$rep)) collect(st$rep) else collect(NULL)
    res[[i]] <- data.frame(value = grid[i], t(comp),
                           maxgrad = if (!is.null(st)) st$grad else NA_real_,
                           check.names = FALSE)
  }

  prof <- do.call(rbind, res)
  prof <- prof[order(prof$value), ]
  prof$runtag    <- format(Sys.time(), "%Y%m%d_%H%M%S")
  prof$parameter <- label

  ## ---- 7. Write the text file (header once; append-safe) ----------------
  if (!dir.exists(rundir)) dir.create(rundir, recursive = TRUE)
  fn <- file.path(rundir, paste0("LPT ", label, ".txt"))
  do_append <- isTRUE(append) && file.exists(fn)
  utils::write.table(prof, file = fn, sep = "\t", row.names = FALSE,
                     col.names = !do_append, append = do_append, quote = FALSE)
  message("Wrote ", fn)

  ## ---- 8. Plot ----------------------------------------------------------
  if (plot) plotLT(fn, rundir = rundir, label = label)

  invisible(prof)
}


#' Plot a likelihood-profile output file produced by profileLT()
#'
#' Shows the change in each NLL component (relative to its own minimum across the
#' profiled range) against the parameter value, with the total in bold. The
#' dashed line at dNLL = 1.92 is the approximate 95% CI for the TOTAL on one
#' parameter (chi-square(1)/2) - valid because this is a true profile.
plotLT <- function(file, rundir = dirname(file), label = NULL,
                   keep_latest = TRUE, ci_line = TRUE) {

  d <- utils::read.delim(file, sep = "\t", check.names = FALSE)
  if (keep_latest && "runtag" %in% names(d))
    d <- d[d$runtag == max(d$runtag), , drop = FALSE]
  d <- d[!duplicated(d$value, fromLast = TRUE), ]
  d <- d[order(d$value), ]
  if (is.null(label))
    label <- if ("parameter" %in% names(d)) as.character(d$parameter[1]) else "parameter"

  comp_cols <- setdiff(names(d), c("value", "maxgrad", "runtag", "parameter"))
  long <- do.call(rbind, lapply(comp_cols, function(cc) {
    y <- suppressWarnings(as.numeric(d[[cc]]))
    if (all(is.na(y)) || diff(range(y, na.rm = TRUE)) < 1e-8) return(NULL)  # drop flat/empty
    data.frame(value = d$value, component = cc, dNLL = y - min(y, na.rm = TRUE))
  }))
  if (is.null(long)) stop("Nothing to plot - all components were flat or NA.")
  long$component <- factor(long$component,
                           levels = c("total", setdiff(unique(long$component), "total")))

  p <- ggplot2::ggplot(long, ggplot2::aes(value, dNLL, colour = component)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_line(data = subset(long, component == "total"),
                       colour = "black", linewidth = 1.3, show.legend = FALSE) +
    ggplot2::labs(x = paste0(label, "  (estimation scale)"),
                  y = expression(Delta * " negative log-likelihood"),
                  colour = NULL,
                  title = paste("Likelihood profile:", label)) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "right")
  if (ci_line && "total" %in% levels(long$component))
    p <- p + ggplot2::geom_hline(yintercept = 1.92, linetype = 2, colour = "grey50")

  fn <- file.path(rundir, paste0("LPT ", label, ".png"))
  ggplot2::ggsave(fn, p, width = 8, height = 5, dpi = 150)
  message("Wrote ", fn)
  invisible(list(data = d, long = long, plot = p, file = fn))
}
