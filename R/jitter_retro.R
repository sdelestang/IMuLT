## ============================================================================
##  Jitter and retrospective diagnostics for IMuLT
##
##  Both drive the existing globals-driven FitModel(). They run it quietly
##  (report = FALSE, PrintNll = FALSE, console output swallowed) and harvest the
##  converged fit from the ProfileReport / ProfileGrad globals that FitModel
##  stashes on its final phase. report = FALSE means no BigSave.lda and no
##  Output.RL are written, so nothing accumulates across the many runs; the only
##  files FitModel still writes are Output/model<phase>.par, which overwrite in
##  place (fixed names) and therefore do not build up.
## ============================================================================

.PARBLOCKS <- c("MainPars","RecruitPars","PuerPowPars","SelPars","RetPars",
                "RecDevs","Qpars","efpars","RecSpatDevs","MovePars","GrowthPars")
.LINKMAP   <- c(MainPars = "MparsLink", RecruitPars = "RecparsLink",
                SelPars  = "SelparsLink", efpars = "EffparsLink",
                MovePars = "MoveparsLink")

## ---- Quiet single fit: returns report, max|grad|, convergence flag ---------
.quiet_fit <- function(grad_thresh = 0.1, ...) {
  GE <- .GlobalEnv
  suppressWarnings(rm(list = intersect(c("ProfileReport","ProfileGrad"), ls(GE)),
                      envir = GE))
  ok <- tryCatch({
    utils::capture.output(
      suppressWarnings(FitModel(report = FALSE, PrintNll = FALSE, ...)),
      file = nullfile())
    TRUE
  }, error = function(e) { message("   fit failed: ", conditionMessage(e)); FALSE })
  if (!ok || !exists("ProfileReport", envir = GE)) return(NULL)
  rep  <- get("ProfileReport", envir = GE)
  grad <- if (exists("ProfileGrad", envir = GE))
            max(abs(get("ProfileGrad", envir = GE))) else NA_real_
  list(rep = rep, grad = grad, converged = is.finite(grad) && grad < grad_thresh)
}


## ============================================================================
#' Jitter analysis for an IMuLT fit
#'
#' Perturbs the starting values of the active (positive-phase, non-mirrored)
#' parameters and re-fits `n` times through the full sandwich, to test whether
#' the model returns to the same optimum (convergence robustness / multimodality).
#'
#' @param n          Number of jittered fits.
#' @param jitter_sd  Perturbation size. With bounds, a fraction of (upr-lwr);
#'                   without bounds, a fraction of max(|init|, min_scale).
#' @param base_seed  Seed base; run i uses set.seed(base_seed + i) and records it.
#' @param grad_thresh max|grad| below which a run is flagged converged.
#' @param min_scale  Floor on the magnitude scale when no bounds are available.
#' @param bounds     If TRUE, look for bound fields in InitialVars and clip to them.
#' @param tol        ΔNLL within which a run counts as having found the best mode.
#' @param rundir,append,plot  Output controls.
#' @param ...        Forwarded to FitModel (phit, lphit, mxph, ...). Not `report`.
#'
#' @return (invisibly) list(summary, pars, best_nll, plot).
#'
#' @examples
#' \dontrun{
#' # 50 jittered re-fits, perturbing each start by 10% of its bound range
#' # (or magnitude where no bounds), then summarise convergence.
#' jit <- JitterFit(n = 50, jitter_sd = 0.1, mxph = MaxPhase)
#' table(jit$summary$converged)
#' jit$plot                       # ranked dNLL-from-best
#'
#' # Wider perturbation, fewer runs, written to a chosen directory.
#' JitterFit(n = 20, jitter_sd = 0.25, base_seed = 100,
#'           rundir = "Output/jitter", mxph = MaxPhase)
#' }
#' @name JitterFit
#' @export
JitterFit <- function(n = 50, jitter_sd = 0.1, base_seed = 1,
                      grad_thresh = 0.1, min_scale = 1, bounds = TRUE,
                      tol = 0.1, rundir = ".", append = FALSE, plot = TRUE, ...) {

  GE <- .GlobalEnv
  if (is.null(get0("InitialVars", envir = GE))) stop("global `InitialVars` not found.")
  if (is.null(get0("Data", envir = GE)))        stop("global `Data` not found.")
  if (!is.function(get0("FitModel", envir = GE))) stop("global `FitModel` not found.")

  IV_orig   <- get("InitialVars", envir = GE)
  Data_orig <- get("Data", envir = GE)
  on.exit({ assign("InitialVars", IV_orig, envir = GE)
            assign("Data", Data_orig, envir = GE) }, add = TRUE)

  Data <- Data_orig; Data$DoProject <- 0L; assign("Data", Data, envir = GE)
  base_IV <- IV_orig

  ## locate optional bound fields on a block
  ## bounds live in $Bnd: an [n, 2] matrix, col 1 = lower, col 2 = upper
  get_bounds <- function(b) {
    bnd <- base_IV[[b]]$Bnd
    if (is.null(bnd) || !is.matrix(bnd) || nrow(bnd) == 0) return(NULL)
    if (nrow(bnd) != length(base_IV[[b]]$Initial)) return(NULL)   # dimension guard
    list(lo = bnd[, 1], hi = bnd[, 2])
  }

  if (bounds && !any(vapply(.PARBLOCKS,
                            function(b) !is.null(get_bounds(b)), logical(1))))
    message("Note: no bound fields found on InitialVars blocks - jitter uses ",
            "magnitude scaling and starts are NOT clipped to bounds. Check the ",
            "bound field names or keep jitter_sd small.")

  ## active, non-mirrored positions per block
  active_idx <- function(b) {
    ph  <- base_IV[[b]]$Phase
    act <- which(!is.na(ph) & ph > 0)
    lv  <- if (b %in% names(.LINKMAP)) Data[[ .LINKMAP[[b]] ]] else NULL
    if (!is.null(lv)) act <- setdiff(act, which(!is.na(lv) & lv > 0))  # drop mirrors
    act
  }

  pull_pars <- function(rep) unlist(rep[.PARBLOCKS])

  ## ---- base (unperturbed) reference fit ----------------------------------
  assign("InitialVars", base_IV, envir = GE)
  message("Jitter: base (unperturbed) fit ...")
  b0 <- .quiet_fit(grad_thresh, ...)
  base_nll <- if (!is.null(b0)) as.numeric(b0$rep$neglogL) else NA_real_

  ## ---- jitter loop --------------------------------------------------------
  rows <- vector("list", n); pars <- vector("list", n)
  for (i in seq_len(n)) {
    seed_i <- base_seed + i
    set.seed(seed_i)
    newIV <- base_IV
    for (b in .PARBLOCKS) {
      act <- active_idx(b); if (!length(act)) next
      init <- base_IV[[b]]$Initial
      bd   <- if (bounds) get_bounds(b) else NULL
      sc   <- if (!is.null(bd)) (bd$hi - bd$lo)[act] else pmax(abs(init[act]), min_scale)
      init[act] <- init[act] + jitter_sd * sc * rnorm(length(act))
      if (!is.null(bd))
        init[act] <- pmin(pmax(init[act], bd$lo[act] + 1e-8), bd$hi[act] - 1e-8)
      newIV[[b]]$Initial <- init
    }
    assign("InitialVars", newIV, envir = GE)
    message(sprintf("Jitter run %d/%d (seed %d) ...", i, n, seed_i))
    fr <- .quiet_fit(grad_thresh, ...)
    rows[[i]] <- data.frame(run = i, seed = seed_i,
                            nll = if (!is.null(fr)) as.numeric(fr$rep$neglogL) else NA_real_,
                            maxgrad   = if (!is.null(fr)) fr$grad else NA_real_,
                            converged = if (!is.null(fr)) fr$converged else FALSE)
    pars[[i]] <- if (!is.null(fr)) pull_pars(fr$rep) else NULL
  }

  summ <- do.call(rbind, rows)
  best_nll <- min(c(base_nll, summ$nll), na.rm = TRUE)
  summ$dNLL <- summ$nll - best_nll
  summ$runtag <- format(Sys.time(), "%Y%m%d_%H%M%S")

  ## ---- write & report -----------------------------------------------------
  if (!dir.exists(rundir)) dir.create(rundir, recursive = TRUE)
  fn <- file.path(rundir, "Jitter.txt")
  do_append <- isTRUE(append) && file.exists(fn)
  utils::write.table(summ, fn, sep = "\t", row.names = FALSE,
                     col.names = !do_append, append = do_append, quote = FALSE)
  message("Wrote ", fn)

  nconv  <- sum(summ$converged, na.rm = TRUE)
  natbest <- sum(summ$dNLL <= tol & summ$converged, na.rm = TRUE)
  message(sprintf("Converged: %d/%d.  Reached best NLL (within %.3g): %d/%d.  Base NLL = %.4f, best = %.4f",
                  nconv, n, tol, natbest, n, base_nll, best_nll))
  if (any(summ$nll < base_nll - tol, na.rm = TRUE))
    message("** A jittered run beat the base fit by > tol - the base fit was NOT the global optimum. **")

  p <- NULL
  if (plot) {
    d <- summ[order(summ$nll), ]; d$rank <- seq_len(nrow(d))
    p <- ggplot2::ggplot(d, ggplot2::aes(rank, dNLL, colour = converged)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "red"),
                                   name = "converged") +
      ggplot2::labs(x = "run (ranked by NLL)",
                    y = expression(Delta * " NLL from best"),
                    title = sprintf("Jitter: %d runs, %d at best mode", n, natbest)) +
      ggplot2::theme_bw()
    pf <- file.path(rundir, "Jitter.png")
    ggplot2::ggsave(pf, p, width = 7, height = 5, dpi = 150)
    message("Wrote ", pf)
  }

  invisible(list(summary = summ, pars = pars, base_nll = base_nll,
                 best_nll = best_nll, plot = p))
}


## ============================================================================
#' Retrospective analysis for an IMuLT fit
#'
#' Sequentially removes the last `npeel` years, re-fits each peel through the
#' full sandwich, overlays a chosen quantity by calendar year, and computes
#' Mohn's rho on the terminal points.
#'
#' @param npeel   Number of years to peel (peels run 0..npeel; 0 = full model).
#' @param rebuild function(peel) that repopulates the globals (Data, InitialVars,
#'                ParOld, MaxPhase) for a model with the last `peel` years
#'                removed - typically a one-line wrapper around FileBuilder() with
#'                the terminal year reduced by `peel`. REQUIRED.
#' @param quantities Names of REPORTed vector quantities to track, e.g.
#'                "MatBio" (mature biomass) and "Recruits".
#' @param grad_thresh max|grad| convergence flag.
#' @param rundir,append,plot Output controls.
#' @param ...     Forwarded to FitModel. Not `report`.
#'
#' @return (invisibly) list(series, rho, plot).
#' @examples
#' \dontrun{
#' # Peel the last 5 years. `rebuild` re-runs FileBuilder with the terminal
#' # year reduced by `peel` - adjust the FileBuilder argument name to your setup.
#' full_terminal <- Data$Year1 + Data$Nyear - 1
#' retro <- RetroFit(
#'   npeel      = 5,
#'   rebuild    = function(peel) FileBuilder(LastYear = full_terminal - peel),
#'   quantities = c("MatBio", "Recruits"),
#'   mxph       = MaxPhase)
#' retro$rho                      # Mohn's rho per quantity
#' retro$plot                     # peel overlay, rho in the facet strip labels
#' }
#' @name RetroFit
#' @export
RetroFit <- function(npeel = 5, rebuild = NULL,
                     quantities = c("MatBio", "Recruits"),
                     grad_thresh = 0.1, rundir = ".", append = FALSE,
                     plot = TRUE, ...) {

  GE <- .GlobalEnv
  owd <- getwd(); on.exit(setwd(owd), add = TRUE)   # restore wd even if a peel errors
  if (!is.function(rebuild))
    stop("Provide `rebuild`: a function(peel) ...")
  if (!is.function(rebuild))
    stop("Provide `rebuild`: a function(peel) that repopulates the globals ",
         "(Data, InitialVars, ParOld, MaxPhase) for a model with the last ",
         "`peel` years removed - e.g. a wrapper around FileBuilder() with the ",
         "terminal year reduced by `peel`.")

  ## snapshot the four globals rebuild() will overwrite, restore on exit
  snap <- mget(c("Data","InitialVars","ParOld","MaxPhase"),
               envir = GE, ifnotfound = list(NULL))
  on.exit(for (nm in names(snap))
            if (!is.null(snap[[nm]])) assign(nm, snap[[nm]], envir = GE), add = TRUE)

  ## assessment-period series (calendar year vs value) from a REPORTed vector
  series <- function(rep, nm, Data) {
    v <- rep[[nm]]; if (is.null(v)) return(NULL)
    v  <- as.numeric(v)
    By <- Data$BurnIn; Y1 <- Data$Year1; Ny <- Data$Nyear
    r  <- (By + 1):(By + Ny)                       # R indices of assessment years
    data.frame(year = Y1 + (r - 1) - By, value = v[r])
  }

  ## ---- peel loop ----------------------------------------------------------
  store <- list()                                   # store[[quantity]][[peel+1]]
  for (q in quantities) store[[q]] <- vector("list", npeel + 1L)

  for (p in 0:npeel) {
    message(sprintf("Retro peel %d/%d ...", p, npeel))
    ok <- tryCatch({ rebuild(p); TRUE },
                   error = function(e) { message("   peel ", p, " rebuild failed: ",
                                                 conditionMessage(e)); FALSE })
    if (!ok) next
    Data <- get("Data", envir = GE); Data$DoProject <- 0L
    assign("Data", Data, envir = GE)                                    # repopulates the globals
    fr <- .quiet_fit(grad_thresh, ...)
    if (is.null(fr)) { message("   peel ", p, " did not return a fit - skipping."); next }
    if (!fr$converged)
      message(sprintf("   warning: peel %d max|grad| = %.4g (>= %.3g)",
                      p, fr$grad, grad_thresh))
    for (q in quantities) {
      s <- series(fr$rep, q, Data)
      if (!is.null(s)) { s$peel <- p; s$quantity <- q; store[[q]][[p + 1L]] <- s }
    }
  }

  setwd(owd)   # back to launch dir: Retro.txt / Retro.png write here, not the last peel folder

  ## ---- assemble long series & Mohn's rho ---------------------------------
  long <- do.call(rbind, unlist(store, recursive = FALSE))
  long$runtag <- format(Sys.time(), "%Y%m%d_%H%M%S")

  mohn <- function(q) {
    byp  <- store[[q]]
    full <- byp[[1]]; if (is.null(full)) return(NA_real_)
    cs <- vapply(seq_len(npeel), function(p) {
      sp <- byp[[p + 1L]]; if (is.null(sp)) return(NA_real_)
      Tp <- max(sp$year)
      xp <- sp$value[sp$year == Tp]
      x0 <- full$value[full$year == Tp]
      if (length(x0) != 1 || x0 == 0) return(NA_real_)
      (xp - x0) / x0
    }, numeric(1))
    mean(cs, na.rm = TRUE)
  }
  rho <- vapply(quantities, mohn, numeric(1))

  ## ---- write & report -----------------------------------------------------
  if (!dir.exists(rundir)) dir.create(rundir, recursive = TRUE)
  fn <- file.path(rundir, "Retro.txt")
  do_append <- isTRUE(append) && file.exists(fn)
  utils::write.table(long, fn, sep = "\t", row.names = FALSE,
                     col.names = !do_append, append = do_append, quote = FALSE)
  message("Wrote ", fn)
  for (q in quantities) message(sprintf("Mohn's rho [%s] = %+.4f", q, rho[[q]]))

  p_obj <- NULL
  if (plot) {
    long$peel <- factor(long$peel)
    labs <- setNames(sprintf("%s  (rho = %+.3f)", quantities, rho[quantities]), quantities)
    p_obj <- ggplot2::ggplot(long, ggplot2::aes(year, value, colour = peel, group = peel)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ quantity, scales = "free_y",
                          labeller = ggplot2::labeller(quantity = labs)) +
      ggplot2::labs(x = "Year", y = NULL, colour = "peel",
                    title = "Retrospective analysis") +
      ggplot2::theme_bw()
    pf <- file.path(rundir, "Retro.png")
    ggplot2::ggsave(pf, p_obj, width = 9, height = 5, dpi = 150)
    message("Wrote ", pf)
  }

  invisible(list(series = long, rho = rho, plot = p_obj))
}
