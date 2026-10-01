# ── Main function ─────────────────────────────────────────────────────────────
#' Fit the IMuLT Stock Assessment Model
#'
#' Fits the IMuLT TMB model using a phased optimisation approach. Parameters are
#' progressively introduced across phases. In the final phase, nlminb is run
#' until the gradient is small enough (pre-sandwich), then sandwich restarts
#' (alternating BFGS on logit-transformed parameters and nlminb) polish the
#' solution. Optional Newton polishing steps can further refine it.
#'
#' @param phit Integer. Maximum number of function evaluations per phase for
#'   all phases except the last. Default 500.
#' @param lphit Integer. Maximum number of function evaluations for the last
#'   phase (per nlminb call). Default 1000.
#' @param mxph Integer. Maximum phase number. If 0, treated as 1. Defaults to
#'   the global \code{MaxPhase}.
#' @param PrintLag Integer. Print and plot progress every \code{PrintLag}
#'   function evaluations (including during BFGS). Default 50.
#' @param report Logical. If \code{TRUE}, produce SD report and save outputs
#'   after the final phase. Default \code{FALSE}.
#' @param nRestarts Logical. If \code{TRUE}, sandwich restarts run in the final
#'   phase. Default \code{TRUE}.
#' @param newtonSteps Integer. Maximum number of damped Newton polishing steps
#'   after the sandwich. Default 5.
#' @param newton_grad_thresh Numeric. Gradient threshold below which Newton
#'   polishing is considered safe; also the sandwich convergence criterion.
#'   Default 0.1.
#' @param sandwich_entry_grad Numeric. Before the sandwich starts, nlminb is
#'   re-run from its last point until max|grad| <= this value (or it stalls, or
#'   \code{max_pre_nlminb} calls are used). Also the floor for the BFGS
#'   gradient-inflation rejection. Default 10.
#' @param max_pre_nlminb Integer. Maximum number of extra nlminb calls before
#'   the sandwich. Default 5.
#' @param bfgs_maxit Integer. Maximum BFGS iterations per sandwich cycle.
#'   Default 200.
#' @param checkpoint Logical. If \code{TRUE}, save the current final-phase
#'   parameters to \code{Output/checkpoint_par.rds} after every nlminb call.
#'   Default \code{TRUE}.
#' @param PrintNll Logical. If \code{TRUE}, display a live trace plot of the
#'   negative log-likelihood. Default \code{TRUE}.
#' @param plateau_k Integer. Window size (in reporting intervals) for the
#'   plateau diagnostic. Default 4.
#' @param plateau_cv Numeric. CV threshold for the plateau diagnostic in
#'   non-final phases. Set to 0 to disable. Default 0 (off) - a steady rate of
#'   slow improvement also has a low CV, so the diagnostic can mislead.
#'
#' @details
#' \strong{Sandwich restarts:} BFGS operates in logit-transformed space, where
#' parameters near a bound have tiny chain-rule gradients and can be pushed onto
#' the bound. A BFGS result is therefore rejected if its gradient is non-finite,
#' or if it inflates max|grad| more than 5-fold (and above
#' \code{sandwich_entry_grad}); nlminb then continues from the pre-BFGS point.
#'
#' \strong{Gradient wrapper:} returns zeros for non-finite gradients to keep
#' optimisers alive; every accepted point is checked with \code{.max_grad()}.
#'
#' @export
FitModel <- function(phit = 500, lphit = 1000, mxph = MaxPhase,
                     PrintLag = 50, report = FALSE,
                     nRestarts = TRUE, newtonSteps = 5,
                     newton_grad_thresh = 0.1,
                     sandwich_entry_grad = 10, max_pre_nlminb = 5,
                     bfgs_maxit = 200, checkpoint = TRUE,
                     PrintNll = TRUE,
                     plateau_k = 4, plateau_cv = 0) {

  ##  globals populated by LoadData() / LoadPars() --------
  need <- c("Data", "InitialVars")
  miss <- need[!vapply(need, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)]
  ## mxph defaults to the MaxPhase global; only require it if not supplied
  if (missing(mxph) && !exists("MaxPhase", envir = .GlobalEnv, inherits = FALSE))
    miss <- c(miss, "MaxPhase")
  if (length(miss)) {
    warning("FitModel: ", paste(miss, collapse = ", "),
            " not found - run LoadData() then LoadPars() first.", call. = FALSE)
    return(invisible(NULL))
  }

  MaxPhase <- ifelse(mxph == 0, 1, mxph)

  ## Force negative phase for mirrored parameters ##
  link_map <- list(MainPars    = Data$MparsLink,
                   RecruitPars = Data$RecparsLink,
                   SelPars     = Data$SelparsLink,
                   efpars      = Data$EffparsLink,
                   GrowthPars  = Data$GrowparsLink)
  for (pname in names(link_map)) {
    lv <- link_map[[pname]]
    if (!is.null(lv)) {
      for (i in seq_along(lv)) {
        if (!is.na(lv[i]) && lv[i] > 0 && InitialVars[[pname]]$Phase[i] > 0) {
          warning(paste(pname, "parameter", i,
                        "is mirrored but has positive phase — forcing negative"))
          InitialVars[[pname]]$Phase[i] <- -abs(InitialVars[[pname]]$Phase[i])
        }
      }
    }
  }

  assign("InitialVars", InitialVars, envir = .GlobalEnv)   # keep global phases consistent with the fit

  ## Track ##
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

  ## Logit transform  (bounded <-> unconstrained) ##
  .to_unbounded <- function(x, lo, hi) {
    x_c <- pmax(pmin(x, hi - 1e-8), lo + 1e-8)
    log((x_c - lo) / (hi - x_c))
  }

  .to_bounded <- function(y, lo, hi) {
    lo + (hi - lo) / (1 + exp(-y))
  }

  # Chain-rule correction: dL/dy = dL/dx * dx/dy
  .chain_grad <- function(g, y, lo, hi) {
    sig  <- 1 / (1 + exp(-y))
    dxdy <- (hi - lo) * sig * (1 - sig)
    g * dxdy
  }

  ## Phase loop ##
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

    ## Plateau diagnostic state (non-final phases; off by default) ##
    delta_history  <- numeric(0)
    plateau_eval   <- NA_integer_
    min_evals_diag <- PrintLag * (plateau_k + 2)

    ## Print suppression flag ##
    # Only set TRUE around Newton Hessian calculations; BFGS and nlminb print progress.
    suppress_print <- FALSE

    ## Store originals then wrap both fn and gr ##
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
          cat(CurrentStage, " ", FnCallNo, " -LogLike: ",
              round(tyy, 3), " | Delta: ", round(abs(delta), 6), "%\n", sep = "")
          LastPrintFn <<- tyy
          .trace_append(TotalEval + FnCallNo, tyy)
          if (PrintNll) {
            dev.hold()
            .plot_trace(TraceDF, tyy)
            dev.flush()
          }

          ## Plateau diagnostic (non-final phases only, if enabled) ##
          if (!is_final && plateau_cv > 0 && FnCallNo >= min_evals_diag) {
            delta_history <<- c(delta_history, abs(delta))
            if (length(delta_history) >= plateau_k && is.na(plateau_eval)) {
              recent <- tail(delta_history, plateau_k)
              cv     <- sd(recent) / abs(mean(recent))
              if (cv < plateau_cv) plateau_eval <<- FnCallNo
            }
          }
        }
      }
      return(tyy)
    }

    # Gradient wrapper: zeros for non-finite gradients (keeps optimisers alive).
    # A zero gradient looks like convergence, so accepted points are checked with .max_grad().
    nonfinite_gr <- 0L
    model$gr <- function(x) {
      g <- tryCatch(model$gr_Orig(x), error = function(e) NULL)
      if (is.null(g) || any(!is.finite(g))) {
        nonfinite_gr <<- nonfinite_gr + 1L
        return(rep(0, length(x)))
      }
      return(g)
    }

    # Max |gradient| at x; Inf if non-finite or error.
    .max_grad <- function(x) {
      g <- tryCatch(model$gr_Orig(x), error = function(e) NA_real_)
      if (length(g) == 0 || any(!is.finite(g))) return(Inf)
      max(abs(g))
    }

    # Save current final-phase parameters so an interrupted run isn't wasted
    .checkpoint <- function(m, label) {
      if (is_final && isTRUE(checkpoint))
        saveRDS(list(par = m$par, pnames = pnames, objective = m$objective,
                     label = label, time = Sys.time()),
                "Output/checkpoint_par.rds")
    }

    # Reset the progress-print state at the start of an optimiser call
    .start_stage <- function(label, start_value) {
      CurrentStage <<- label
      FnCallNo     <<- 0
      BestFn       <<- start_value
      LastPrintFn  <<- start_value
      .trace_append(TotalEval, start_value)
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

    if (!is.finite(.max_grad(mout$par)))
      cat("  WARNING: initial nlminb ended at a point with a non-finite gradient (",
          nonfinite_gr, " non-finite gradient calls so far)\n", sep = "")

    ## Post-hoc plateau diagnostic (non-final phases only) ##
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
    .checkpoint(mout, "initial nlminb")

    ## Pre-sandwich nlminb: keep going until the gradient is small enough for the sandwich ##
    if (is_final && isTRUE(nRestarts)) {
      entry_grad <- .max_grad(mout$par)
      for (k in seq_len(max_pre_nlminb)) {
        if (!is.finite(entry_grad) || entry_grad <= sandwich_entry_grad) break
        cat("\n  Pre-sandwich nlminb", k, "- starting max|grad| =", round(entry_grad, 3), "\n")
        prev_obj <- mout$objective
        .start_stage(paste0("Pre-sandwich nlminb ", k), prev_obj)
        mout_k <- nlminb(mout$par, model$fn, model$gr,
                         lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                         control = list(iter.max = MaXeVaL, eval.max = MaXeVaL,
                                        rel.tol = 1e-12, x.tol = 1e-12, abs.tol = 0))
        TotalEval <<- TotalEval + FnCallNo
        new_grad <- .max_grad(mout_k$par)
        if (!is.finite(new_grad)) {
          cat("  Pre-sandwich nlminb ended on a non-finite gradient — keeping previous point\n")
          break
        }
        mout <- mout_k
        .report_fit(mout, model, pnames, prev_obj, label = paste("  Pre-sandwich nlminb", k))
        .checkpoint(mout, paste("pre-sandwich nlminb", k))
        if (prev_obj - mout$objective < 1e-6 && new_grad >= entry_grad) {
          cat("  Pre-sandwich nlminb stalled — moving to sandwich\n")
          entry_grad <- new_grad
          break
        }
        entry_grad <- new_grad
      }
    }

    ## Sandwich restarts (final phase only) ##
    if (is_final && isTRUE(nRestarts)) {

      sandwich_done <- FALSE
      no_improve    <- 0L                    # consecutive cycles with no grad improvement
      max_sandwich  <- 20L                   # safety ceiling
      prev_grad     <- .max_grad(mout$par)   # gradient at sandwich entry

      for (restart in seq_len(max_sandwich)) {

        cat("\n--- Sandwich restart", restart, "---\n")

        bfgs_reltol  <- ifelse(restart == 1, 1e-12, 1e-15)
        nlminb_rtol  <- ifelse(restart == 1, 1e-12, 1e-15)
        nlminb_xtol  <- ifelse(restart == 1, 1e-12, 1e-15)

        ## Step A: BFGS in logit-transformed unconstrained space ##
        cat("  Step A: BFGS (logit-transformed, max", bfgs_maxit, "iterations)\n")

        # last.par.best can be a point with a non-finite gradient - only use it if fn and gradient are finite
        bfgs_start     <- model$env$last.par.best
        fn_check       <- model$fn(bfgs_start)
        pre_bfgs_grad  <- .max_grad(bfgs_start)
        if (!is.finite(fn_check) || !is.finite(pre_bfgs_grad)) {
          cat("  WARNING: last.par.best has non-finite objective or gradient, falling back to mout$par\n")
          bfgs_start    <- mout$par
          fn_check      <- model$fn(bfgs_start)
          pre_bfgs_grad <- .max_grad(bfgs_start)
        }
        .start_stage(paste0("Restart ", restart, " \u2013 BFGS"), fn_check)

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
                            control = list(maxit  = bfgs_maxit,
                                           reltol = bfgs_reltol))

          fit_bfgs$par <- .to_bounded(fit_bfgs$par, lo, hi)

        } else {
          fit_bfgs <- optim(bfgs_start, model$fn, model$gr,
                            method  = "BFGS",
                            control = list(maxit  = bfgs_maxit,
                                           reltol = bfgs_reltol))
        }

        bfgs_grad <- .max_grad(fit_bfgs$par)
        cat("  BFGS complete: obj =", round(fit_bfgs$value, 6),
            "| max|grad| =", round(bfgs_grad, 6),
            "| convergence:", fit_bfgs$convergence,
            "| fn evals:", FnCallNo, "\n")
        TotalEval <<- TotalEval + FnCallNo

        ## BFGS acceptance guards ##
        if (!is.finite(bfgs_grad)) {
          # Non-finite gradient: nlminb can't polish from here (gr wrapper returns zeros)
          cat("  BFGS ended at a point with a non-finite gradient — reverting to pre-BFGS parameters\n")
          fit_bfgs$par   <- bfgs_start
          fit_bfgs$value <- fn_check
          bfgs_grad      <- pre_bfgs_grad

        } else if (bfgs_grad > 5 * pre_bfgs_grad && bfgs_grad > sandwich_entry_grad) {
          # Gradient badly inflated - typically parameters pushed onto bounds in logit space
          cat("  BFGS inflated max|grad| (", round(pre_bfgs_grad, 3), "->", round(bfgs_grad, 3),
              ") — rejecting and keeping pre-BFGS parameters\n", sep = "")
          fit_bfgs$par   <- bfgs_start
          fit_bfgs$value <- fn_check
          bfgs_grad      <- pre_bfgs_grad

        } else {
          bfgs_nll_improved  <- fit_bfgs$value < fn_check
          bfgs_grad_improved <- bfgs_grad <= pre_bfgs_grad

          if (!bfgs_nll_improved && !bfgs_grad_improved) {
            cat("  BFGS degraded both NLL (", round(fn_check, 4), "->",
                round(fit_bfgs$value, 4), ") and gradient (",
                round(pre_bfgs_grad, 4), "->", round(bfgs_grad, 4),
                ") — reverting to pre-BFGS parameters\n", sep = "")
            fit_bfgs$par   <- bfgs_start
            fit_bfgs$value <- fn_check
            bfgs_grad      <- pre_bfgs_grad
          } else {
            if (bfgs_nll_improved && !bfgs_grad_improved) {
              cat("  BFGS improved NLL (", round(fn_check, 4), "->",
                  round(fit_bfgs$value, 4), ") but degraded gradient (",
                  round(pre_bfgs_grad, 4), "->", round(bfgs_grad, 4),
                  ") — keeping, nlminb will polish\n", sep = "")
            } else if (!bfgs_nll_improved && bfgs_grad_improved) {
              cat("  BFGS improved gradient (", round(pre_bfgs_grad, 4), "->",
                  round(bfgs_grad, 4), ") but degraded NLL (",
                  round(fn_check, 4), "->", round(fit_bfgs$value, 4),
                  ") — keeping\n", sep = "")
            } else {
              cat("  BFGS improved both NLL (", round(fn_check, 4), "->",
                  round(fit_bfgs$value, 4), ") and gradient (",
                  round(pre_bfgs_grad, 4), "->", round(bfgs_grad, 4),
                  ")\n", sep = "")
            }
            pre_bfgs_grad <- bfgs_grad
          }
        }

        ## Step B: nlminb from BFGS solution ##
        cat("  Step B: nlminb\n")
        .start_stage(paste0("Restart ", restart, " \u2013 nlminb"), fit_bfgs$value)
        initBestFn <- fit_bfgs$value

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

        cur_grad <- .max_grad(mout$par)

        ## nlminb non-finite gradient guard ##
        if (!is.finite(cur_grad)) {
          cat("  WARNING: nlminb ended at a point with a non-finite gradient (",
              nonfinite_gr, " non-finite gradient calls this phase) — reverting to BFGS point and ending sandwich\n", sep = "")
          mout$par       <- fit_bfgs$par
          mout$objective <- fit_bfgs$value
          cur_grad       <- .max_grad(mout$par)
          .checkpoint(mout, paste("sandwich", restart, "(reverted)"))
          break
        }
        .checkpoint(mout, paste("sandwich", restart))

        cat("  Restart", restart, "complete: max|grad| =",
            round(cur_grad, 6), "\n")

        ## Exit ##
        # 1. Gradient low enough for Newton — clean exit
        if (cur_grad < newton_grad_thresh) {
          cat("  Gradient <", newton_grad_thresh,
              "\u2014 sandwich converged, proceeding to Newton.\n")
          sandwich_done <- TRUE
          break
        }

        # 2. Gradient stall — two consecutive cycles with no meaningful gradient reduction
        grad_improvement <- prev_grad - cur_grad
        if (!is.finite(grad_improvement) || grad_improvement < newton_grad_thresh * 0.01) {
          no_improve <- no_improve + 1L
          cat("  Gradient not improving (delta =",
              round(grad_improvement, 6), ") — strike", no_improve, "of 2\n")
          if (no_improve >= 2L) {
            cat("  Sandwich stalled on gradient — exiting after",
                restart, "cycles.\n")
            break
          }
        } else {
          no_improve <- 0L
        }
        prev_grad <- cur_grad
      }

      ## Newton polishing ##
      if (newtonSteps > 0) {

        cur_grad <- .max_grad(mout$par)

        if (!is.finite(cur_grad) || cur_grad >= newton_grad_thresh) {
          cat("\n  WARNING: gradient =", round(cur_grad, 4),
              ">= newton_grad_thresh (", newton_grad_thresh, ")",
              "— skipping Newton (Hessian likely ill-conditioned).\n")
        } else {
          cat("\n--- Newton polishing steps ---\n")
          CurrentStage <<- "Newton polish"
          newton_par   <- model$env$last.par.best
          if (!is.finite(.max_grad(newton_par))) newton_par <- mout$par

          for (ns in seq_len(newtonSteps)) {
            improved <- FALSE
            tryCatch({
              suppress_print <- TRUE                     # silence fn prints during Hessian
              H    <- optimHess(newton_par, model$fn, model$gr)
              suppress_print <- FALSE
              g    <- as.vector(model$gr_Orig(newton_par))
              step <- solve(H, g)

              base_obj   <- model$fn(newton_par)
              step_size  <- 1.0
              for (backstep in seq_len(10)) {
                candidate <- newton_par - step_size * step
                if (has_bounds) {
                  candidate <- pmax(candidate, RunSpecs$lowBnd)
                  candidate <- pmin(candidate, RunSpecs$uppBnd)
                }
                if (model$fn(candidate) < base_obj && is.finite(.max_grad(candidate))) {
                  newton_par <- candidate
                  improved   <- TRUE
                  break
                }
                step_size <- step_size / 2
              }
              if (!improved) {
                if (max(abs(g)) < newton_grad_thresh) {
                  cat("  Newton step", ns,
                      "\u2014 already at optimum, no further improvement needed.\n")
                } else {
                  cat("  Newton step", ns,
                      "\u2014 line search failed, stopping Newton.\n")
                }
              } else {
                ng <- .max_grad(newton_par)
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
              if (exists("H", inherits = TRUE) && exists("g", inherits = TRUE)) {
                .diagnose_singular_hessian(H, g, pnames, e = e)
              } else {
                cat("    (H not yet computed — failure occurred before solve())\n")
              }
            })
            suppress_print <- FALSE
            if (!improved) break
          }

          # Final nlminb from Newton-polished parameters
          cat("  Final nlminb (post-Newton)\n")
          start_obj <- model$fn(newton_par)
          .start_stage("Post-Newton nlminb", start_obj)
          initBestFn <- start_obj

          ctrl <- list(iter.max = MaXeVaL, eval.max = MaXeVaL,
                       rel.tol = 1e-12, x.tol = 1e-12, abs.tol = 0)
          mout_newton <- nlminb(newton_par, model$fn, model$gr,
                                lower = RunSpecs$lowBnd, upper = RunSpecs$uppBnd,
                                control = ctrl)
          if (is.finite(.max_grad(mout_newton$par))) {
            mout <- mout_newton
          } else {
            cat("  WARNING: post-Newton nlminb ended on a non-finite gradient — keeping Newton point\n")
            mout$par       <- newton_par
            mout$objective <- model$fn(newton_par)
          }
          .report_fit(mout, model, pnames, initBestFn,
                      label = "  Final nlminb (post-Newton)")
          TotalEval <<- TotalEval + FnCallNo
          .checkpoint(mout, "post-Newton nlminb")
        }
      }

      ## Bound diagnostics ##
      .check_bounds(mout$par, RunSpecs$lowBnd, RunSpecs$uppBnd, pnames, model)
    }

    if (nonfinite_gr > 0)
      cat("  NOTE: Phase", CurrPhase, "had", nonfinite_gr,
          "gradient evaluation(s) with non-finite values (returned as zeros to the optimiser).\n")

    ## Store and save parameters ##
    pars        <- mout$par
    names(pars) <- pnames
    ParOld      <- mout$par

    pout <- unlist(parameters)
    pout[names(pout) %in% names(pars)] <- pars
    suffix <- ifelse(is_final, " final", CurrPhase)
    write.table(pout, paste0("Output/model", suffix, ".par"),
                sep = "\t", col.names = c("name\test"), quote = FALSE)

    if (is_final) {
      assign("ProfileReport", model$report(),               envir = .GlobalEnv)
      assign("ProfileGrad",   abs(model$gr_Orig(mout$par)), envir = .GlobalEnv)
    }

    ## Report (final phase only) ##
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
