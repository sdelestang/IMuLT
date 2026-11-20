#' Global goodness-of-fit diagnostics across all data types
#'
#' Scores the fit of every CPUE series, length-composition series (Francis
#' mean lengths), larval series and tag release-recapture series in an IMuLT
#' report, using the same set of tests for each, and flags series that fail
#' against absolute (not ranked) criteria, so that "nothing flagged" is a
#' possible outcome. Results go to a "Diagnostics" tab via
#' \code{addplot}/\code{addtable}.
#'
#' Tests applied to each series:
#' \describe{
#'   \item{CCC}{Lin's concordance correlation (obs vs pred; log scale for CPUE
#'     and larval, mean length for comps, sqrt counts for tags), split into
#'     precision (Pearson r) and accuracy (C_b). Flagged when below the 5th
#'     percentile of CCC values from \code{nsim} data sets simulated from the
#'     predictions with the stated error, i.e. worse than a correct model would
#'     typically give for that series' contrast and noise.}
#'   \item{SDNR}{SD of normalised residuals, using the error the likelihood
#'     assumes (CPUE: input CV x self-tuned Sigma_used; comps: input N; larval:
#'     reported SD; tags: overall quasi-Poisson dispersion). Flagged when above
#'     sqrt(qchisq(0.95, n-1)/(n-1)) (Francis 2011).}
#'   \item{Runs}{One-sided Wald-Wolfowitz runs test on residual signs (too few
#'     runs = systematic misfit). Flagged when p < 0.05 (Carvalho et al. 2021).}
#'   \item{Trend}{Slope of residuals over the last \code{trend_n} points.
#'     Flagged when p < 0.05.}
#' }
#' Status: green = no flags; amber = one flag; red = two or more flags.
#' \code{Conflict} (reported, not used in the status) marks a failed runs test
#' whose mean residual has the same sign as a flagged series of a different
#' data type in the same area.
#'
#' @param dat Character vector of the report lines (as passed to \code{findNclean}).
#' @param rundir Run directory for outputs.
#' @param fleets Fleet table with columns \code{fleet} and \code{newarea}.
#' @param areanames Area labels in model area order (used for larval and tag
#'   data, which are reported by area number). Should match \code{fleets$newarea}.
#' @param lbins Length-bin mid-points for the composition data.
#' @param nsim Number of simulated data sets for the calibrated CCC threshold.
#' @param trend_n Number of most recent points used for the residual-trend test.
#' @param min_n Minimum number of points for a series to be tested.
#' @param seed Random seed for the simulations.
#'
#' @return Invisibly, the data frame of diagnostics (one row per series).
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_tile scale_colour_manual scale_fill_manual facet_grid labs theme_bw theme element_text
#' @importFrom stats cor cov var sd lm qchisq pnorm quantile rnorm rpois rnbinom aggregate
#' @importFrom utils tail
#' @export
GlobalDiagnostics <- function(dat, rundir, fleets,
                              areanames = c("C_shallow","C_deep","B_shallow","B_deep",
                                            "A_shallow","A_deep","BN_shallow","BN_deep"),
                              lbins = seq(44.5, by = 3, length.out = 31),
                              nsim = 500, trend_n = 10, min_n = 6, seed = 1) {

  print("Making global diagnostics")
  set.seed(seed)

  ## ---------- helpers ----------
  # Read a report block: lines after the header (and its column-name line) until next '#' or blank
  rl_block <- function(header, skip = 2) {
    i <- grep(paste0("^", header), dat)[1]
    if (is.na(i)) return(NULL)
    j <- i + skip
    out <- character(0)
    while (j <= length(dat) && nzchar(trimws(dat[j])) && !startsWith(dat[j], "#")) {
      out <- c(out, dat[j]); j <- j + 1 }
    if (!length(out)) return(NULL)
    strsplit(trimws(out), "\\s+")
  }
  ccc_parts <- function(x, y) {
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
    if (length(x) < 3 || var(x) == 0 || var(y) == 0) return(c(ccc = NA, r = NA, cb = NA))
    r   <- cor(x, y)
    ccc <- 2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
    c(ccc = ccc, r = r, cb = ccc / r)
  }
  runs_p <- function(res) {                      # one-sided: too few runs
    s <- sign(res[is.finite(res) & res != 0])
    n1 <- sum(s > 0); n2 <- sum(s < 0); n <- n1 + n2
    if (n1 == 0 || n2 == 0 || n < min_n) return(NA)
    runs <- 1 + sum(diff(s) != 0)
    mu <- 2 * n1 * n2 / n + 1
    v  <- 2 * n1 * n2 * (2 * n1 * n2 - n) / (n^2 * (n - 1))
    pnorm((runs - mu) / sqrt(v))
  }
  trend_p <- function(t, res) {
    ok <- is.finite(res); t <- t[ok]; res <- res[ok]
    if (length(res) < min_n) return(NA)
    k <- tail(order(t), trend_n)
    if (length(k) < 5) return(NA)
    summary(lm(res[k] ~ t[k]))$coefficients[2, 4]
  }
  # one series -> one row. sim_fun() returns one simulated obs vector on the same scale as x
  score <- function(type, id, area, t, x, y, res, sim_fun) {
    n  <- sum(is.finite(res))
    cp <- ccc_parts(x, y)
    if (n >= min_n && is.finite(cp["ccc"])) {
      sims <- replicate(nsim, ccc_parts(sim_fun(), y)["ccc"])
      ccc_lim <- quantile(sims, 0.05, na.rm = TRUE)
    } else ccc_lim <- NA
    sdnr <- if (n >= min_n) sd(res, na.rm = TRUE) else NA
    sdnr_lim <- if (n >= min_n) sqrt(qchisq(0.95, n - 1) / (n - 1)) else NA
    data.frame(Type = type, Series = id, Area = area, N = n,
               CCC = cp["ccc"], Precision_r = cp["r"], Accuracy_Cb = cp["cb"],
               CCC_limit = ccc_lim, SDNR = sdnr, SDNR_limit = sdnr_lim,
               Runs_p = runs_p(res[order(t)]), Trend_p = trend_p(t, res),
               MeanRes = mean(res, na.rm = TRUE), row.names = NULL)
  }
  fleet_area <- function(f) as.character(fleets$newarea[match(f, fleets$fleet)])

  out <- list()

  ## ---------- CPUE ----------
  b <- rl_block("#Index data")
  if (!is.null(b)) {
    cp <- as.data.frame(do.call(rbind, b), stringsAsFactors = FALSE)
    names(cp) <- c("Data_set","Fleet","Sex","Year","Time_step","Observed","Relative_CV","Predicted","Residual")
    cp[] <- lapply(cp, as.numeric)
    tune <- rl_block("#Data_set N Sigma", skip = 1)
    sig  <- if (!is.null(tune)) setNames(as.numeric(sapply(tune, `[`, 4)), as.numeric(sapply(tune, `[`, 1))) else NULL
    for (ds in sort(unique(cp$Data_set))) {
      d  <- cp[cp$Data_set == ds & cp$Observed > 0 & cp$Predicted > 0, ]
      if (nrow(d) < 3) next
      # Error the likelihood actually assumes: input CV x self-tuned Sigma_used (report Residual = log(O/P)/CV)
      smult <- if (!is.null(sig) && is.finite(sig[as.character(ds)])) sig[as.character(ds)] else 1
      sd_i <- d$Relative_CV * smult
      lo <- log(d$Observed); lp <- log(d$Predicted)
      id <- paste0("DS", ds, " F", d$Fleet[1],
                   ifelse(d$Sex[1] == 0, " F", ifelse(d$Sex[1] == 1, " M", "")),
                   " ts", paste(unique(d$Time_step), collapse = ""))
      out[[length(out) + 1]] <- score("CPUE", id, fleet_area(d$Fleet[1]),
                                      d$Year + (d$Time_step - 1) / 4, lo, lp, (lo - lp) / sd_i,
                                      function() lp + rnorm(length(lp), 0, sd_i))
    }
  }

  ## ---------- Length compositions (Francis mean length) ----------
  b <- rl_block("#Obs/Pred Fleet Sex Year Step", skip = 1)
  if (!is.null(b)) {
    op  <- sapply(b, `[`, 1)
    hdr <- t(sapply(b, function(z) as.numeric(z[2:7])))   # Fleet Sex Year Step Nsamp EffN
    prp <- lapply(b, function(z) as.numeric(z[-(1:7)]))
    nb  <- min(length(lbins), min(lengths(prp)))
    L <- lbins[1:nb]
    iO <- which(op == "O"); iP <- which(op == "P")
    key <- function(i) apply(hdr[i, 1:4, drop = FALSE], 1, paste, collapse = "_")
    iP <- iP[match(key(iO), key(iP))]
    ml <- data.frame(Fleet = hdr[iO, 1], Sex = hdr[iO, 2], Year = hdr[iO, 3], Step = hdr[iO, 4],
                     N = hdr[iO, 5],
                     Obar = sapply(iO, function(i) sum(prp[[i]][1:nb] * L) / sum(prp[[i]][1:nb])),
                     Pbar = sapply(iP, function(i) sum(prp[[i]][1:nb] * L) / sum(prp[[i]][1:nb])),
                     Pvar = sapply(iP, function(i) { p <- prp[[i]][1:nb] / sum(prp[[i]][1:nb])
                                                     sum(p * L^2) - sum(p * L)^2 }))
    ml$se <- sqrt(ml$Pvar / pmax(ml$N, 1e-6))        # Francis: se of mean length from input N
    for (g in unique(paste(ml$Fleet, ml$Sex))) {
      d <- ml[paste(ml$Fleet, ml$Sex) == g & is.finite(ml$Obar) & is.finite(ml$Pbar), ]
      if (nrow(d) < 3) next
      id <- paste0("F", d$Fleet[1], ifelse(d$Sex[1] == 1, " F", " M"))
      out[[length(out) + 1]] <- score("LengthComp", id, fleet_area(d$Fleet[1]),
                                      d$Year + (d$Step - 1) / 4, d$Obar, d$Pbar, (d$Obar - d$Pbar) / d$se,
                                      function() d$Pbar + rnorm(nrow(d), 0, d$se))
    }
  }

  ## ---------- Larval ----------
  b <- rl_block("#Larval data")
  if (!is.null(b)) {
    lv <- as.data.frame(do.call(rbind, b), stringsAsFactors = FALSE)
    names(lv) <- c("Area","Year","Observed","SD","Predicted","Residual")
    lv[] <- lapply(lv, as.numeric)
    for (a in sort(unique(lv$Area))) {
      d <- lv[lv$Area == a & lv$Observed > 0 & lv$Predicted > 0, ]
      if (nrow(d) < 3) next
      sdl <- sqrt(log(1 + (d$SD / d$Observed)^2))
      lo <- log(d$Observed); lp <- log(d$Predicted)
      out[[length(out) + 1]] <- score("Larval", paste("Larval area", a), areanames[a],
                                      d$Year, lo, lp, (lo - lp) / sdl,
                                      function() lp + rnorm(length(lp), 0, sdl))
    }
  }

  ## ---------- Tags (recaptures by release x recapture area, per year) ----------
  b <- rl_block("#Tagging data - Tag numbers", skip = 2)
  if (!is.null(b)) {
    tg <- as.data.frame(do.call(rbind, b), stringsAsFactors = FALSE)
    names(tg) <- c("Sex","RelArea","RecArea","Year","Tstep","RepSplit","Obs","Est")
    tg[] <- lapply(tg, as.numeric)
    tg <- aggregate(cbind(Obs, Est) ~ RelArea + RecArea + Year, data = tg, FUN = sum)
    # Tag counts are overdispersed relative to Poisson; scale by the overall dispersion so series
    # are judged against typical tag misfit (quasi-Poisson), and simulate negative binomial with it
    phi <- with(tg[tg$Est > 0, ], max(1, mean((Obs - Est)^2 / Est)))
    rtag <- function(E) if (phi > 1.0001) rnbinom(length(E), size = E / (phi - 1), mu = E) else rpois(length(E), E)
    for (g in unique(paste(tg$RelArea, tg$RecArea))) {
      d <- tg[paste(tg$RelArea, tg$RecArea) == g & tg$Est > 0, ]
      if (nrow(d) < 3 || sum(d$Obs) == 0) next
      id <- paste0("Tag ", areanames[d$RelArea[1]], " -> ", areanames[d$RecArea[1]])
      out[[length(out) + 1]] <- score("Tags", id, areanames[d$RelArea[1]],
                                      d$Year, sqrt(d$Obs), sqrt(d$Est), (d$Obs - d$Est) / sqrt(phi * d$Est),
                                      function() sqrt(rtag(d$Est)))
    }
  }

  if (!length(out)) { print("No data found for global diagnostics"); return(invisible(NULL)) }
  diag <- do.call(rbind, out)

  ## ---------- Flags and status ----------
  diag$Flag_CCC   <- diag$CCC  < diag$CCC_limit
  diag$Flag_SDNR  <- diag$SDNR > diag$SDNR_limit
  diag$Flag_Runs  <- diag$Runs_p  < 0.05
  diag$Flag_Trend <- diag$Trend_p < 0.05
  fl <- c("Flag_CCC","Flag_SDNR","Flag_Runs","Flag_Trend")
  diag$Nflags <- rowSums(diag[, fl] == TRUE, na.rm = TRUE)
  # cross-data conflict: failed runs test, and a flagged series of another type in the same area with same residual sign
  diag$Conflict <- sapply(seq_len(nrow(diag)), function(i) {
    if (!isTRUE(diag$Flag_Runs[i]) || is.na(diag$Area[i])) return(FALSE)
    o <- diag$Area == diag$Area[i] & diag$Type != diag$Type[i] & diag$Nflags > 0 &
         sign(diag$MeanRes) == sign(diag$MeanRes[i])
    any(o, na.rm = TRUE) })
  diag$Status <- ifelse(diag$Nflags >= 2, "red",
                        ifelse(diag$Nflags == 1, "amber", "green"))
  diag <- diag[order(-match(diag$Status, c("green","amber","red")), -diag$Nflags, diag$CCC), ]
  diag$Type <- factor(diag$Type, levels = c("CPUE","LengthComp","Larval","Tags"))
  stcol <- c(green = "forestgreen", amber = "orange", red = "red")

  ## ---------- Table ----------
  tab <- diag
  numc <- sapply(tab, is.numeric); tab[numc] <- lapply(tab[numc], round, 3)
  row.names(tab) <- NULL
  addtable(intable = tab, filen = "Global_Diagnostics.csv", rundir = rundir, category = "Diagnostics",
           caption = paste("Goodness-of-fit diagnostics for every CPUE, length-composition (Francis mean length),",
                           "larval and tag series. Flags use absolute criteria: CCC below the 5th percentile of",
                           "CCC from data simulated from the fit; SDNR above the Francis (2011) chi-square limit;",
                           "one-sided runs test p < 0.05; residual trend over the last", trend_n, "points p < 0.05.",
                           "Red = 2+ flags; amber = 1 flag. Conflict = failed runs test with same-sign misfit in another data type in the same area."))

  ## ---------- Plot 1: rug of CCC by data type ----------
  nflag <- sum(diag$Status != "green")
  filename <- filenametopath(rundir, "Diagnostics CCC rug.png")
  plotprep(width = 9, height = 5, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1, 1))
  p <- ggplot(diag, aes(x = CCC, y = Type, colour = Status)) +
    geom_point(shape = "|", size = 9) +
    geom_text(data = diag[diag$Status != "green", ], aes(label = Series),
              angle = 40, hjust = 0, vjust = 0, nudge_y = 0.12, size = 2.7, show.legend = FALSE) +
    scale_colour_manual(values = stcol, drop = FALSE) +
    labs(x = "Concordance correlation (obs vs pred)", y = NULL,
         title = if (nflag == 0) "No series flagged" else paste(nflag, "of", nrow(diag), "series flagged")) +
    theme_bw() + theme(axis.text = element_text(size = 11), legend.position = "bottom")
  print(p)
  addplot(filen = filename, rundir = rundir, category = "Diagnostics",
          caption = paste("CCC of every series by data type, coloured by status. Only flagged (amber/red)",
                          "series are labelled; a low CCC is not flagged if it is within what the stated error allows."))

  ## ---------- Plot 2: heat map of flagged series ----------
  fd <- diag[diag$Status != "green", ]
  if (nrow(fd) > 0) {
    hm <- do.call(rbind, lapply(fl, function(f)
      data.frame(Series = paste0(fd$Series, " (", fd$Type, ")"), Test = sub("Flag_", "", f),
                 Result = ifelse(is.na(fd[[f]]), "n/a", ifelse(fd[[f]], "fail", "pass")))))
    hm$Series <- factor(hm$Series, levels = rev(unique(paste0(fd$Series, " (", fd$Type, ")"))))
    hm$Test <- factor(hm$Test, levels = c("CCC","SDNR","Runs","Trend"))
    filename <- filenametopath(rundir, "Diagnostics heat map.png")
    plotprep(width = 7, height = min(2 + 0.22 * nrow(fd), 11), filename = filename, cex = 0.9, verbose = FALSE)
    parset(plots = c(1, 1))
    p <- ggplot(hm, aes(x = Test, y = Series, fill = Result)) +
      geom_tile(colour = "white") +
      scale_fill_manual(values = c(fail = "red", pass = "grey85", `n/a` = "white")) +
      labs(x = NULL, y = NULL) + theme_bw() + theme(axis.text.y = element_text(size = 8))
    print(p)
    addplot(filen = filename, rundir = rundir, category = "Diagnostics",
            caption = "Tests failed by each flagged series (series with no flags are in the table only).")
  }

  ## ---------- Plot 3: flags by area and data type ----------
  ad <- aggregate(cbind(Flagged = Status != "green", Total = 1) ~ Area + Type,
                  data = diag[!is.na(diag$Area), ], FUN = sum)
  ad$Label <- paste0(ad$Flagged, "/", ad$Total)
  filename <- filenametopath(rundir, "Diagnostics by area.png")
  plotprep(width = 7, height = 4.5, filename = filename, cex = 0.9, verbose = FALSE)
  parset(plots = c(1, 1))
  p <- ggplot(ad, aes(x = Type, y = Area, fill = Flagged / Total)) +
    geom_tile(colour = "white") + geom_text(aes(label = Label), size = 4) +
    ggplot2::scale_fill_gradient(low = "grey95", high = "red", limits = c(0, 1), name = "Prop.\nflagged") +
    labs(x = NULL, y = NULL) + theme_bw() + theme(axis.text = element_text(size = 11))
  print(p)
  addplot(filen = filename, rundir = rundir, category = "Diagnostics",
          caption = paste("Flagged/total series by area and data type. An area flagged across several data types",
                          "points to a cause in that area (e.g. movement, M, growth) rather than one data set's weighting."))

  invisible(diag)
}
