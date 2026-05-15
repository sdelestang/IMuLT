#' Compare Legal Biomass Across IMuLT Model Runs
#'
#' Scans a summary folder for model run subfolders, extracts Legal Biomass
#' by area from each Output.RL file, and produces comparison plots.
#'
#' @param summary_dir Character. Path to the Summary folder containing
#'   model run subfolders. Default \code{"Output/Summary"}.
#' @param by_area Logical. If \code{TRUE}, plot separate panels per area.
#'   If \code{FALSE}, sum across areas and plot total biomass. Default \code{FALSE}.
#' @param relative Logical. If \code{TRUE}, plot biomass relative to virgin
#'   (B/B0). If \code{FALSE}, plot absolute biomass. Default \code{FALSE}.
#' @param runs Character vector or \code{NULL}. Specific subfolder names to
#'   include. If \code{NULL}, all subfolders containing Output.RL are used.
#' @param scale Numeric. Divisor for absolute biomass (e.g. 1000 for tonnes).
#'   Ignored when \code{relative = TRUE}. Default 1.
#' @param ylab Character or \code{NULL}. Y-axis label. If \code{NULL}, chosen
#'   automatically based on \code{relative}.
#'
#' @return A data frame (invisibly) with columns: Year, Area, Biomass, Virgin,
#'   Relative, Run.
#'
#' @examples
#' \dontrun{
#' compare_legal_biomass("Output/Summary")
#' compare_legal_biomass("Output/Summary", relative = TRUE)
#' compare_legal_biomass("Output/Summary", by_area = TRUE)
#' compare_legal_biomass("Output/Summary", by_area = TRUE, relative = TRUE)
#' }
#'
#' @export
compare_legal_biomass <- function(Vers=NULL, summary_dir = "Output/Summary",
                                  by_area = FALSE,
                                  relative = FALSE,
                                  runs = NULL,
                                  scale = 1,
                                  ylab = NULL) {

  if (is.null(ylab)) ylab <- if (relative) "B/B0" else "Legal Biomass (t)"

  # ── Locate run folders ──────────────────────────────────────
  if (!dir.exists(summary_dir))
    stop("Summary directory not found: ", summary_dir, call. = FALSE)

  all_dirs <- list.dirs(summary_dir, full.names = TRUE, recursive = FALSE)

  if(!is.null(Vers)) {
    VersLong <- paste0("Output/Summary/",Vers)
    Matches <- match(tolower(VersLong), tolower(all_dirs)); Matches <- Matches[!is.na(Matches)]
    if (length(Matches) == 0){
      stop("None of the specified runs found in ", summary_dir, call. = FALSE)}
    if (length(Matches) != length(Vers)){
      warning(paste0("Only ",length(Matches), " of your ",length(Vers)," specified runs were found in ", summary_dir,'\nEntered: ',paste(Vers[Matches], collapse=' '),'\nAvailable: ',paste(gsub("Output/Summary/", "",all_dirs), collapse=' ')))}
    all_dirs <- all_dirs[Matches]
      }

  if (!is.null(runs)) {
    all_dirs <- all_dirs[basename(all_dirs) %in% runs]
    if (length(all_dirs) == 0)
      stop("None of the specified runs found in ", summary_dir, call. = FALSE)
  }

  # Find Output.RL in each folder (may be nested one level deeper)
  rl_paths <- vapply(all_dirs, function(d) {
    found <- list.files(d, pattern = "^Output\\.RL$",
                        recursive = TRUE, full.names = TRUE)
    if (length(found) > 0) found[1] else NA_character_
  }, character(1), USE.NAMES = FALSE)

  keep     <- !is.na(rl_paths)
  all_dirs <- all_dirs[keep]
  rl_paths <- rl_paths[keep]

  if (length(all_dirs) == 0)
    stop("No Output.RL files found in subfolders of ", summary_dir, call. = FALSE)

  # Full folder names and short labels for plotting
  full_names <- basename(all_dirs)
  run_labels <- vapply(full_names, function(nm)
    substr(nm, max(1, nchar(nm) - 20), nchar(nm)),
    character(1), USE.NAMES = FALSE)
  if (anyDuplicated(run_labels)) run_labels <- full_names

  cat("Found", length(all_dirs), "model runs:",
      paste(run_labels, collapse = ", "), "\n")

  # ── Extract Legal Biomass from each run ─────────────────────
  all_data <- list()

  for (i in seq_along(all_dirs)) {
    rl_path <- rl_paths[i]
    dat <- read.table(rl_path, comment.char = "?", fill = TRUE,
                      blank.lines.skip = FALSE, stringsAsFactors = FALSE,
                      col.names = paste0("V", 1:200))

    # Legal Biomass by area: returns long format (area, Year, est, se)
    lb <- findNclean(c("#Legal", "Biomass", "by"), dat, 2, convert = 2)

    if (identical(lb, NA) || is.null(lb)) {
      warning("Could not extract Legal Biomass from: ", full_names[i],
              call. = FALSE)
      next
    }

    # Identify columns by name (area, Year, est, se)
    cnames <- tolower(colnames(lb))
    col_area <- which(cnames == "area")[1]
    col_year <- which(cnames == "year")[1]
    col_est  <- which(cnames == "est")[1]

    if (is.na(col_area) || is.na(col_year) || is.na(col_est)) {
      # Fallback: assume columns are area=1, Year=2, est=3
      col_area <- 1; col_year <- 2; col_est <- 3
    }

    areas   <- as.numeric(lb[, col_area])
    years   <- as.numeric(lb[, col_year])
    biomass <- as.numeric(lb[, col_est])

    # Extract virgin legal biomass (one value per area)
    virgin_raw <- findNclean(c("#Virgin", "Legal"), dat, 1)

    if (identical(virgin_raw, NA) || is.null(virgin_raw)) {
      if (relative) {
        warning("Could not extract Virgin Legal Biomass from: ", full_names[i],
                " - skipping for relative plot", call. = FALSE)
        next
      }
      virgin_lookup <- NULL
    } else {
      # virgin_raw: vector or data.frame, one value per area
      if (is.data.frame(virgin_raw)) {
        virgin_vals <- as.numeric(virgin_raw[, 1])
      } else {
        virgin_vals <- as.numeric(virgin_raw)
      }
      # Build lookup: area index (1, 2, ...) → virgin value
      virgin_lookup <- setNames(virgin_vals, seq_along(virgin_vals))
    }

    for (j in seq_len(nrow(lb))) {
      a   <- areas[j]
      vir <- if (!is.null(virgin_lookup) && as.character(a) %in% names(virgin_lookup)) {
        virgin_lookup[as.character(a)]
      } else {
        NA_real_
      }

      all_data[[length(all_data) + 1]] <- data.frame(
        Year     = years[j],
        Area     = a,
        Biomass  = biomass[j] / scale,
        Virgin   = as.numeric(vir) / scale,
        Relative = biomass[j] / as.numeric(vir),
        Run      = run_labels[i],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(all_data) == 0)
    stop("No Legal Biomass data extracted from any run.", call. = FALSE)

  result <- do.call(rbind, all_data)
  rownames(result) <- NULL

  # ── Plot ────────────────────────────────────────────────────
  if (by_area) {
    .plot_by_area(result, ylab, relative)
  } else {
    .plot_total(result, ylab, relative)
  }

  invisible(result)
}


# ── Internal: plot total biomass (summed across areas) ────────

.plot_total <- function(df, ylab, relative) {

  runs   <- unique(df$Run)
  n_runs <- length(runs)
  cols   <- .run_colours(n_runs)

  if (relative) {
    # Sum biomass and virgin across areas per year, then ratio
    agg <- aggregate(cbind(Biomass, Virgin) ~ Year + Run, data = df, FUN = sum)
    agg$PlotVal <- agg$Biomass / agg$Virgin
    main_title  <- "Relative Legal Biomass (B/B0) \u2014 Model Comparison"
  } else {
    agg <- aggregate(Biomass ~ Year + Run, data = df, FUN = sum)
    agg$PlotVal <- agg$Biomass
    main_title  <- "Total Legal Biomass \u2014 Model Comparison"
  }

  ylim <- c(0, max(agg$PlotVal, na.rm = TRUE) * 1.1)
  xlim <- range(agg$Year)

  par(mfrow = c(1, 1), mar = c(4, 5, 3, 1))
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = "Year", ylab = ylab,
       main = main_title, las = 1)

  if (relative) {
    abline(h = 1.0, col = "grey70", lty = 2)
    abline(h = 0.4, col = "green3", lwd = 1.5, lty = 2)
    abline(h = 0.35, col = "orange", lwd = 1.5, lty = 2)
    abline(h = 0.2, col = "red", lwd = 1.5, lty = 2)
  }

  for (r in seq_along(runs)) {
    sub <- agg[agg$Run == runs[r], ]
    sub <- sub[order(sub$Year), ]
    lines(sub$Year, sub$PlotVal, col = cols[r], lwd = 2.5)
    points(sub$Year, sub$PlotVal, col = cols[r], pch = 16, cex = 0.6)
  }

  legend("topright", legend = runs, col = cols, lwd = 2.5,
         bty = "n", cex = 0.85, ncol = ceiling(n_runs / 5))
}


# ── Internal: plot by area (faceted) ──────────────────────────

.plot_by_area <- function(df, ylab, relative) {

  areas   <- sort(unique(df$Area))
  n_areas <- length(areas)
  runs    <- unique(df$Run)
  n_runs  <- length(runs)
  cols    <- .run_colours(n_runs)

  ncol_p <- min(n_areas, 3)
  nrow_p <- ceiling(n_areas / ncol_p)
  par(mfrow = c(nrow_p, ncol_p), mar = c(4, 5, 3, 1), oma = c(2, 0, 2, 0))

  xlim <- range(df$Year)

  for (a in areas) {
    sub_a <- df[df$Area == a, ]
    sub_a$PlotVal <- if (relative) sub_a$Relative else sub_a$Biomass

    ylim <- c(0, max(sub_a$PlotVal, na.rm = TRUE) * 1.1)

    plot(NULL, xlim = xlim, ylim = ylim,
         xlab = "Year", ylab = ylab,
         main = paste("Area", a), las = 1)

    if (relative) {
      abline(h = 0.4, col = "green3", lwd = 1, lty = 2)
      abline(h = 0.35, col = "orange", lwd = 1, lty = 2)
      abline(h = 0.2, col = "red", lwd = 1, lty = 2)
    }

    for (r in seq_along(runs)) {
      sub <- sub_a[sub_a$Run == runs[r], ]
      sub <- sub[order(sub$Year), ]
      lines(sub$Year, sub$PlotVal, col = cols[r], lwd = 2)
      points(sub$Year, sub$PlotVal, col = cols[r], pch = 16, cex = 0.5)
    }
  }

  # Shared legend in outer margin
  main_title <- if (relative) {
    "Relative Legal Biomass (B/B0) by Area"
  } else {
    "Legal Biomass by Area \u2014 Model Comparison"
  }
  mtext(main_title, outer = TRUE, cex = 1.1, font = 2)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
      mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", legend = runs, col = cols, lwd = 2.5,
         horiz = FALSE, ncol = min(n_runs, 5), bty = "n", cex = 0.85)
}


# ── Internal: colour palette ─────────────────────────────────

.run_colours <- function(n) {
  base <- c("steelblue", "tomato", "seagreen", "darkorange",
            "mediumpurple", "goldenrod", "deeppink", "cyan4")
  if (n <= length(base)) return(base[1:n])
  grDevices::colorRampPalette(base)(n)
}


