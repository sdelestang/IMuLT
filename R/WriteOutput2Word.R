#' Write IMuLT Output to Word Document
#'
#' Reads resultTable.csv from an IMuLT model run and assembles all plots
#' and tables into a formatted Word document, mirroring the HTML tab structure.
#' When called without arguments from the model working directory, it auto-detects
#' the most recent run in Output/Summary.
#'
#' @param summary_dir Path to the summary directory containing run folders.
#'   Defaults to "Output/Summary".
#' @param run Character string specifying which run folder to use (e.g. "run2").
#'   If NULL (default), uses the most recent run based on resultTable.csv timestamp.
#' @param output_file Name of the output Word file. Defaults to "IMuLT_Report.docx"
#'   saved in the run directory.
#' @param img_width Width of images in inches (default 6.5 for portrait A4 with margins).
#' @param img_height Height of images in inches (default 5). Set to NULL for auto-scaling
#'   based on the actual image aspect ratio.
#' @param title Optional title for the front page. If NULL, derives from the model
#'   directory name.
#'
#' @return Invisible path to the output file.
#'
#' @examples
#' \dontrun{
#' # From the model working directory - auto-detects latest run
#' WriteWord()
#'
#' # Specify a particular run
#' WriteWord(run = "run2")
#'
#' # Auto aspect ratio
#' WriteWord(img_height = NULL)
#'
#' # Custom summary directory
#' WriteWord(summary_dir = "C:/path/to/Output/Summary", run = "run3")
#' }
#'
#' @name WriteWord
#' @export
WriteWord <- function(summary_dir = "Output/Summary",
                      run         = NULL,
                      output_file = NULL,
                      img_width   = 6.5,
                      img_height  = 5,
                      title       = NULL) {

  # --- Check dependencies ---
  if (!requireNamespace("officer", quietly = TRUE))
    stop("Package 'officer' is required. Install with: install.packages('officer')")
  if (!requireNamespace("flextable", quietly = TRUE))
    stop("Package 'flextable' is required. Install with: install.packages('flextable')")

  library(officer)
  library(flextable)

  # --- Locate run directory ---
  if (!dir.exists(summary_dir))
    stop("Summary directory not found: ", summary_dir, call. = FALSE)

  all_dirs <- list.dirs(summary_dir, full.names = TRUE, recursive = FALSE)

  # Filter to directories that contain a resultTable.csv
  has_rt <- vapply(all_dirs, function(d)
    file.exists(file.path(d, "resultTable.csv")), logical(1))
  all_dirs <- all_dirs[has_rt]

  if (length(all_dirs) == 0)
    stop("No resultTable.csv found in any subfolder of: ", summary_dir, call. = FALSE)

  cat("Available runs:", paste(basename(all_dirs), collapse = ", "), "\n")

  if (!is.null(run)) {
    # User specified a run
    run_dir <- all_dirs[basename(all_dirs) == run]
    if (length(run_dir) == 0)
      stop("Run '", run, "' not found in ", summary_dir,
           "\n  Available: ", paste(basename(all_dirs), collapse = ", "),
           call. = FALSE)
    run_dir <- run_dir[1]
  } else {
    # Auto-detect most recent run based on resultTable.csv modification time
    rt_times <- file.mtime(file.path(all_dirs, "resultTable.csv"))
    run_dir  <- all_dirs[which.max(rt_times)]
  }

  cat("Using run:", basename(run_dir), "\n")

  # --- Read and validate resultTable ---
  rt_path <- file.path(run_dir, "resultTable.csv")
  res <- read.csv(rt_path, stringsAsFactors = FALSE, strip.white = TRUE)

  # Clean trailing empty columns (the CSV has a trailing comma)
  res <- res[, !grepl("^X", names(res)) & names(res) != ""]

  # Extract just the filename from the full path
  res$basename <- basename(res$file)

  # --- Set up output path ---
  if (is.null(output_file)) {
    output_file <- file.path(run_dir, "IMuLT_Report.docx")
  } else if (!grepl("[/\\\\]", output_file)) {
    output_file <- file.path(run_dir, output_file)
  }

  # --- Derive title ---
  if (is.null(title)) {
    # Walk up to the model directory name (e.g. "7Area1AgeRun96_25")
    title <- basename(dirname(dirname(run_dir)))
  }

  # --- Helper: get image dimensions for auto aspect ratio ---
  get_img_dims <- function(img_path, target_width) {
    info <- tryCatch(png::readPNG(img_path, info = TRUE),
                     error = function(e) NULL)
    if (is.null(info)) return(list(width = target_width, height = target_width * 0.75))
    dims <- attr(info, "dim")  # height, width
    aspect <- dims[1] / dims[2]
    list(width = target_width, height = target_width * aspect)
  }

  # --- Prettify category names for headings ---
  pretty_category <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("(^|\\s)(\\w)", "\\1\\U\\2", x, perl = TRUE)
    x
  }

  # --- Build the document ---
  doc <- read_docx()

  # Title page
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(doc, paste("Run:", basename(run_dir)), style = "Normal")
  doc <- body_add_par(doc, paste("Generated:", Sys.time()), style = "Normal")
  doc <- body_add_break(doc, "page")

  # Table of contents
  doc <- body_add_par(doc, "Contents", style = "heading 1")
  doc <- body_add_toc(doc, level = 2)
  doc <- body_add_break(doc, "page")

  # --- Loop through categories (in order of appearance) ---
  categories <- unique(res$category)
  n_items    <- nrow(res)
  item_count <- 0

  for (i in seq_along(categories)) {

    cat_name <- categories[i]
    cat_data <- res[res$category == cat_name, ]

    # Section heading
    doc <- body_add_par(doc, pretty_category(cat_name), style = "heading 1")

    for (j in seq_len(nrow(cat_data))) {

      row   <- cat_data[j, ]
      fpath <- file.path(run_dir, row$basename)
      item_count <- item_count + 1

      if (!file.exists(fpath)) {
        message("  [", item_count, "/", n_items, "] File not found, skipping: ",
                row$basename)
        next
      }

      if (row$type == "plot") {
        # --- Insert image ---
        cat("  [", item_count, "/", n_items, "] ", row$basename, "\n")

        if (is.null(img_height)) {
          dims <- get_img_dims(fpath, img_width)
          w <- dims$width
          h <- dims$height
        } else {
          w <- img_width
          h <- img_height
        }

        doc <- body_add_img(doc, src = fpath, width = w, height = h,
                            style = "centered")

        # Caption below the image
        if (!is.na(row$caption) && nchar(trimws(row$caption)) > 0) {
          doc <- body_add_par(doc, trimws(row$caption), style = "Normal")
        }

        doc <- body_add_par(doc, "", style = "Normal")  # spacer

      } else if (row$type == "table") {
        # --- Insert table from CSV ---
        cat("  [", item_count, "/", n_items, "] ", row$basename, " (table)\n")

        tbl_data <- tryCatch(
          read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE),
          error = function(e) {
            message("  Could not read table: ", row$basename, " - ", e$message)
            NULL
          }
        )

        if (!is.null(tbl_data) && nrow(tbl_data) > 0) {

          # Sub-heading with the table name
          tbl_name <- tools::file_path_sans_ext(row$basename)
          tbl_name <- gsub("[_.]", " ", tbl_name)
          doc <- body_add_par(doc, tbl_name, style = "heading 2")

          # Build flextable
          ft <- flextable(tbl_data)
          ft <- autofit(ft)
          ft <- fontsize(ft, size = 8, part = "all")
          ft <- font(ft, fontname = "Arial", part = "all")
          ft <- bold(ft, part = "header")
          ft <- theme_box(ft)

          # Handle wide tables - set to page width
          ft <- set_table_properties(ft, layout = "autofit", width = 1)

          doc <- body_add_flextable(doc, value = ft)

          # Caption
          if (!is.na(row$caption) && nchar(trimws(row$caption)) > 0) {
            doc <- body_add_par(doc, trimws(row$caption), style = "Normal")
          }

          doc <- body_add_par(doc, "", style = "Normal")  # spacer
        }
      }
    }

    # Page break between categories (except after the last one)
    if (i < length(categories)) {
      doc <- body_add_break(doc, "page")
    }
  }

  # --- Save ---
  print(doc, target = output_file)
  cat("\nReport written to:", output_file, "\n")

  invisible(output_file)
}
