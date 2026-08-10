
find_model_dir <- function(pattern = "Run", up = 3L) {
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (i in 0:up) {
    hits <- list.files(d, pattern = pattern)
    if (length(hits) > 0) return(d)
    parent <- dirname(d)
    if (parent == d) break          # hit the drive root
    d <- parent
  }
  stop("Could not find any items matching pattern '", pattern,
       "' in the working directory or its ", up, " parents.")
}


find_model_file <- function(filename = "ModelStructure.xlsx", up = 3L, highest = FALSE) {
  d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  best <- NULL
  for (i in 0:up) {
    cand <- file.path(d, filename)
    if (file.exists(cand)) {
      if (!highest) return(cand)   # old behaviour: nearest match, unchanged
      best <- cand                  # highest behaviour: keep overwriting -> ends up topmost hit
    }
    parent <- dirname(d)
    if (parent == d) break          # hit the drive root
    d <- parent
  }
  if (!is.null(best)) return(best)
  stop("Could not find '", filename, "' in the working directory or its ", up, " parents.")
}

#3 Loads ModelStucture whether it is open or not
load_model_structure <- function(filename = "ModelStructure.xlsx",
                                 tries = 5, wait = 0.25) {
  path <- find_model_file(filename)
  tmp  <- file.path(tempdir(), sprintf("RtmpCopy_%d_%s", Sys.getpid(), basename(path)))
  on.exit(unlink(tmp), add = TRUE)

  for (i in seq_len(tries)) {
    ok <- tryCatch(fs::file_copy(path, tmp, overwrite = TRUE), error = function(e) e)
    if (!inherits(ok, "error")) break
    if (i == tries)
      stop("Could not copy '", basename(path), "': ", conditionMessage(ok))
    Sys.sleep(wait)
  }

  wb <- loadWorkbook(tmp)               # loads a local, unsynced copy
  attr(wb, "source_path") <- path
  wb
}


#' Get path to ModelStructure.xlsx template
#'
#' @export
get_model_structure_path <- function() {
  system.file("extdata", "ModelStructure.xlsx", package = "IMuLT")
}

#' Copy ModelStructure.xlsx to current directory
#'
#' @param dest_path Destination path (default: current directory)
#' @export
copy_model_structure <- function(dest_path = ".") {
  src <- get_model_structure_path()
  dest <- file.path(dest_path, "ModelStructure.xlsx")
  f <- list.files(dest_path, include.dirs = FALSE, full.names = TRUE, recursive = TRUE, pattern = "ModelStructure.xlsx")
  Archive <- 'N'
  if(length(f)>0)  Archive <- toupper(svDialogs::dlg_input(paste0("ModelStructure.xlsx already exists. Archive it first? (Y or N)")
  )$res)
  if (Archive == 'Y') {
    invisible(file.rename("ModelStructure.xlsx", "ModelStructureArchive.xlsx"))
    print("Archived old ModelStructure.xlsx")
  }
  file.copy(src, dest, overwrite = FALSE)
  message("ModelStructure.xlsx copied to: ", dest)
}

