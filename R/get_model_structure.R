

find_model_file <- function(filename = "ModelStructure.xlsx", up = 2L, down = 2L) {
  here  <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  roots <- Reduce(function(p, .) dirname(p), seq_len(up), accumulate = TRUE, init = here)
  roots <- unique(roots)

  pat  <- glob2rx(filename)
  hits <- character(0)
  for (r in roots) {
    f <- list.files(r, pattern = pat, recursive = TRUE,
                    full.names = TRUE, ignore.case = TRUE)
    if (!length(f)) next
    f   <- normalizePath(f, winslash = "/", mustWork = FALSE)
    rel <- substring(f, nchar(r) + 2L)            # path below root, sans leading "/"
    dep <- lengths(strsplit(rel, "/")) - 1L       # folders between root and file
    hits <- c(hits, f[dep <= down])
  }
  hits <- unique(hits)
  if (!length(hits))
    stop(sprintf("Could not find '%s' within %d folder(s) up/down of:\n  %s",
                 filename, max(up, down), here))
  if (length(hits) > 1L)
    warning(sprintf("Multiple copies of '%s' found; using the first:\n%s",
                    filename, paste(" -", hits, collapse = "\n")))
  hits[[1]]
}

#3 Laods ModelStucture whther it is open or not
load_model_structure <- function(filename = "ModelStructure.xlsx",
                                 tries = 8, wait = 0.25) {
  path <- find_model_file(filename)

  retry <- function(expr, what) {
    for (i in seq_len(tries)) {
      ok <- tryCatch(expr, error = function(e) e)
      if (!inherits(ok, "error") && !identical(ok, FALSE)) return(ok)
      Sys.sleep(wait)
    }
    stop(sprintf("'%s' still failing after %d tries - OneDrive may be mid-sync. ",
                 basename(path), tries),
         "Right-click the folder > 'Always keep on this device', or pause OneDrive and retry.")
  }

  # copy to LOCAL temp (not the synced folder) so the read path never touches OneDrive
  tmp <- file.path(tempdir(), sprintf("RtmpCopy_%d_%s", Sys.getpid(), basename(path)))
  on.exit(unlink(tmp), add = TRUE)

  retry(file.copy(path, tmp, overwrite = TRUE), "copy")  # forces hydration; retries past sync locks
  wb <- retry(loadWorkbook(tmp), "load")                 # loads a local, unsynced file
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
  if(length(f)>0)  Archive <- toupper(dlg_input(paste0("ModelStructure.xlsx already exists. Archive it first? (Y or N)")
  )$res)
  if (Archive == 'Y') {
    invisible(file.rename("ModelStructure.xlsx", "ModelStructureArchive.xlsx"))
    print("Archived old ModelStructure.xlsx")
  }
  file.copy(src, dest, overwrite = FALSE)
  message("ModelStructure.xlsx copied to: ", dest)
}

