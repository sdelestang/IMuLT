.onLoad <- function(libname, pkgname) {
  # Only load DLL if we're in an installed package (not devtools::load_all)
  if(!is.null(libname)) {
    tryCatch(
      library.dynam("IMuLT", pkgname, libname),
      error = function(e) {
        # Silently fail - devtools loads it differently
      }
    )
  }

  # Source workflow scripts
  r_files <- system.file("R_files", package = "IMuLT")
  if(dir.exists(r_files)) {
    tryCatch({
      source(file.path(r_files, "BackgroundLoadOutput.R"), local = FALSE)
    }, error = function(e) {
      message("Note: Workflow scripts not yet available")
    })
  }
}

.onUnload <- function(libpath) {
  tryCatch(
    library.dynam.unload("IMuLT", libpath),
    error = function(e) invisible()
  )
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("IMuLT package loaded successfully.")
}
