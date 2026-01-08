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
