.onLoad <- function(libname, pkgname) {
  library.dynam("IMuLT", pkgname, libname)
}

.onUnload <- function(libpath) {
  library.dynam.unload("IMuLT", libpath)
}

.onAttach <- function(libname, pkgname) {
  # Load required packages
  require(TMB, quietly = TRUE)
  require(readxl, quietly = TRUE)
  require(reshape2, quietly = TRUE)
  require(svDialogs, quietly = TRUE)
  require(stats4, quietly = TRUE)
  require(makehtml, quietly = TRUE)
  require(hplot, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(magrittr, quietly = TRUE)
  require(tidyr, quietly = TRUE)
  require(ggplot2, quietly = TRUE)
  require(openxlsx, quietly = TRUE)

  # Load the TMB DLL
  dllpath <- system.file("libs", package = "IMuLT", lib.loc = libname)
  dyn.load(TMB::dynlib(file.path(dllpath, .Platform$r_arch, "IMuLT")))

  # Auto-source all workflow files
  r_files <- system.file("R_files", package = "IMuLT")
  source(file.path(r_files, "Utilities.R"), local = FALSE)
  source(file.path(r_files, "ReadMaterial.R"), local = FALSE)
  source(file.path(r_files, "WriteMaterial.R"), local = FALSE)
  source(file.path(r_files, "Estimate.R"), local = FALSE)
  source(file.path(r_files, "WriteDat.R"), local = FALSE)
  source(file.path(r_files, "MakeOutPut.R"), local = FALSE)
  source(file.path(r_files, "LoadOutputData.R"), local = FALSE)
  source(file.path(r_files, "FileBuilder.R"), local = FALSE)

  packageStartupMessage("IMuLT package loaded. All workflow functions available.")
}
