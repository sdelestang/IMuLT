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
  file.copy(src, dest, overwrite = FALSE)
  message("ModelStructure.xlsx copied to: ", dest)
}
