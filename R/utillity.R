

#' @title get files in the package
#'
#' @keywords internal
pkg_file <- function(..., mustWork = TRUE) {
    system.file(..., package = "biozhuoer", mustWork = mustWork)
}


