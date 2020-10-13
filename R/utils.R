#' Loads dataset from package
#'
#' @export
#'
#' @importFrom utils data
#' @param name_data name of the data as a string
#' @param pkg name of the package where the data is as a string
load_named_data <- function(name_data, pkg) {
  # check if package is installed
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(get(utils::data(list = name_data, package = pkg, envir = environment())))
  } else {
    stop("Install package '", pkg, "'.")
  }
}
