#' Print.mars
#'
#' @param marsobject
#'
#' @return print the fitted model and the coefficients
#'
#' @export
#'
#'
#' @examples
#' print.mars(mars)
#'
#' Call:
#' ...
#'
#' Coefficients:
#' ...
#'
#'
#'
print.mars <- function(marsobject) {

  cat("\n Call: \n")
  print.default(format(marsobject$call))

  cat("\n Coefficients: \n")
  print.default(format(marsobject$coefficients))
}
