#' Summary function of the mars object
#'
#' @param marsobject
#'
#' @return
#' Print the summary of the mars object using print.mars
#'
#' @export
#'
#' @examples
#' summary.mars(obj)
#'
#' Call
#' ...
#'
#' Coefficients
#' ...
#'
#' Residuals errors
#' ....
#'
summary.mars <- function(marsobject) {

  print.mars(marsobject)


  print.default("\n Residuals errors")
  residuals.err <- (marsobject$residuals - mean(marsobject$residuals))^2/length(marsobject$residuals)

  cat(residuals.err)



}
