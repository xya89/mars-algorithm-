library(rpart)
library(earth)







#' Multivariate Adaptive Regression Splines (MARS)
#'
#' This function takes a formula and a data set and uses MARS algorithm to
#' create a model for data analysis.
#'
#' The MARS model is a form of linear regression model that automatically models
#' nonlinearities and interactions between variables.
#'
#' @usage mars(formula, data, control)
#'
#' @param formula am R formula
#' @param data a data frame containing the data fro the model
#' @param control an object of class 'mars.control'
#'
#' @details
#' The function returns a model which is the linear combination of the base functions.
#' The model is fitted through a forward and backwards procedure.
#'
#'
#' @return an object of class 'mar.control'
#'
#'
#' @aliases Sean Yang
#'
#' @references
#' Friedman, Jerome H. "Multivariate Adaptive Regression Splines". The Annals of Statistics, Mar., 1991, Vol. 19, No.1 (Mar., 1991), pp. 1-67
#'
#' @import stats
#' @import ISLR
#'
#' @seealso
#' plot.mars() to plot a mars object
#' predict.mars() to make predictions for new data using a mars object
#' print.mars() to print information about a mars object
#' summary.mars() to print summary of a mars object.
#'
#'
#' @examples
#'
#' mm <- mars(wage ~ age, data = ISLR::Wage)
#'
#'
#' @export
#'
mars <- function(formula, data, control = NULL...) {

  cc <- match.call() # save the call
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  # drop the first column "intercept"
  x <- x[,-1]


  fwd <- fwd_stepwise(y, x, control)
  bwd <- bwd_stepwise(fwd, control)


  model <- lm(y ~.-1, data.frame(y=y, bwd$B))
  out <- c(list(call = cc, formula = formula, y=y, B=bwd$B, Bfuncs = bwd$Bfuncs,
                x_names = colnames(x)), model)
  class(out) <- c("mars", class(model))
  return(out)
}



# Constructor for mars.control
new_mars.control <- function(Mmax, d, trace) {
  structure(list(Mmax = Mmax, d = d, trace = trace), class = "mars.control")
}

# Validator
validate_mars.control <- function(control) {
  if(control$Mmax < 2) stop("Mmax needs to be greater than 2")
  if(class(control$d) != "numeric") stop("d is not numeric")
  if(class(control$trace) != "logical") stop("trace is not logical")
}

# Helper
mars.control <- function (Mmax = 2, d = 3, trace = FALSE) {

  control <- new_mars.control(Mmax, d, trace)
  validate_mars.control(control)
  return(control)

}




fwd_stepwise <- function(y, x, control){

  # Initialize:

  Mmax <- control$Mmax
  N <- length(y)
  n <- ncol(x)
  B <- init_B(N,Mmax)
  Bfuncs <- vector(mode = 'list', length = Mmax + 1)


  #---------------------------------------------------
  # Looping for forward selection:
  for(i in 1:(Mmax/2)) { # contrast to indexing 2...Mmax in Friedman

    M <- 2*i-1
    lof_best <- Inf

    for(m in 1:M) { # choose a basis function to split
      Xv <- setdiff(1:n, Bfuncs[[m]][,"v"])

      for(v in Xv){ # select a variable to split on
        tt <- split_points(x[,v], B[,m])

        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             Btem1 = B[,m]*h(x[,v], +1, t),
                             Btem2 = B[,m]*h(x[,v], -1, t))

          gdat <- data.frame(y = y,Bnew)

          lof <- LOF(y~., gdat, control)

          if(lof < lof_best) {
            lof_best <- lof
            Best_split <- c(m=m, v=v, t=t)
          }

        } # end loop over splits

      } # end loop over variables

    } # end loop over basis functions to split

    m <- Best_split["m"]; v <- Best_split["v"]; t <- Best_split["t"]


    B[ ,M+1] <- B[,m]*h(x[,v], -1, t)
    B[ ,M+2] <- B[,m]*h(x[,v], +1,  t)

    Bfuncs[[M+1]] <- rbind(Bfuncs[[m]], c(s=-1, v, t))
    Bfuncs[[M+2]] <- rbind(Bfuncs[[m]], c(s=+1, v, t))

  } # end loop over M
  colnames(B) <- paste0("B", (0:(ncol(B) - 1)))
  return(list(y=y, B=B,Bfuncs=Bfuncs))
}

bwd_stepwise <- function(fwd, control) {

  # Initialize
  Mmax = length(fwd$B) - 1
  Jstar = 2:(Mmax + 1)
  Kstar <- Jstar

  dat <- data.frame(y=fwd$y, fwd$B)

  LOFstar <- LOF(y~.-1, dat, control)

  # Loop

  for (M in (Mmax+1):2) {

    b <- Inf
    L <- Kstar

    for (m in L) {
      K = setdiff(L, m)

      dat.temp <- data.frame(y=fwd$y, fwd$B[,K])
      lof <- LOF(y~., dat.temp, control)

      if(lof < b){
        b <- lof
        Kstar <- K
      }

      if(lof < LOFstar){
        LOFstar <- lof
        Jstar <- K
      }
    }

  }

  Jstar = c(1, Jstar)
  return(list(y=fwd$y, B=fwd$B[,Jstar], Bfuncs = fwd$Bfunc[Jstar]))
}

LOF <- function(formula, data, control) {
  ff <- lm(formula, data)
  RSS <- sum(residuals(ff)^2)
  N <- nrow(data)
  M <- length(coefficients(ff)) - 1
  CM <- sum(hatvalues(ff))
  CM_Tilde <- CM + (control$d*M)
  GCV <- (RSS*N/(N-CM_Tilde)^2)
  return(GCV)
}


h <- function(x, s, t){
  return(pmax(0, s*(x-t)))
}

init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}









