
#' Plot method to plot a mars object
#'
#' @param marsobject
#' @param ...
#'
#' @return
#' the one way interaction and two way interaction plot for mars object output
#'
#' @export
#'
#' @examples
#' plot.mars(mars)
#'
#' out:
#' ...
#' ...
#'
plot.mars <- function(marsobject, ...) {


  data <- eval(marsobject$call$data)
  tt <- terms(marsobject$formula, data = data)
  tt <- delete.response(tt)
  mf <- model.frame(tt, data)
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)[,-1]


  Bfunc <- marsobject$Bfuncs
  singleB <- which(sapply(Bfunc, function(x) NROW(x) == 1))
  doubleB <- which(sapply(Bfunc, function(x) NROW(x) == 2))


  totalfuncs <- length(singleB) + length(doubleB)

  for(i in singleB) {
    vv <- Bfunc[[i]][1,"v"]

    xx <- seq(from=min(X[,vv]),to=max(X[,vv]),length=100)

    bb <- h(Bfunc[[i]][1,"s"],xx,Bfunc[[i]][1,"t"])

    plot(xx,bb,type="l",...)
  }


  for(i in doubleB) {
    vv1 <- Bf[[i]][1,"v"]
    varname1 <- x$x_names[[vv1]]

    vv2 <- Bf[[i]][2,"v"]
    varname2 <- x$x_names[[vv2]]

    xx <- seq(from=min(X[,vv1]),to=max(X[,vv1]),length=100)
    yy <- seq(from=min(X[,vv2]),to=max(X[,vv2]),length=100)

    ff <- function(x,y) {
      return (h(Bf[[i]][1,"s"],x,Bf[[i]][1,"t"])*
                h(Bf[[i]][2,"s"],y,Bf[[i]][2,"t"]))

    }

    zz <- outer(xx, yy, FUN=ff)

    persp(xx, yy, zz, xlab=varname1, ylab=varname2, zlab="",
          main = paste0("B", i), theta = -30, phi=30,
          col="lightblue", lwd=.1)

  }


}
