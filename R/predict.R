#' Predict method using mars object
#'
#' @param object
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
#' predict.mars(mars.object)
#'
#'
predict.mars <- function(object,newdata) {
  
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  
  beta <- object$coefficients
  drop(B %*% beta)
}


make_B = function(data,Bfuncs){
  B = matrix(1,nrow=dim(data)[1],ncol=length(Bfuncs))
  
  for(i in 2:length(Bfuncs)){
    
    temp = Bfuncs[[i]]
    
    for(j in 1:(dim(temp)[1])){
      
      B[,i] = B[,i]*h(data[,temp[j,"v"]],temp[j,"s"],temp[j,"t"])
    }
  }
  return(B)
}
