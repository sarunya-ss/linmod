#' Linear Regression
#' @param formula = y and x
#' @param data = data
#' @export

linmod <- function(formula, data=list()) {
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  # compute QR-decomposition of x
  qx <- qr(x)
  # compute (x’x)^(-1) x’y
  coef <- solve.qr(qx, y)
  ## degrees of freedom and standard deviation of residuals
  df <- nrow(x)-ncol(x)
  sigma2 <- sum((y - x%*%coef)^2)/df
  # compute sigma^2 * (x’x)^-1
  vcov <- sigma2 * chol2inv(qx$qr)
  est <- list(coefficients = coef, vcov = vcov, sigma = sqrt(sigma2), df = df)
  est$fitted.values <- as.vector(x %*% est$coefficients)
  est$residuals <- y - est$fitted.values
  est$call <- match.call()
  est$formula <- formula
  class(est) <- "linmod"
  est
}

#' Summary Linear Regression
#' @param object = object
#' @param ... = parameter
#' @export

summary.linmod <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  res
}
