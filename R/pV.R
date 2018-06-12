#' @name pV
#' @title Present value of the future payments vector.
#' @description Calculate the present value of the future payments vector.
#' @param object Output from \code{separateProvision} or \code{glmProvision} functions.
#' @param yieldCurve Vector with the interest rate for the following years.
#' @return A list of 4 elements with:
#' \itemize{
#'   \item \code{actualpayments} Present value of the payment per origin year.
#'   \item \code{totactpay} Present value of the total payment.
#'   \item \code{yieldCurve} Input vector with interest rates.
#'   \item \code{fpv} Input vector with the future payments per origin year.
#' }
#' @export
#' @import stats

pV <- function(object, yieldCurve){
  # 1. Check input ----
  if (!(class(object) %in% c("glmprov", "sepprov"))){
    stop("Object is not an output from either 'glmProvision' or 'separateProvision' function.")
  }
  if (!(length(object$params$fpv) == length(yieldCurve))){
    stop("Dimensions of fpv and yieldCurve are not the same")
  }

  for (i in 1:length(yieldCurve)){
    if(!(is.numeric(yieldCurve[i]))){
      stop("yieldCurve is not a number")
    }
  }

  # 2. Actual value of future payments (discrete payments) ----
  discretepayments <- NULL
  for (i in 1:length(object$params$fpv)){
    discretepayments[i] <- object$params$fpv[i] * (1 + yieldCurve[i])^(-i)
  }
  totdiscpay <- sum(discretepayments)

  # 3. Actual value of future payments (continuous payments) ----
  fcont <- numeric(length(yieldCurve))
  for (i in 1:length(yieldCurve)){
    fcont[i] <- (1 - ((1+yieldCurve[i])^(-1)))/log(1+yieldCurve[i])
  }
  continuouspayments <- rep(1,length(yieldCurve))
  for (i in 2:length(yieldCurve)){
    continuouspayments[i]<- (1+yieldCurve[i])^(-i+1)
  }
  continuouspayments <- object$fpv * continuouspayments * fcont
  totcontpay <- sum(continuouspayments)

  # 4. Results object ----
  res <- list("discretepayments" = discretepayments,
              "totdiscpay" = totdiscpay,
              "continuouspayments" = continuouspayments,
              "totcontpay" = totcontpay,
              "yieldCurve" = yieldCurve,
              "fpv" = object$params$fpv)
  return(res)
}
