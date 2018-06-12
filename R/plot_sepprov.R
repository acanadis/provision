#' @name plot_sepprov
#' @title Plots for separateProvision objects
#' @description Displays plots of a separateProvision object
#' @param object separateProvision object to plot
#' @param ... other arguments ignored (for compatibility with generic)
#' @export
#' @importFrom dplyr group_by summarise
#' @importFrom reshape2 melt
#' @importFrom graphics boxplot hist par points
#' @include utilities.R

plot_sepprov <- function(object, ...){
  par(mfrow = c(2,2))
  # Histogram - IBNR
  sum.test <- NULL
  test <- NULL
  for (boots in 1:dim(object$params$increm.tri)[3]){
    aux <- object$params$increm.tri[, , boots]
    for (i in 1:nrow(aux)){
      test <- c(test, aux[i, ncol(aux)] - aux[i, (ncol(aux)-i+1)])
    }
    sum.test <- c(sum.test, sum(test))
    test <- NULL
  }
  hist(x = sum.test, main = "Histogram of total IBNR", xlab = "IBNR")
  rm(test, sum.test, boots)
  # Boxplot - Ultimate
  ultimate <- matrix(0, ncol = nrow(object$params$increm.tri),
                     nrow = dim(object$params$increm.tri)[3])
  for (i in 1:dim(object$params$increm.tri)[3]){
    ultimate[i, ] <- object$params$increm.tri[, ncol(object$params$increm.tri), i]
  }
  aux <- reshape2::melt(ultimate)
  colnames(aux) <- c("B", "oy", "ultimate")
  aux$oy <- factor(aux$oy)
  auxoy <- aux %>% dplyr::group_by(oy) %>% dplyr::summarise(mu = mean(ultimate))
  boxplot(aux$ultimate ~ aux$oy,
          main = "Simulated ultimate claim cost",
          xlab = "Origin year", ylab = "Claim cost",
          pch = 20, border = "grey", outline = FALSE)
  points(auxoy$oy, auxoy$mu, type = "p", pch = 20, col = "red")
  rm(ultimate, i)
  # Boxplot - Latest - OY
  latest <- matrix(0, ncol = nrow(object$params$increm.tri),
                   nrow = dim(object$params$increm.tri)[3])
  for (i in 1:dim(object$params$increm.tri)[3]){
    latest[i, ] <- ObtainMDiagonal(object$params$increm.tri[, , i])
  }
  aux <- reshape2::melt(latest)
  colnames(aux) <- c("B", "oy", "latest")
  aux$oy <- factor(aux$oy)
  boxplot(aux$latest ~ aux$oy,
          main = "Simulated latest per origin year",
          xlab = "Origin year", ylab = "Latest",
          pch = 20, border = "grey", outline = FALSE)
  rm(latest, i)
  # Boxplot - FPV - CY
  boxplot(object$params$fpvfutureboot,
          main = "Simulated future payments per calendar year",
          xlab = "Calendar year", ylab = "Future payments",
          pch = 20, border = "grey", outline = FALSE)
  points(colMeans(object$params$fpvfutureboot), type = "p", pch = 20, col = "red")
}
