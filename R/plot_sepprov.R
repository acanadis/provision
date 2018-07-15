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
  hist(x = sum.test, main = "Histogram of total claims cost", xlab = "Total claims cost")
  rm(test, sum.test, boots)

  # Boxplot - Latest - OY
  latest <- matrix(0, ncol = nrow(object$params$increm.tri),
                   nrow = dim(object$params$increm.tri)[3])
  for (i in 1:dim(object$params$increm.tri)[3]){
    latest[i, ] <- ObtainMDiagonal(object$params$triangboot[, , i])
  }
  aux <- reshape2::melt(latest)
  colnames(aux) <- c("B", "oy", "latest")
  aux$oy <- factor(aux$oy)
  if (log == TRUE){
    aux$latest <- log(aux$latest)
    boxplot(aux$latest ~ aux$oy,
            main = "Latest actual incremental claim against simulated values",
            xlab = "Origin year", ylab = "log(Latest incremental claims)",
            pch = 20, border = "grey", outline = FALSE)
  } else {
    boxplot(aux$latest ~ aux$oy,
            main = "Latest actual incremental claim against simulated values",
            xlab = "Origin year", ylab = "Latest incremental claims",
            pch = 20, border = "grey", outline = FALSE)
  }
  aux <- data.frame(Value = ObtainMDiagonal(object$triangle),
                    Index = c(0:(ncol(object$triangle)-1)),
                    stringsAsFactors = FALSE)
  if (log == TRUE) {
    aux$Value  <- log(aux$Value)
  }
  points(x = aux$Index, y = aux$Value, type = "p", pch = 20, col = "red")
  rm(latest, i, aux)

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
          main = "Simulated ultimate claims cost",
          xlab = "Origin year", ylab = "Ultimate claims costs",
          pch = 20, border = "grey", outline = FALSE)
  points(auxoy$oy, auxoy$mu, type = "p", pch = 20, col = "red")
  rm(ultimate, i, aux)

  # Boxplot - FPV - CY
  aux <- object$params$fpvfutureboot
  colnames(aux) <- c(as.numeric(colnmes(aux) + ncol(object$params$increm.tri)))
  boxplot(aux,
          main = "Actual future payments per calendar year against simulated values",
          xlab = "Calendar year", ylab = "Future payments",
          pch = 20, border = "grey", outline = FALSE)
  aux <- data.frame(Value = object$params$fpv,
                    Index = c(1:ncol(object$params$fpvfutureboot)),
                    stringsAsFactors = FALSE)
  points(x = aux$Index, y = aux$Value, type = "p", pch = 20, col = "red")
  rm(aux)
}
