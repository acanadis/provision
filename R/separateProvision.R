#' @name separateProvision
#' @title Calculate provisions using separation methods of Taylor.
#' @description Calculate provisions using separation methods of Taylor.
#' @usage separateProvision(lossData, freqData,
#' modelSep = "arithmetic", lambdaK = 0,
#' B = 1000, seed = NULL)
#' @param lossData Matrix of incremental losses \eqn{Cij},
#' for \eqn{i = 1,...,t} origin years (rows) and for \eqn{j = 1,...,t}
#' development years (columns); filled with \code{NAs} for \eqn{i + j > t}.
#' @param freqData Vector with the claim numbers by origin year \eqn{i}
#' for \eqn{i = 1,...t}.
#' @param modelSep Model to be used, can be \code{arithmetic} or \code{geometric}. Defaults to "arithmetic".
#' @param lambdaK Percentage of the trend in the inflation index. Defaults to 0.
#' @param B Number of iterations to perform in the bootstrapping procedure. Defaults to 1000.
#' @param seed Seed to make bootstrap reproducible.
#' @return A list of 5 elements with:
#' \itemize{
#'   \item \code{triangle} Input data (lossData).
#'   \item \code{glm.triangle} Fitted values by \code{glm} modeling.
#'   \item \code{model} \code{glm} model.
#'   \item \code{glm.triangle.bootstrap} Three-dimensional array of lower triangles from bootstrap procedure.
#'   \item \code{summary} Summary table.
#'   \item \code{params} List of output parameters.
#' }
#' @examples
#' data("TaylorData")
#' separateProvision(lossData = TaylorData$lossData,
#'                   freqData = TaylorData$freqData)
#' @export
#' @import statmod
#' @importFrom stats glm rgamma
#' @include utilities.R

devtools::use_package("statmod")

separateProvision <- function(lossData, freqData,
                              modelSep = "arithmetic",
                              lambdaK = 0,
                              B = 1000,
                              seed = NULL){
  # 1. Check input ----
  if (! (class(lossData) == "matrix")){
    stop("Object lossData is not a matrix.")
  }
  check_Data(lossData)
  if (! (class(freqData) == c("numeric"))){
    stop("Object freqData is not a numeric vector.")
  }
  # Function arguments
  if (! (modelSep %in% c("arithmetic", "geometric"))){
    stop("Arg modelSep must be either one of 'arithmetic', 'geometric'.")
  }

  # 2. Common variables ----
  lambdaK <- lambdaK/100

  call <- match.call()
  t <- nrow(lossData)
  oy <- rep(1:t,t:1)
  dy <- sequence(t:1)
  labelsoy <- paste0("oy", 0:(t-1))
  labelsdy <- paste0("dy", 0:(t-1))
  colnames(lossData) <- labelsdy
  rownames(lossData) <- labelsoy

  # Vector of average losses
  freq.data.oy <- rep(freqData, t:1)
  v.lossData <- as.vector(t(lossData))
  v.lossData <- v.lossData[!is.na(v.lossData)]
  average.loss.data <- v.lossData/freq.data.oy

  # Triangle of average incremental losses
  m.average.loss.data <- matrix(0, t, t)
  for (i in 1:t){
    m.average.loss.data[i,] <- c(average.loss.data[oy == i],
                                 rep(times = i-1, NA))
  }

  labelsoy <- paste0("oy", 0:(t-1))
  labelsdy <- paste0("dy", 0:(t-1))
  colnames(m.average.loss.data) <- labelsdy
  rownames(m.average.loss.data) <- labelsoy

  oy <- as.factor(oy); dy <- as.factor(dy)
  cy <- as.numeric(oy) + as.numeric(dy); cy <- as.factor(cy)

  # 3. Triangle (GLM) ----
  # 3.1. Arithmetic ----
  if (modelSep == "arithmetic"){
    glmArit <- glm(average.loss.data ~ dy + cy,
                   family = statmod::tweedie(var.power = 1, link.power = 0))
    # Triangle of incremental losses fitted
    m.average.fitted.loss.data <- matrix(0, t, t)
    for (i in 1:t){
      m.average.fitted.loss.data[i, ] <- c(glmArit$fitted.values[oy == i],
                                           rep(times = i-1, NA))
    }
    m.fitted.loss.data <- m.average.fitted.loss.data * freqData

    triangglm <- matrix(0, t, t)
    for (i in 2:t){
      aux <- glmArit$fitted.values[((oy == (t-i+1)) & (dy == i))] *
        freqData[(t-i+2):t] * (1+lambdaK)^(1:(i-1))
      triangglm[, i] <- c(lossData[1:(t-i+1), i], aux)
      rm(aux)
    }
    triangglm[, 1] <- lossData[, 1]
  }
  # 3.2. Geometric ----
  if (modelSep == "geometric"){
    glmGeom <- glm(log(average.loss.data) ~ dy + cy,
                   family = statmod::tweedie(var.power = 0,link.power = 1))

    # Triangle of incremental losses fitted
    m.average.fitted.loss.data <- matrix(0, t, t)
    for (i in 1:t){
      m.average.fitted.loss.data[i, ] <- c(glmGeom$fitted.values[oy == i],
                                           rep(times = i-1, NA))
    }
    m.fitted.loss.data <- exp(m.average.fitted.loss.data) * freqData
    triangglm <- matrix(0, t, t)

    for (i in 2:t){
      aux <- exp(glmGeom$fitted.values[((oy == (t-i+1)) & (dy == i))])*
        freqData[(t-i+2):t] * (1 + lambdaK)^(1:(i-1))
      triangglm[, i] <- c(lossData[1:(t-i+1), i], aux)
      rm(aux)
    }
    triangglm[, 1] <- lossData[, 1]
  }
  # 4. OY Reserve ----
  oyres <- rep(0, t-1)
  for (i in 2:t){
    oyres[i-1] <- sum(triangglm[i, (t-i+2):t])
  }
  # Total reserve
  totres <- sum(oyres)
  # Vector of future payments
  fpv <- rep(0, dim(triangglm)[1] - 1)
  for (k in 1:dim(triangglm)[1] - 1) {
    future <- row(triangglm) + col(triangglm) - 1 == dim(triangglm)[1] + k
    fpv[k] <- sum(triangglm[future])
  }

  # 5. Lambda ----
  # 5.1. Arithmetic ----
  if (modelSep == "arithmetic"){
    Xij.mat <- matrix(0, t, t)
    for (k in 1:length(v.lossData)){
      Xij.mat[oy[k], dy[k]] <- v.lossData[k]
      ni <- matrix(rep(freqData, each = t), nrow = t, ncol = t, byrow = TRUE)
      sij <- Xij.mat/ni
    }
    sum.col <- colSums(sij)

    sum.diag <- rep(0, t)
    for (k in 1:t){
      diag.sij <- row(sij) + col(sij) - 1 == k
      sum.diag[k] <- sum(sij[diag.sij])
    }
    rm(diag.sij)
    lambda <- rep(0, t)
    r <- rep(0, t)
    lambda[t] <- sum.diag[t]
    r[t] <- sum.col[t]/lambda[t]
    for (i in 1:(t-1)){
      lambda[t-i] <- sum.diag[t-i]/(1-sum(r[(t-i+1):t]))
      r[t-i] <- sum.col[t-i]/sum(lambda[(t-i):t])
    }
    lambdafut <- lambda[t] * (1+lambdaK)^(1:(t-1))
    lambdatot <- c(lambda, lambdafut)
  }
  # 5.2. Geometric ----
  if (modelSep == "geometric"){
    Xij.mat <- matrix(0, t, t)
    for (k in 1:length(v.lossData)){
      Xij.mat[oy[k], dy[k]] <- v.lossData[k]
      ni <- matrix(rep(freqData, each = t), nrow = t, ncol = t, byrow = TRUE)
      sij <- Xij.mat/ni
    }
    prod.col <- rep(0, t)
    for (k in 1:t){
      prod.col[k] <- prod(sij[(1:(t+1-k)), k])
    }
    prod.diag <- rep(1, t)
    for (k in 1:t){
      diag.sij <- row(sij) + col(sij) - 1 == k
      prod.diag[k] <- prod(sij[diag.sij])
    }
    lambda <- rep(0, t)
    r <- rep(0, t)
    lambda[t] <- prod.diag[t]^(1/t)
    r[t] <- prod.col[t]/lambda[t]
    for (i in 1:(t-1)){
      lambda[t-i] <- (prod.diag[t-i]*prod(r[(t-i+1):t]))^(1/(t-i))
      r[t-i] <- (prod.col[t-i]/prod(lambda[(t-i):t]))^(1/(i+1))
    }
    lambdafut <- lambda[t] * (1+lambdaK)^(1:(t-1))
    lambdatot <- c(lambda, lambdafut)
  }

  # 6. Bootstraping ----
  matriz.lambda <- matrix(0, t, t)
  for (i in 1:t){
    matriz.lambda[, i] <- lambdatot[i:(t+i-1)]
  }
  r.vector <- rep(r, each = t)
  r.mat <- matrix(r.vector, t, t)
  cijboot <- array(dim = c(t, t, B), data = 0)

  if (modelSep == "arithmetic"){
    phi <- with(glmArit, sum(residuals^2/df.residual))
  }
  if (modelSep == "geometric"){
    phi <- with(glmGeom, sum(residuals^2/df.residual))
  }

  for (i in 1:t){
    for (j in 1:t){
      cijboot[i, j, ] <- rgamma(B, ni[i, j]/phi,
                                1/(r.mat[i,j] * matriz.lambda[i,j] * phi))
    }
  }

  # Triangles with bootstrapped incremental losses
  triangboot <- array(dim = c(t, t, B), data = 0)
  for (boots in 1:B){
    if(! is.null(seed)){
      set.seed(seed)
    }
    for(i in 1:t){
      triangboot[i, , boots] <- c(cijboot[i, 1:(t-i+1), boots],
                                  rep(times = i-1, NA))
    }
  }

  # Future payments bootstrapped
  futureboot <- array(dim = c(t, t, B), data = NA)
  for (boots in 1:B){
    if(! is.null(seed)){
      set.seed(seed)
    }
    for (i in 2:t){
      futureboot[i, , boots] <- c(rep(times = (t-i+1), NA),
                                  cijboot[i, (t-i+2):t, boots])
    }
  }

  m.bootData <- matrix(0, t, t)
  lossDataboot <- matrix(NA, t, t)
  triangglmboot <- array(data = 0, dim = c(t, t, B))
  oyresboot <- matrix(0, B, t-1)
  totresboot <- rep(0, times = B)
  fpvboot <- matrix(0, B, t-1)

  if (modelSep == "arithmetic"){
    # 6.1. Arithmetic on bootstrap triangles  ----
    for (boots in 1:B) {
      if(! is.null(seed)){
        set.seed(seed)
      }
      m.bootData <- t(triangboot[, , boots])
      v.bootData <- as.vector(m.bootData)
      v.bootData <- v.bootData[!is.na(v.bootData)]

      average.boot.data <- v.bootData/freq.data.oy

      glmAritboot <- glm(average.boot.data ~ dy + cy,
                         family = statmod::tweedie(var.power = 1, link.power = 0))

      for (i in 2:t){
        aux <- glmAritboot$fitted.values[((oy == (t-i+1)) & (dy == i))] *
          freqData[(t-i+2):t] * (1+lambdaK)^(1:(i-1))
        triangglmboot[, i, boots] <- c(lossDataboot[1:(t-i+1), i], aux)
        rm(aux)
      }
      triangglmboot[, 1, boots] <- lossDataboot[, 1]
    }
  }

  if (modelSep == "geometric"){
    # 6.2. Geometric on bootstrap triangles  ----
    for (boots in 1:B) {
      if(! is.null(seed)){
        set.seed(seed)
      }
      m.bootData <- t(triangboot[, , boots])
      v.bootData <- as.vector(m.bootData)
      v.bootData <- v.bootData[!is.na(v.bootData)]

      average.boot.data <- v.bootData/freq.data.oy

      glmGeomboot <- glm(log(average.boot.data) ~ dy + cy,
                         family = statmod::tweedie(var.power = 1,link.power = 0))
      for (i in 2:t){
        aux <- exp(glmGeomboot$fitted.values[((oy==(t-i+1))&(dy==i))]) *
          freqData[(t-i+2):t] * (1+lambdaK)^(1:(i-1))
        triangglmboot[, i, boots] <- c(lossDataboot[1:(t-i+1), i], aux)
        rm(aux)
      }
      triangglmboot[, 1, boots] <- lossDataboot[,1]
    }
  }

  # 7. Prediction errors ----
  for (boots in 1:B){
    if(! is.null(seed)){
      set.seed(seed)
    }
    # Origin year reserve
    for (i in 2:t){
      oyresboot[boots, i-1] <- sum(triangglmboot[i, (t-i+2):t, boots])
    }
    # Total reserve
    totresboot[boots] <- sum(oyresboot[boots, ])
    # Vector of future payments
    fpv1 <- rep(0, dim(triangglmboot)[1] - 1)
    for (k in 1:dim(triangglmboot)[1] - 1) {
      future <- row(triangglmboot[, , 1]) + col(triangglmboot[, , 1]) - 1 == nrow(triangglmboot) + 1
      fpv1[k] <- sum(triangglmboot[future])
      fpvboot[boots, ] <- fpv1
    }
  }

  # Computations with future payments bootstrapped
  oyresfutureboot <- matrix(0, B, t-1)
  totresfutureboot <- rep(0, times = B)
  fpvfutureboot <- matrix(0, B, t-1)

  for (boots in 1:B){
    if(! is.null(seed)){
      set.seed(seed)
    }
    m.bootfutureData <- futureboot[, , boots]
    # Origin year reserve
    for (i in 2:t){
      oyresfutureboot[boots, i-1] <- sum(m.bootfutureData[i, (t-i+2):t])
    }
    # Total reserve
    totresfutureboot[boots] <- sum(oyresfutureboot[boots, ])
    # Vector of future payments
    a <- ncol(m.bootfutureData) + 1
    fpv2 <- NULL
    for (z in 1:(nrow(m.bootfutureData)-1)){
      aux <- 0
      for (i in 1:nrow(m.bootfutureData)){
        for (j in 1:ncol(m.bootfutureData)){
          if ((i+j-1) == a){
            aux <- sum(aux, m.bootfutureData[i, j], na.rm = TRUE)
          }
        }
      }
      fpv2[z] <- aux
      a <- a + 1
    }
    fpvfutureboot[boots, ] <- fpv2
  }

  # Computation of prediction errors
  peorigin <- matrix(0, B, t-1)
  pefpv <- matrix(0, B, t-1)
  for (boots in 1:B){
    if(! is.null(seed)){
      set.seed(seed)
    }
    for(i in 1:(t-1)){
      peorigin[boots, ] <- (oyresfutureboot[boots, ] - oyresboot[boots, ])
      pefpv[boots, ] <- fpvfutureboot[boots, ]-fpvboot[boots, ]
    }
  }
  petot <- totresfutureboot - totresboot

  # Predictive distributions
  oyresmatrix <- matrix(data = oyres, B, t-1, byrow = TRUE)
  fpvmatrix<-matrix(data=fpv, B, t-1, byrow = TRUE)
  predoyres <- oyresmatrix + peorigin
  predfpv <- fpvmatrix + pefpv
  predtotres <- totres + petot

  # 8. Summaries ----
  # 8.1. Summary creation - OY ----
  out.sum <- matrix(NA, ncol = 10, nrow = t+1,
                    dimnames = list(c(rownames(lossData), "TOTAL.oy"),
                                    c("Latest", "dev.to.date", "Ultimate",
                                      "IBNR", "IBNR.mean", "PredErr.Abs", "CV",
                                      "IBNR.quantile.75", "IBNR.quantile.95",
                                      "IBNR.quantile.995")))
  aux <- array(0, c(t, t, B))
  auxup <- triangboot; auxdn <- triangglmboot
  auxup[is.na(auxup)] <- 0; auxdn[is.na(auxdn)] <- 0
  for (i in 1:B){
    aux[, , i] <- auxup[, , i] + auxdn[, , i]
  }
  rm(auxup, auxdn, i)
  increm.tri <- Increm.triangle(aux)
  ## Latest
  latest <- matrix(NA, ncol = nrow(lossData), nrow = B)
  for (i in 1:B){
    latest[i, ] <- c(ObtainMDiagonal(increm.tri[, , i]))
  }
  out.sum[, 1] <- c(colMeans(latest), sum(colMeans(latest)))
  rm(latest, i)
  ## Ultimate
  ultimate <- NULL
  for (i in 1:nrow(aux)){
    ultimate <- c(ultimate, mean(increm.tri[i, ncol(increm.tri), ]))
  }
  out.sum[, 3] <- c(ultimate, sum(ultimate))
  rm(ultimate, i)
  ## IBNR
  out.sum[, 4]<-c(0, oyres, sum(oyres))
  ## meanIBNR
  out.sum[, 5] <- c(0, colMeans(oyresboot), mean(totresboot))
  ## dev.to.date
  out.sum[, 2] <- out.sum[, 1]/out.sum[, 3]
  ## PE
  out.sum[, 6] <- c(0, colMeans(abs(peorigin)), mean(abs(petot)))
  ## CV
  out.sum[ , 7] <- out.sum[, 6]/out.sum[, 4]
  ## Quantiles
  out.sum[, 8] <- c(0, apply(oyresboot, 2, quantile, 0.75),
                    quantile(totresboot, 0.75))
  out.sum[, 9] <- c(0, apply(oyresboot, 2, quantile, 0.95),
                    quantile(totresboot, 0.95))
  out.sum[, 10] <- c(0, apply(oyresboot, 2, quantile, 0.995),
                     quantile(totresboot, 0.995))
  out.sum[is.nan(out.sum)] <- 0

  # 8.2. Summary - CY ----
  labelscy <- paste0("cy", (t+1):(t+ncol(lossData)-1))
  out.sum2 <- matrix(NA, ncol = 7, nrow = t,
                     dimnames = list(c(labelscy, "TOTAL.cy"),
                                     c("IBNR", "IBNR.mean", "PredErr.Abs", "CV",
                                       "IBNR.quantile.75", "IBNR.quantile.95",
                                       "IBNR.quantile.995")))
  out.sum2[, 1] <- c(fpv, sum(fpv))
  out.sum2[, 2] <- c(apply(fpvfutureboot,2,mean), sum(apply(fpvfutureboot,2,mean)))
  out.sum2[, 3] <- c(apply(abs(pefpv),2,mean), sum(apply(abs(pefpv),2,mean)))
  out.sum2[, 4] <- out.sum2[,3]/out.sum2[,2]
  out.sum2[, 5] <- c(apply(fpvfutureboot, 2, quantile, 0.75, na.rm = TRUE),
                     quantile(totresfutureboot, 0.75, na.rm = TRUE))
  out.sum2[, 6] <- c(apply(fpvfutureboot, 2, quantile, 0.95, na.rm = TRUE),
                     quantile(totresfutureboot, 0.95, na.rm = TRUE))
  out.sum2[, 7] <- c(apply(fpvfutureboot, 2, quantile, 0.995, na.rm = TRUE),
                     quantile(totresfutureboot, 0.995, na.rm = TRUE))

  out.sum2[is.nan(out.sum2)] <- 0

  # 9. Results object ----
  if (modelSep == "arithmetic"){
    res <- list(call = match.call(expand.dots = FALSE),
                triangle = lossData,
                freqData = freqData,
                glm.triangle = triangglm,
                model = glmArit,
                glm.triangle.bootstrap = triangglmboot,
                summary.oy = out.sum,
                summary.cy = out.sum2,
                params = list("modelSep" = modelSep,
                              "lambdaK" = lambdaK,
                              "B" = B,
                              "seed" = seed,
                              "fpv" = fpv,
                              "fpvfutureboot" = fpvfutureboot,
                              "increm.tri" = increm.tri,
                              "triangboot" = triangboot,
                              "lambdafut" = lambdafut))
  }
  if (modelSep == "geometric"){
    res <- list(call = match.call(expand.dots = FALSE),
                triangle = lossData,
                freqData = freqData,
                glm.triangle = triangglm,
                model = glmGeom,
                glm.triangle.bootstrap = triangglmboot,
                summary.oy = out.sum,
                summary.cy = out.sum2,
                params = list("modelSep" = modelSep,
                              "lambdaK" = lambdaK,
                              "B" = B,
                              "seed" = seed,
                              "fpv" = fpv,
                              "fpvfutureboot" = fpvfutureboot,
                              "increm.tri" = increm.tri,
                              "triangboot" = triangboot,
                              "lambdafut" = lambdafut))
  }
  class(res) <- "sepprov"
  return(res)
}


