#' @name PDR
#' @title Claim development results.
#' @description Calculate claim development results.
#' @usage PDR(object, probs = c(0.75, 0.95, 0.995))
#' @param object Output from \code{separateProvision} or \code{glmProvision} functions.
#' @param probs Quantiles of the CDR distribution.
#' @return A list of 6 elements with:
#' \itemize{
#'   \item \code{meanCDRi} CDR per origin year.
#'   \item \code{sdCDRi} Standard deviation of CDR per origin year.
#'   \item \code{meanCDRT} CDR.
#'   \item \code{sdCDRT} Standard deviation of CDR.
#'   \item \code{qCDRTi} Quantiles of CDR distribution.
#'   \item \code{summary} Summary table.
#' }
#' @examples
#' data("TaylorData")
#' res <- separateProvision(TaylorData$lossData, TaylorData$freqData)
#' PDR(res)
#' @export
#' @importFrom stats quantile
#' @include utilities.R

PDR <- function(object, probs = c(0.75, 0.95, 0.995)){
  # 1. Check input ----
  if (!(class(object) %in% c("glmprov", "sepprov"))){
    stop(paste0("Object is not an output from either 'glmProvision' or",
                " 'separateProvision' function."))
  }
  for(i in 1:length(probs)){
    if(!(probs[i] >= 0 & probs[i] <= 1)){
      stop("Probabilities must be between 0 and 1")
    }
  }
  if (class(object) == "glmprov"){
    if (is.null(object$bootstrap.losses)){
      stop("The input is not a bootstrap object")
    }
  }

  origin <- object$triangle
  if (class(object) == "sepprov"){
    freqData <- object$freqData
  }

  if (class(object) == "glmprov"){
    payments <- object$bootstrap.losses
  } else {
    payments <- object$glm.triangle.bootstrap
  }

  # 2. Compute new triangle, (I+1) and the incremental one ----
  trianglem1 <- NewTriangle(origin, payments)
  trianglem1Inc <- Increm.triangle(trianglem1$Ntriangle)

  # 3. Compute the new triangles using the new information ----
  newtriangles <- array(NA, dim(trianglem1$Ntriangle))

  if (class(object) == "glmprov"){
    # 3.1. glmProvision ----
    for (z in 1:dim(trianglem1$Ntriangle)[3]){
      lossData <- trianglem1$Ntriangle[, , z]
      t <- nrow(lossData)
      oy <- c(rep(c(1,2), each = t), rep(3:t, (t-1):2)); oy <- as.factor(oy)
      dy <- c(rep(1:t,2), sequence((t-1):2)); dy <- as.factor(dy)
      v.lossData <- as.vector(t(lossData))
      cij <- v.lossData[!is.na(v.lossData)]

      # Incremental losses matrix
      glmPro <- glm(cij ~ oy + dy,
                    family = statmod::tweedie(as.numeric(object$params$fam),
                                              as.numeric(object$params$link)))
      # Results function to compute the missing values
      triangglm <- results(triang = lossData, glm.res = glmPro)$newtriang
      newtriangles[, , z]<-triangglm
    }
  }

  if (class(object) == "sepprov"){
    # 3.2. separateProvision ----
    modelSep <- object$params$modelSep

    for (z in 1:dim(trianglem1$Ntriangle)[3]){
      if (modelSep == "arithmetic"){
        # 3.2.1 Method: arithmetic ----
        lossData <- trianglem1$Ntriangle[, , z]
        t <- nrow(lossData)
        oy <- c(rep(c(1,2),each=t),rep(3:t, (t-1):2)); oy <- as.factor(oy)
        dy <- c(rep(1:t,2),sequence((t-1):2)); dy <- as.factor(dy)
        cy <- as.numeric(oy) + as.numeric(dy); cy <- as.factor(cy)
        v.lossData <- as.vector(t(trianglem1$Ntriangle[, , z]))
        v.lossData <- v.lossData[!is.na(v.lossData)]

        freq.data.oy <- c(rep(freqData[c(1,2)], each = t), rep(freqData[3:t], (t-1):2))

        average.loss.data <- v.lossData/freq.data.oy

        glmArit <- glm(average.loss.data ~ dy + cy,
                       family = statmod::tweedie(var.power = 1, link.power = 0))

        # Triangle of incremental losses fitted
        m.average.fitted.loss.data <- matrix(0, t, t)
        aux <- c(2*t, (2*t + cumsum((t-1):2)))
        glmArit$fitted.values <- glmArit$fitted.values[-aux]
        rm(aux)

        oy <- rep(1:t, t:1)
        dy <- sequence(t:1)
        for (i in 1:t){
          m.average.fitted.loss.data[i, ] <- c(glmArit$fitted.values[oy == i],
                                               rep(times = i-1, NA))
        }
        m.fitted.loss.data <- m.average.fitted.loss.data * freqData
        triangglm <- matrix(0, t, t)
        for (i in 2:t){
          aux <- glmArit$fitted.values[((oy == (t-i+1)) & (dy == i))] *
            freqData[(t-i+2):t] * (1 + object$params$lambdaK)^(1:(i-1))
          triangglm[, i] <- c(lossData[1:(t-i+1), i], aux)
          rm(aux)
        }
        triangglm[, 1] <- lossData[, 1]
      }

      if (modelSep == "geometric"){
        # 3.2.2 Method: geometric ----
        lossData <- trianglem1$Ntriangle[, , z]
        t <- nrow(lossData)
        oy <- c(rep(c(1,2), each = t),rep(3:t, (t-1):2)); oy <- as.factor(oy)
        dy <- c(rep(1:t,2), sequence((t-1):2)); dy <- as.factor(dy)
        cy <- as.numeric(oy) + as.numeric(dy); cy <- as.factor(cy)
        v.lossData <- as.vector(t(trianglem1$Ntriangle[, , z]))
        v.lossData <- v.lossData[!is.na(v.lossData)]
        freq.data.oy <- c(rep(freqData[c(1,2)], each = t), rep(freqData[3:t], (t-1):2))

        average.loss.data <- v.lossData/freq.data.oy

        for(i in 1:length(average.loss.data)){
          if(average.loss.data[i] == 0){
            average.loss.data[i] <- 1
          }
        }

        glmGeom <- glm(log(average.loss.data) ~ dy + cy,
                       family = statmod::tweedie(var.power = 0,link.power = 1))

        # Triangle of incremental losses fitted
        m.average.fitted.loss.data <- matrix(0, t, t)
        aux <- c(2*t, (2*t + cumsum((t-1):2)))
        glmGeom$fitted.values <- glmGeom$fitted.values[-aux]
        rm(aux)

        oy <- rep(1:t,t:1)
        dy <- sequence(t:1)
        for (i in 1:t){
          m.average.fitted.loss.data[i, ] <- c(glmGeom$fitted.values[oy == i],
                                               rep(times = i-1, NA))
        }
        m.fitted.loss.data <- exp(m.average.fitted.loss.data) * freqData
        triangglm <- matrix(0, t, t)
        for (i in 2:t){
          aux <- exp(glmGeom$fitted.values[((oy == (t-i+1)) & (dy == i))])*
            freqData[(t-i+2):t] * (1 + object$params$lambdaK)^(1:(i-1))
          triangglm[, i] <- c(lossData[1:(t-i+1), i], aux)
          rm(aux)
        }
        triangglm[, 1] <- lossData[, 1]
      }
      newtriangles[, , z] <- triangglm
    }
  }
  # 4. Create cumulative triangles ----
  if (class(object) == "glmprov"){
    CIij <- Increm.triangle(object$bootstrap.losses)
  }
  if (class(object) == "sepprov"){
    aux <- object$glm.triangle
    aux <- array(aux, dim = c(nrow(aux), ncol(aux), dim(trianglem1$Ntriangle)[3]))
    CIij <- Increm.triangle(aux)
    rm(aux)
  }
  CI1ij <- Increm.triangle(newtriangles)

  # 5. Calculate CDRi ----
  CDRi <- matrix(NA, dim(CIij)[3], ncol(CIij)-1)
  for(z in 1:dim(CIij)[3]){
    for(i in 2:nrow(CIij)){
      CDRi[z, i-1] <- CIij[i, ncol(CIij), z] - CI1ij[i, ncol(CIij), z]
    }
  }

  # 6. Calculate the CDR for each origin year and the total CDR ----
  meanCDRi <- c(0, apply(CDRi, 2, mean))
  sdCDRi <- c(0, apply(CDRi, 2, sd))
  qCDRi <- apply(CDRi, 2, quantile, probs, na.rm = TRUE)
  CDRTi <- c(0, rowSums(CDRi, na.rm = TRUE))
  meanCDRT <- mean(CDRTi)
  sdCDRT <- sd(CDRTi)
  qCDRTi <- quantile(CDRTi, probs, na.rm = TRUE)

  # 7. Results object ----
  aux <- matrix(NA, ncol = 5 + length(probs), nrow = nrow(origin)+1,
                dimnames = list(c(rownames(origin), "TOTAL"),
                                c("IBNR", "IBNR.mean", "PredErr.Abs",
                                  "CDR", "CDR.sd",
                                  paste0("CDR.q_", probs))))
  aux[, 1:3] <- object$summary.oy[, 4:6]
  aux[, 4] <- c(meanCDRi, meanCDRT)
  aux[, 5] <- c(sdCDRi, sdCDRT)
  q.aux <- rbind(0, t(qCDRi), t(qCDRTi))
  aux[, 6:ncol(aux)] <- q.aux
  rm(q.aux)

  res <- list(meanCDRi = meanCDRi,
              sdCDRi = sdCDRi,
              meanCDRT = meanCDRT,
              sdCDRT = sdCDRT,
              qCDRTi = qCDRTi,
              summary = aux)
  class(res) <- "pdr"
  return(res)
}
