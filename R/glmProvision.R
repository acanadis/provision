#' @name glmProvision
#' @title Calculate provisions using glm modeling.
#' @description Calculate provision using glm modeling.
#' @usage glmProvision(lossData, peMethod = "formula", fam = 1, link = 0, B = 1000, seed = NULL)
#' @param lossData Matrix of incremental losses \eqn{Cij},
#' for \eqn{i = 1,...,t} origin years (rows) and for \eqn{j = 1,...,t}
#' development years (columns); filled with \code{NAs} for \eqn{i + j > t}.
#' @param peMethod Method to be used, can be \code{formula} or \code{bootstrap}. Defaults to "formula".
#' @param fam Index of power variance function as defined in \code{tweedie}. Defaults to Poisson.
#' @param link Index of power link function as defined in \code{tweedie}. Defaults to logarithmic link.
#' @param B Number of iterations to perform in the bootstrapping procedure. Defaults to 1000.
#' @param seed Seed to make bootstrap reproducible.
#' @return A list of 5 elements with:
#' \itemize{
#'   \item \code{triangle} Input data (lossData).
#'   \item \code{glm.triangle} Fitted values by \code{glm} modeling.
#'   \item \code{model} \code{glm} model.
#'   \item \code{bootstrap.losses} (Only with \code{peMethod} = "bootsrap") Three-dimensional
#'    array of lower triangles from bootstrap procedure.
#'   \item \code{summary} Summary table.
#'   \item \code{params} List of output parameters.
#' }
#' @examples
#' data("TaylorData")
#' glmProvision(lossData = TaylorData$lossData)
#' @export
#' @import statmod
#' @import tweedie
#' @importFrom stats glm coef fitted xtabs model.matrix vcov sd
#' @include utilities.R

devtools::use_package("statmod")
devtools::use_package("tweedie")

glmProvision <- function(lossData,
                         peMethod = "formula",
                         fam = 1, link = 0,
                         B = 1000,
                         seed = NULL){
  # 1. Check input ----
  if (! (class(lossData) == "matrix")){
    stop("Object lossData is not a matrix.")
  }
  check_Data(lossData)
  # Function arguments
  if (! (peMethod %in% c("formula", "bootstrap"))){
    stop("Arg peMethod must be either one of 'formula', 'bootstrap'.")
  }
  if (peMethod == "bootstrap"){
    if (! (is.numeric(B) & B > 0)){
      stop(paste0("Selected method is bootstrap but the number of iterations ",
                  "(B) is not a positive integer."))
    }
  }
  # 2. Save common variables ----
  t <- nrow(lossData)
  oy <- rep(1:t, t:1); oy <- as.factor(oy)
  dy <- sequence(t:1); dy <- as.factor(dy)
  v.lossData <- as.vector(t(lossData))
  cij <- v.lossData[!is.na(v.lossData)]

  ## Incremental losses matrix
  glmPro <- glm(cij ~ oy + dy, family = statmod::tweedie(fam, link))
  # Results function to compute the missing values
  triangglm <- results(triang = lossData, glm.res = glmPro)$newtriang
  payments <- results(triang = lossData, glm.res = glmPro)$payments

  ## Origin year reserve
  oyres <- rep(0, t-1)
  for (i in 2:t){
    oyres[i-1] <- sum(triangglm[i, (t-i+2):t])
  }

  ## Total reserve
  totres <- sum(oyres)

  ## Vector of future payments
  fpv <- rep(0, dim(triangglm)[1] - 1)
  for (k in 1:dim(triangglm)[1] - 1){
    future <- row(triangglm) + col(triangglm) - 1 == dim(triangglm)[1] + k
    fpv[k] <- sum(triangglm[future])
  }

  ## Pearson residuals
  ro <- fam
  coefs <- exp(as.numeric(coef(glmPro)))
  alpha <- c(1, coefs[2:t])*coefs[1]
  beta <- c(1,coefs[(t+1):(2*t-1)])
  orig.fits <- alpha %*% t(beta)
  future <- row(orig.fits) + col(orig.fits) - 1 > t
  Prs.resid <- (cij - fitted(glmPro)) / sqrt(fitted(glmPro)^ro)

  ## Computation of n, p and phi.p
  n <- (t * (t + 1)) / 2
  p <- 2 * t - 1
  phi.P <- sum(Prs.resid^2) / (n-p)

  ## Scaled residuals
  Adj.Prs.resid <- Prs.resid * sqrt(n / (n-p))

  ## Variance and covariance matrices
  cij.l <- xtabs(cij ~ oy + dy)
  cij.v <- as.vector(cij.l)
  ii <- row(cij.l); jj <- col(cij.l)
  futurebis <- as.numeric(ii+jj-1 > k)
  ii <- as.factor(ii); jj <- as.factor(jj)
  Cov.eta <-  model.matrix(cij.v ~ ii + jj) %*% vcov(glmPro) %*%
    t(model.matrix(cij.v ~ ii + jj))
  mu.hat <- as.vector(orig.fits * future)

  # 3. Method: formula ----
  if(peMethod == "formula"){ # Prediction error with formula.
    ## Total reserve
    PEfor <- sqrt(phi.P * sum(mu.hat^ro) + t(mu.hat) %*% Cov.eta %*% mu.hat)

    ## Origin year reserve
    mu.hat.m <- matrix(mu.hat, t, t)
    mu.hat.orig <- matrix(0, t, t)
    mu.hat.vec <- numeric(t * t)
    PEfororig <- numeric(t-1)

    for(orig in 1:(t-1)){
      mu.hat.orig <- rbind(matrix(0, orig, t),
                           mu.hat.m[orig+1, ],
                           matrix(0, t-(orig+1), t))
      mu.hat.vec <- as.vector(mu.hat.orig)
      PEfororig[orig] <- sqrt(phi.P*sum(mu.hat.vec^ro) +
                                t(mu.hat.vec) %*% Cov.eta %*% mu.hat.vec)
    }

    ## Calendar year future payments
    mu.hat.m <- matrix(mu.hat, t, t)
    mu.hat.d <- matrix(0, t, t)
    mu.hat.vec <- numeric(t*t)
    PEforcal <- numeric(t-1)
    for(cal in 2:t){
      mu.hat.d[row(mu.hat.m) + col(mu.hat.m) ==
                 t + cal] <- mu.hat.m[row(mu.hat.m) + col(mu.hat.m) == t + cal]
      mu.hat.vec <- as.vector(mu.hat.d)
      PEforcal[cal-1] <- sqrt(phi.P * sum(mu.hat.vec^ro) +
                                t(mu.hat.vec) %*% Cov.eta %*% mu.hat.vec)
      mu.hat.d <- matrix(0,t,t)
    }
  }

  # 4. Method: bootstrap ----
  if(peMethod == "bootstrap"){
    nBoot <- B
    # Total reserve
    if(! is.null(seed)){
      set.seed(seed)
    }
    payments <- reserves <- numeric(nBoot)
    for(boots in 1:nBoot){
      Ps.cij <- sample(Adj.Prs.resid, n, replace = TRUE)
      Ps.cij <- Ps.cij * sqrt(fitted(glmPro)) + fitted(glmPro)
      Ps.cij <- pmax(Ps.cij, 0)
      Ps.CL <- glm(Ps.cij ~ oy + dy , family = statmod::tweedie(fam, link))
      coefs <- exp(as.numeric(coef(Ps.CL)))
      Ps.alpha <-c (1, coefs[2:t]) * coefs[1]
      Ps.beta <- c(1, coefs[(t+1):(2*t-1)])
      Ps.fits <- Ps.alpha %*% t(Ps.beta)
      Ps.reserve <- sum(Ps.fits[future])
      Ps.totpayments <- tweedie::rtweedie(1, mu=Ps.reserve, phi=phi.P, power=fam)
      reserves[boots] <- Ps.reserve
      payments[boots] <- Ps.totpayments
    }
    PEbs <- sqrt(phi.P * sum(orig.fits[future]^ro)+sd(reserves)^2)

    # Statistics
    cv <- (sd(payments) / mean(payments)) * 100 # CV in percentage
    pp <- (payments - mean(payments)) / sd(payments)
    # sum(pp^3) / (nBoot-1)  # Skewness estimation
    # sum(pp^4) / (nBoot-1) - 3 # Kurtosis esitmation
    # hist(payments, breaks = 21, prob = TRUE, main = "Predictive distribution of total reserve")
    # lines(density(payments), lty="dashed")
    # curve(dnorm(x, mean = mean(payments), sd = sd(payments)),
    #       lty = "dotted", add = TRUE)

    ## Origin year reserve
    if(! is.null(seed)){
      set.seed(seed)
    }
    payments <- reservesorig <- matrix(0, nBoot, t-1)
    for (boots in 1:nBoot){
      Ps.cij <- sample(Adj.Prs.resid, n, replace = TRUE)
      Ps.cij <- Ps.cij * sqrt(fitted(glmPro)) + fitted(glmPro)
      Ps.cij <- pmax(Ps.cij, 0)
      Ps.CL <- glm(Ps.cij ~ oy + dy, family = statmod::tweedie(fam,link))
      coefs <- exp(as.numeric(coef(Ps.CL)))
      Ps.alpha <- c(1, coefs[2:t]) * coefs[1]
      Ps.beta <- c(1, coefs[(t+1):(2*t-1)])
      Ps.fits <- Ps.alpha %*% t(Ps.beta)
      provor <- numeric(t-1)
      payori <- numeric(t-1)
      for(orig in 1:(t-1)){
        provor[orig] <- sum(Ps.fits[orig+1, (t-(orig-1)):t])
        payori[orig] <- tweedie::rtweedie(1, mu = provor[orig],
                                          phi = phi.P, power = fam)
      }
      reservesorig[boots, ] <- provor
      payments[boots, ] <- payori
    }

    PEbsorig <- numeric(t-1)
    for (orig in 1:(t-1)){
      PEbsorig[orig] <- sqrt(phi.P * sum((orig.fits[orig+1, (t-(orig-1)):t])^ro) +
                               sd(reservesorig[, orig])^2)
    }

    ## Calendar year future payments
    if(! is.null(seed)){
      set.seed(seed)
    }
    payments <- reservescal <- matrix(0, nBoot, t-1)

    for (boots in 1:nBoot){
      Ps.cij <- sample(Adj.Prs.resid, n, replace = TRUE)
      Ps.cij <- Ps.cij * sqrt(fitted(glmPro)) + fitted(glmPro)
      Ps.cij <- pmax(Ps.cij, 0)
      Ps.CL <- glm(Ps.cij ~ oy + dy, family = statmod::tweedie(fam,link))
      coefs <- exp(as.numeric(coef(Ps.CL)))
      Ps.alpha <- c(1, coefs[2:t]) * coefs[1]
      Ps.beta <- c(1, coefs[(t+1):(2*t-1)])
      Ps.fits <- Ps.alpha %*% t(Ps.beta)

      matres.m <- matrix(Ps.fits, t, t)
      matres.d <- matrix(0, t, t)
      provc <- numeric(t-1)
      paycal <- numeric(t-1)
      for(cal in 2:t){
        matres.d[row(matres.m) + col(matres.m) ==
                   t + cal] <- matres.m[row(matres.m) + col(matres.m) == t + cal]
        provc[cal-1] <- sum(matres.d)
        paycal[cal-1] <- tweedie::rtweedie(1,mu=provc[cal-1],phi=phi.P,power=fam)
        matres.d <- matrix(0,t,t)
      }
      reservescal[boots, ] <- provc
      payments[boots, ] <- paycal
    }

    PEbscal <- numeric(t-1)
    cal.d <- matrix(0, t, t)
    for(cal in 1:t-1){
      cal.d[row(triangglm) + col(triangglm) ==
              t+cal+1] <- triangglm[row(triangglm) + col(triangglm) == t+cal+1]
      PEbscal[cal] <- sqrt(phi.P*sum(cal.d^ro) + sd(reservescal[, cal])^2)
      cal.d <- matrix(0, t, t)
    }

    # Incremental triangles with bootstrap
    if(! is.null(seed)){
      set.seed(seed)
    }
    bootlosses <- array(0, dim = c(t, t, nBoot))
    for (boots in 1:nBoot){
      Ps.cij <- sample(Adj.Prs.resid, n, replace = TRUE)
      Ps.cij <- Ps.cij*sqrt(fitted(glmPro)) + fitted(glmPro)
      Ps.cij <- pmax(Ps.cij, 0)
      Ps.CL <- glm(Ps.cij ~ oy + dy, family = statmod::tweedie(fam,link))
      coefs <- exp(as.numeric(coef(Ps.CL)))
      Ps.alpha <- c(1, coefs[2:t])*coefs[1]
      Ps.beta <- c(1, coefs[(t+1):(2*t-1)])
      Ps.fits <- Ps.alpha %*% t(Ps.beta)
      matres.d <- matrix(0, t, t)
      matres.m <- matrix(Ps.fits, t, t)
      for (i in 1:t){
        for (j in 1:t){
          matres.d[i, j] <- tweedie::rtweedie(1, mu = matres.m[i, j],
                                              phi = phi.P, power = fam)
        }
      }
      bootlosses[, , boots] <- matres.d
    }
  }
  # 5. Results object ----
  if (peMethod == "bootstrap"){
    ## Summary creation - BOOTSTRAP
    increm.tri <- Increm.triangle(bootlosses)
    ## Origin year
    labelscy <- paste0("cy", (t+1):(t+ncol(lossData)-1))
    out.sum <- matrix(NA, ncol = 10, nrow = t+1,
                      dimnames = list(c(rownames(lossData), "TOTAL"),
                                      c("Latest.mean", "dev.to.date", "Ultimate.mean",
                                        "IBNR", "IBNR.mean", "PredErr", "CV",
                                        "IBNR.quantile.75", "IBNR.quantile.95",
                                        "IBNR.quantile.99")))

    ## Latest
    latest <- matrix(NA, ncol = nrow(lossData), nrow = B)
    for (i in 1:B){
      latest[i, ] <- c(ObtainMDiagonal(increm.tri[, , i]))
    }
    out.sum[, 1] <- c(colMeans(latest), sum(colMeans(latest)))
    rm(latest, i)
    ## Ultimate
    ultimate <- NULL
    for (i in 1:nrow(bootlosses)){
      ultimate <- c(ultimate, mean(increm.tri[i, ncol(increm.tri), ]))
    }
    out.sum[, 3] <- c(ultimate, sum(ultimate))
    rm(ultimate, i)
    ## dev.to.date
    out.sum[, 2] <- out.sum[, 1]/out.sum[, 3]

    out.sum[, 4] <- c(0, oyres, sum(oyres))
    out.sum[, 5] <- c(0, apply(reservesorig,2,mean), sum(apply(reservesorig,2,mean)))
    out.sum[, 6] <- c(0, abs(PEbsorig), PEbs)
    out.sum[, 7] <- out.sum[,6]/out.sum[,4]
    out.sum[, 8] <- c(0, apply(reservesorig, 2, quantile, 0.75, na.rm = TRUE),
                     quantile(sum(reservesorig), 0.75, na.rm = TRUE))
    out.sum[, 9] <- c(0, apply(reservesorig, 2, quantile, 0.95, na.rm = TRUE),
                      quantile(sum(reservesorig), 0.95, na.rm = TRUE))
    out.sum[, 10] <- c(0, apply(reservesorig, 2, quantile, 0.995, na.rm = TRUE),
                      quantile(sum(reservesorig), 0.995, na.rm = TRUE))
    out.sum[is.nan(out.sum)] <- 0

    ## Calendar year
    labelscy <- paste0("cy", (t+1):(t+ncol(lossData)-1))
    out.sum2 <- matrix(NA, ncol = 7, nrow = t,
                       dimnames = list(c(labelscy, "TOTAL.cy"),
                                       c("IBNR", "IBNR.mean", "PredErr", "CV",
                                         "IBNR.quantile.75", "IBNR.quantile.95",
                                         "IBNR.quantile.99")))
  out.sum2[, 1] <- c(fpv, sum(fpv))
  out.sum2[, 2] <- c(apply(reservescal,2,mean), sum(apply(reservescal,2,mean)))
  out.sum2[, 3] <- c(abs(PEbscal), sum(abs(PEbscal)))
  out.sum2[, 4] <- out.sum2[,3]/out.sum2[,1]
  out.sum2[, 5] <- c(apply(reservescal, 2, quantile, 0.75, na.rm = TRUE),
                     quantile(sum(reservescal), 0.75, na.rm = TRUE))
  out.sum2[, 6] <- c(apply(reservescal, 2, quantile, 0.95, na.rm = TRUE),
                     quantile(sum(reservescal), 0.95, na.rm = TRUE))
  out.sum2[, 7] <- c(apply(reservescal, 2, quantile, 0.995, na.rm = TRUE),
                     quantile(sum(reservescal), 0.995, na.rm = TRUE))

  out.sum2[is.nan(out.sum2)] <- 0
  }

  if (peMethod == "formula"){
    # Summary creation - FORMULA
    ## Origin year
    increm.tri <- Increm.triangle(array(triangglm,
                                        dim = c(ncol(triangglm), ncol(triangglm), 1)))
    out.sum <- matrix(NA, ncol = 6, nrow = t+1,
                      dimnames = list(c(rownames(lossData), "TOTAL"),
                                      c("Latest", "dev.to.date", "Ultimate",
                                        "IBNR", "IBNR.PredErr", "CV")))
    ## Latest
    out.sum[, 1] <- c(ObtainMDiagonal(increm.tri[, ,1]),
                      sum(ObtainMDiagonal(increm.tri[, , 1])))
    ## Ultimate
    out.sum[, 3] <- c(increm.tri[, ncol(increm.tri), 1],
                      sum(increm.tri[, ncol(increm.tri), 1]))
    ## dev.to.date
    out.sum[, 2] <- out.sum[, 1]/out.sum[, 3]
    ## IBNR
    out.sum[, 4] <- c(0, oyres, totres)
    ## PE
    out.sum[, 5] <- c(0, abs(PEfororig), PEfor)
    ## CV
    out.sum[, 6] <- out.sum[, 5]/out.sum[, 4]
    out.sum[is.nan(out.sum)] <- 0

    #Calendar year
    labelscy <- paste0("cy", (t+1):(t+ncol(lossData)-1))
    out.sum2 <- matrix(NA, ncol = 3, nrow = t,
                       dimnames = list(c(labelscy, "TOTAL.cy"),
                                       c("IBNR", "PredErr", "CV")))
    ## IBNR
    out.sum2[, 1] <- c(fpv, sum(fpv))
    ## Pred. Error
    out.sum2[, 2] <- c(abs(PEforcal),sum(abs(PEforcal)))
    ## CV
    out.sum2[, 3] <- out.sum2[, 2]/out.sum2[, 1]
    out.sum2[is.nan(out.sum2)] <- 0
  }
  # Build result object
  if (peMethod == "bootstrap"){
    res <- list(triangle = lossData,
                glm.triangle = triangglm,
                glm.model = glmPro,
                bootstrap.losses = bootlosses,
                summary.oy = out.sum,
                summary.cy = out.sum2,
                params = list("method" = peMethod,
                              "fam" = fam,
                              "link" = link,
                              "B" = B,
                              "seed" =  seed,
                              "fpv" = fpv))
  }
  if (peMethod == "formula"){
    res <- list(triangle = lossData,
                glm.triangle = triangglm,
                glm.model = glmPro,
                summary.oy = out.sum,
                summary.cy = out.sum2,
                params = list("method" = peMethod,
                              "fam" = fam,
                              "link" = link,
                              "fpv" = fpv))
  }
  class(res) <- "glmprov"
  return(res)
}
