# Utilities file ====

# Check input data -----
#' @name checkData
#' @noRd
#' @keywords internal
check_Data <- function(data){
  if (ncol(data) != nrow(data)){ stop("Matrix has not symmetrical dimensions.") }
  if (rowSums(apply(data, 1, is.na)) != c(0:(ncol(data)-1)) ||
      rowSums(data == 0, na.rm = TRUE) != c(0:(ncol(data)-1))){
    stop("Matrix is not filled with 0 or NAs")
  }
  if (any(is.na(data[(col(data) +  row(data) %in% c(2:(ncol(data) + 1)))]))){
    stop ("Data has missing values")
  }
}

# Compute lower triangle prediction ----
#' @name results
#' @noRd
#' @keywords internal
results <- function(triang, glm.res){
  t <- nrow(triang)
  cij <- as.vector(t(triang))[!is.na(triang)]
  payments <- rep(0, length(cij))
  a <- 1
  newtriang <- triang
  for (i in 2:t) {
    for (j in (t-i+2):t) {
      oy.new <- i; oy.new <- as.factor(oy.new)
      dy.new <- j; dy.new <- as.factor(dy.new)
      newdata <- data.frame(oy = oy.new, dy = dy.new)
      payments[a] <- stats::predict(glm.res, newdata, type = "response")
      a <- a + 1
    }
  }

  a <- 1
  for (i in 2:t) {
    for (j in (t-i+2):t) {
      newtriang[i, j] <- payments[a]
      a <- a + 1
    }
  }
  res <- list(newtriang = newtriang,
              payments = payments)
  return(res)
}


# Create new triangle with missing data ----
#' @name newNA
#' @noRd
#' @keywords internal
newNA <- function(triangle){
  a <- 1
  for(i in 3:ncol(triangle)){
    for(j in 1:a){
      triangle[i, (ncol(triangle)-a+1):ncol(triangle)] <- NA
    }
    a <- a + 1
  }
  rm(a)
  return(triangle)
}

# Create new triangle ----
#' @name NewTriangle
#' @description Function to create a new triangle with the m+1 diagonal data available
#' @noRd
#' @keywords internal
NewTriangle <- function(origin, payments){
  #payments triangle inferior, origin triangle amb dades originals
  R <- dim(payments)[3]
  Ntriangle <- array(NA, dim = c(nrow(origin), ncol(origin), R))
  Latest <- matrix(NA, R, ncol(origin)-1)
  for(z in 1:R){
    aux <- NULL
    Ntriangle[ , , z] <- origin
    for(i in 2:nrow(origin)){
      Ntriangle[i, ncol(origin)-i+2, z] <- payments[i, ncol(origin)-i+2, z]
      Latest[z, i-1] <- payments[i, ncol(origin)-i+2, z]
    }
  }
  res <- list(Ntriangle = Ntriangle,
              Latest = Latest)
  return(res)
}

# Create incremental data triangle ----
#' @name Increm.triangle
#' @description Function to create a incremental data triangle
#' @noRd
#' @keywords internal
Increm.triangle <- function(triangle){
  R <- dim(triangle)[3]
  res <- array(NA, c(nrow(triangle), ncol(triangle), R))
  for(z in 1:R){
    for(i in 1:nrow(triangle)){
      for(j in 1:ncol(triangle)){
        if(!is.na(triangle[i, j, z])){
          res[i, j, z]<-sum(triangle[i, 1:j, z])
        }
      }
    }
  }
  return(res)
}

# Create lower data triangle ----
#' @name TriangInf
#' @description Function to create a lower data triangle
#' @noRd
#' @keywords internal
TriangInf <- function(triangle){
  a <- ncol(triangle)
  for(i in 1:ncol(triangle)){
    for(j in 1:a){
      triangle[i, j] <- NA
    }
    a <- a - 1
  }
  rm(a)
  return(triangle)
}

# Get diagonal ----
#' @name ObtainMDiagonal
#' @description Function to extract diagonal from triangle
#' @noRd
#' @keywords internal
ObtainMDiagonal <- function(triangle){
  aux <- NULL
  for(i in 1:nrow(triangle)){
    aux <- c(aux, triangle[i, (ncol(triangle)-i+1)])
  }
  return(aux)
}

# Create complete triangle ----
#' @name complete.tri
#' @description Function to complete lower data triangle
#' @noRd
#' @keywords internal
complete.tri <- function(triangle, g){
  for(z in 1:dim(triangle)[3]){
    for (i in 2:ncol(triangle)){
      for (j in 1:ncol(triangle)){
        if(is.na(triangle[i, j, z])){
          triangle[i, j, z] <- triangle[i, j-1, z] * g[z, j-1]
        }
      }
    }
  }
  return(triangle)
}


# Summary separateProvision ----
#' @noRd
#' @description Displays a useful description of a glmProvision object
#' @param object glmPprovision object to summarise
#' @param output output style: console or latex
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @method summary sepprov
#' @export
#' @examples
#' res <- separateProvision(TaylorData$lossData, TaylorData$freqData)
#' summary(res)
#' summary(res, output = "latex")

devtools::use_package("xtable")

summary.sepprov  <- function(object, output = "console", ...){
  if (output == "console"){
    cat("==== Summary of separateProvision ====\n\n")
    cat(paste0("modelSep = ", object$params$modelSep,"\n"))
    cat(paste0("lambdaK = ", object$params$lambdaK,"\n"))
    cat(paste0("B = ", object$params$B,"\n"))
    if (! is.null(object$params$seed)){
      cat(paste0("Seed = ", object$params$seed,"\n"))
    } else {
      cat(paste0("No seed used.\n"))
    }
    cat("\n")
    cat("== Summary table per origin year == \n")
    print(object$summary.oy)
    cat("\n")
    cat("== Summary table per calendar year == \n")
    print(object$summary.cy)
    cat("\n")
  }
  if (output == "latex"){
    print(xtable::xtable(object$summary.oy,
                   caption = paste0("Summary statistics per origin year of ",
                                    object$params$modelSep,
                                    " modeling with $\\lambda$ = ",
                                    object$params$lambdaK),
                   digits = 2))
    print(xtable::xtable(object$summary.cy,
                   caption = paste0("Summary statistics per calendar year of ",
                                    object$params$modelSep,
                                    " modeling with $\\lambda$ = ",
                                    object$params$lambdaK, "."),
                   digits = 2))
  }
}

# Summary glmProvision ----
#' @noRd
#' @description Displays a useful description of a glmProvision object
#' @param object glmPprovision object to summarise
#' @param output output style: console or latex
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @method summary glmprov
#' @export
#' @examples
#' res <- glmProvision(TaylorData$lossData)
#' summary(res)
#' summary(res, output = "latex")

devtools::use_package("xtable")

summary.glmprov <- function(object, output = "console", ...){
  fam <- unname(ifelse(object$params$fam == "0", "Normal",
                       ifelse(object$params$fam == "1", "Poisson",
                              ifelse(object$params$fam == "2", "Gamma",
                                     ifelse(object$params$fam == "3",
                                            "Inverse gaussian",
                                            paste0("variance power ", object$params$fam))))))
  link <- ifelse(object$params$link == 0, "logarithmic", object$params$link)
  if (output == "console"){
    cat("==== Summary of glmProvision ====\n\n")
    cat(paste0("peMethod = ", object$params$method,"\n"))
    cat(paste0("GLM: family = ", fam, " and link = ", link, "\n"))
    if (object$params$method == "bootstrap"){
      cat(paste0("B = ", object$params$B,"\n"))
      if (! is.null(object$params$seed)){
        cat(paste0("Seed = ", object$params$seed,"\n"))
      } else {
        cat(paste0("No seed used.\n"))
      }
    }
    cat("\n")
    cat("== Summary table per origin year == \n")
    print(object$summary.oy)
    cat("\n")
    cat("== Summary table per calendar year == \n")
    print(object$summary.cy)
    cat("\n")
  }
  if (output == "latex"){
    print(xtable::xtable(object$summary.oy,
                   caption = paste0("Summary statistics per origin year of ",
                                    object$params[2],
                                    " GLM modeling with ", fam,
                                    " distribution and ", link, " link function."),
                   digits = 2))
    print(xtable::xtable(object$summary.cy,
                   caption = paste0("Summary statistics per calendar year of ",
                                    object$params[2],
                                    " GLM modeling with ", fam,
                                    " distribution and ", link, " link function."),
                   digits = 2))
  }
}

# Summary PDR ----
#' @noRd
#' @description Displays a useful description of a PDR object
#' @param object glmPprovision object to summarise
#' @param output output style: console or latex
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @method summary pdr
#' @examples
#' res <- glmProvision(TaylorData$lossData)
#' res <- PDR(res)
#' summary(res)
devtools::use_package("xtable")

summary.pdr <- function(object, output = "console", ...){
  if (output == "console"){
    cat("==== Summary of PDR ====\n\n")
    print(object$summary)
    cat("\n")
  }
  if (output == "latex"){
    print(xtable::xtable(object$summary,
                         caption = paste0("Summary statistics for PDR."),
                         digits = 2))
  }
}

# Plot separateProvision ----
#' @noRd
#' @rdname plot.sepprov
#' @description Displays plots of a separateProvision object
#' @param object separateProvision object to plot
#' @param log Apply logarithm to latest claims cost values
#' @param plot Plot the final results or retrieve the ggplot object
#' @param ... other arguments ignored (for compatibility with generic)
#' @keywords internal
#' @method plot sepprov
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom reshape2 melt
#' @import ggplot2
#' @include utilities.R

devtools::use_package("ggplot2")

plot.sepprov <- function(x, which.plot = 1:4,
                         log.latest = FALSE, onepage = FALSE, ...){
  if (!is.numeric(which.plot) || any(which.plot < 1) || any(which.plot > 4)){
    stop("The value of which.plot must be in 1:4\n")
  }
  object <- x
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
  q1 <- ggplot(data = data.frame(x = sum.test), aes(x)) +
    geom_histogram(bins = 25, fill = "lightgrey") +
    labs(title = "Histogram of total claims cost",
         x = "Total claims cost", y = "Frequency") +
    theme_bw() +
    theme(panel.grid = element_blank())
  rm(test, sum.test, boots)

  # Boxplot - Latest - OY
  latest <- matrix(0, ncol = nrow(object$params$increm.tri),
                   nrow = dim(object$params$increm.tri)[3])
  for (i in 1:dim(object$params$increm.tri)[3]){
    latest[i, ] <- ObtainMDiagonal(object$params$triangboot[, , i])
  }
  aux <- reshape2::melt(latest)
  colnames(aux) <- c("B", "oy", "latest")
  aux$oy <- aux$oy-1
  aux$oy <- factor(aux$oy)

  aux.pt <- data.frame(Value = ObtainMDiagonal(object$triangle),
                       Index = c(1:(ncol(object$triangle))),
                       stringsAsFactors = FALSE)
  if (log.latest == TRUE){
    aux$latest <- log(aux$latest)
    aux.pt$Value  <- log(aux.pt$Value)
  }
  ylims <- c(min(aux.pt$Value, aux$latest), max(aux.pt$Value, aux$latest))
  q2 <- ggplot(data = aux, aes(x = oy, y = latest)) +
    geom_boxplot(color = "lightgrey") +
    ylim(ylims) +
    geom_point(data = aux.pt, mapping = aes(x = Index, y = Value),
               inherit.aes = FALSE,
               color = "red") +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (log.latest == TRUE){
    q2 <- q2 +
      labs(title = "Latest actual incremental claim against simulated values",
           x = "Origin year", y = "log(Latest incremental claims)")
  } else {
    q2 <- q2 +
      labs(title = "Latest actual incremental claim against simulated values",
           x = "Origin year", y = "Latest incremental claims")
  }
  rm(latest, i, aux, aux.pt, ylims)

  # Boxplot - Ultimate
  ultimate <- matrix(0, ncol = nrow(object$params$increm.tri),
                     nrow = dim(object$params$increm.tri)[3])
  for (i in 1:dim(object$params$increm.tri)[3]){
    ultimate[i, ] <- object$params$increm.tri[, ncol(object$params$increm.tri), i]
  }
  aux <- reshape2::melt(ultimate)
  colnames(aux) <- c("B", "oy", "ultimate")
  aux$oy <- aux$oy - 1
  aux$oy <- factor(aux$oy)
  auxoy <- aux %>% dplyr::group_by(oy) %>% dplyr::summarise(mu = mean(ultimate))

  q3 <- ggplot(data = aux, aes(x = oy, y = ultimate)) +
    geom_boxplot(color = "lightgrey") +
    geom_point(data = auxoy, mapping = aes(x = oy, y = mu),
               inherit.aes = FALSE,
               color = "red") +
    labs(title = "Simulated ultimate claims cost",
         x = "Origin year", y = "Ultimate claims costs") +
    theme_bw() +
    theme(panel.grid = element_blank())
  rm(ultimate, i, aux, auxoy)

  # Boxplot - FPV - CY
  aux <- object$params$fpvfutureboot
  colnames(aux) <- c(as.numeric(1:ncol(aux) + ncol(object$params$increm.tri)))
  aux <- reshape2::melt(aux)
  aux.pt <- data.frame(Value = object$params$fpv,
                       Index = c(1:ncol(object$params$fpvfutureboot) +
                                   ncol(object$params$increm.tri)),
                       stringsAsFactors = FALSE)
  q4 <- ggplot(data = aux, aes(x = as.factor(Var2), y = value)) +
    geom_boxplot(color = "lightgrey") +
    geom_point(data = aux.pt, mapping = aes(x = as.factor(Index), y = Value),
               inherit.aes = FALSE,
               color = "red") +
    labs(title = "Actual future payments per calendar year against simulated values",
         x = "Calendar year", y = "Future payments") +
    theme_bw() +
    theme(panel.grid = element_blank())
  rm(aux, aux.pt)
  if (onepage == TRUE){
    res <- gridExtra::arrangeGrob(q1, q2, q3, q4, nrow = 2, ncol = 2)
    grid::grid.draw(res)
  } else {
    if (any(which.plot == 1)) { print(q1) }
    if (any(which.plot == 2)) { print(q2) }
    if (any(which.plot == 3)) { print(q3) }
    if (any(which.plot == 4)) { print(q4) }
  }
}

# SetClasses ----
#' @noRd
#' @description Set classes to perform custom plotting and summary.
#' @keywords internal
#' @import methods

devtools::use_package("methods")
# Classes
"glmprov" <- methods::setClass("glmprov")
"sepprov" <- methods::setClass("sepprov")
"pdr" <- methods::setClass("pdr")

# Summary method
methods::setMethod(f = "summary", signature = "glmprov",
                   definition = summary.glmprov)
methods::setMethod(f = "summary", signature = "sepprov",
                   definition = summary.sepprov)
methods::setMethod(f = "summary", signature = "pdr",
                   definition = summary.pdr)
# Plot method
methods::setMethod(f = "plot", signature = "sepprov",
                   definition = plot.sepprov)
