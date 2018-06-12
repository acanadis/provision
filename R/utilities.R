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
# !!!! la fem servir?
#' @name newNA
#' @noRd
#' @keywords internal
newNA <- function(triangle){
  a <- 1
  for(i in 3:ncol(triangle)){ # !!!! 3? 2+1? verificar.
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
  link <- ifelse(object$params$link == 0, "logit", object$params$link)
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
#' @export
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

