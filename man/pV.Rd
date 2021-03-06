% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pV.R
\name{pV}
\alias{pV}
\title{Present value of the future payments vector.}
\usage{
pV(object, yieldCurve)
}
\arguments{
\item{object}{Output from \code{separateProvision} or \code{glmProvision} functions.}

\item{yieldCurve}{Vector with the interest rate for the following years.}
}
\value{
A list of 4 elements with:
\itemize{
  \item \code{actualpayments} Present value of the payment per origin year.
  \item \code{totactpay} Present value of the total payment.
  \item \code{yieldCurve} Input vector with interest rates.
  \item \code{fpv} Input vector with the future payments per origin year.
}
}
\description{
Calculate the present value of the future payments vector.
}
