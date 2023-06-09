% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanGP.R
\name{meanGP.GibbsZoop}
\alias{meanGP.GibbsZoop}
\title{Posterior Mean Surfaces of the GP's}
\usage{
meanGP.GibbsZoop(
  model,
  coords,
  date = NULL,
  oneY2 = NULL,
  oneT2 = NULL,
  process = c("eta", "alpha", "lambda", "phi"),
  pred.coords,
  background,
  xlim = c(-61771.464, 1228.536),
  ylim = c(-34392.144, 13607.86),
  zlim = NULL
)
}
\arguments{
\item{model}{parameters of \code{\link{GibbsZoop}}}

\item{coords}{\eqn{N \times 2}}

\item{date}{\eqn{N \times 1}}

\item{oneY2}{binary wether \code{Y2} exists or not: \eqn{N \times 1}}

\item{oneT2}{binary wether \code{T2} exists or not: \eqn{N \times 1}}

\item{process}{one of \code{"eta"} (zoop GP's), \code{"alpha"} (temp 
calibration GP's), \code{"lambda"} (zoop calibration GP's), or 
\code{"phi"} (temp GP)}

\item{pred.coords}{coords where to predict}

\item{background}{background of CCB}

\item{xlim, ylim, zlim}{x y z limits in the map}
}
\value{
Saved png
}
\description{
This function saves the map of the posterior mean of the
  GP's in \code{\link{GibbsZoop}} when \code{calibration = "corTime"}
}
\author{
Jorge Castillo-Mateo
}
