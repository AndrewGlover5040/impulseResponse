% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InfluenceFunctions.R
\name{Influence}
\alias{Influence}
\title{Influcence Function}
\usage{
Influence(params, startDay, t_p = 0)
}
\arguments{
\item{params}{A numerical vector of length five. It encodes
[c(p_0, k_1, k_2, tau_1, tau_2)].}

\item{startDay}{The day that training begins at}

\item{t_p}{The day that the performance is supposed to happen}
}
\value{
A vector of length t_p-startday
}
\description{
Influcence Function
}
\examples{
lk;j
}
