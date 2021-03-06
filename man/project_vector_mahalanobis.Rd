% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newMLEfunction.R
\name{project_vector_mahalanobis}
\alias{project_vector_mahalanobis}
\title{A function for projecting a vector such that it will satisfy a linear constraint
based on a mahalanobis metric}
\usage{
project_vector_mahalanobis(x, target, selected, mean_weights, mahal_const,
  mahal_vec)
}
\arguments{
\item{x}{the vector to be projected}

\item{target}{what the linear function of x should equal after the projection}

\item{selected}{which coordinates of x were selected}

\item{mean_weights}{we require that \code{sum(x * mean_weights)} will equal \code{target}}
}
\description{
A function for projecting a vector such that it will satisfy a linear constraint
based on a mahalanobis metric
}
