% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNAfingerprint.R
\name{CNAfingerprint}
\alias{CNAfingerprint}
\title{CNAfingerprint prediction function}
\usage{
CNAfingerprint(cnf, target = "OXA")
}
\arguments{
\item{cnf}{CNA features}

\item{target}{Model Prediction Targets, only support "Oxa" now}
}
\value{
CNAfingerprint prediction result
}
\description{
CNAfingerprint prediction function
}
\examples{
\donttest{
features <- readRDS(system.file("extdata", "examplefeatures.rds",package = "CNAfingerprint", mustWork = TRUE))
score <- CNAfingerprint(features,target="Oxa")
}
}
