% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emultir.R
\name{multir.ib.est}
\alias{multir.ib.est}
\title{The eMultiR estimate with information borrow}
\usage{
multir.ib.est(Prop_multir, Prop_scores, trt_ind, main_outcome)
}
\arguments{
\item{Prop_multir}{The calibration score}

\item{Prop_scores}{The information borrow score}

\item{trt_ind}{Treatment assignment (eg., 0,1)}

\item{main_outcome}{The main outcome in binary}
}
\value{
A vector of estimates
}
\description{
This function can be used to obtain the eMultiR estimate with information borrow.
}
