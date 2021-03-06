% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{equivTestFixedMarginCombPlot}
\alias{equivTestFixedMarginCombPlot}
\title{Provide a combined plot for equivvalence test}
\usage{
equivTestFixedMarginCombPlot(et)
}
\arguments{
\item{et}{the list returned by \code{equivTestFixedMargin}}
}
\description{
Provide a combined plot for equivalence test, including both scatter plot of the sample data and a bar plot indicating the test result, where the null hypothesis is rejected if the red line representing the mean value of the test product lies within a grey rectangle centered at a blue line representing the mean value of the reference product.
}
\examples{
vecR = rnorm(20,0,1)
vecT = rnorm(20,-1.5,1)
et = equivTestFixedMargin(vecT,vecR)
equivTestFixedMarginCombPlot(et)
}
