\name{penCPD}
\alias{penCPD}
\title{
  Penalized change point detection
}
\description{
  The implementation of penalized change point detection algorithm on a signal. 
}
\usage{
  penCPD(
    x                  ,
    plot = FALSE       ,
    threshold = 0.01   ,
    method = "enet"    ,
    criteria = "AICc"  ,
    ...
  )
}

\arguments{
  \item{x}{
    Vector. Vector of the data
  }
  \item{plot}{
    Logical flag. Set to true returns some plots
  }
  \item{threshold}{
    Single positive value. Threshold to filter the change points. The effect sizes below the threshold would be pushed to zero
  }
  \item{method}{
    Case sensitive string. Choose between 'enet' and 'alasso' for lasso and adaptive-lasso change point detection. Defaults is 'enet'.
  }
  \item{criteria}{
    Case sensitive string. Choose between 'Cp', 'AICc', 'BIC' and 'GCV' for the model selection criteria in penalized algorithm. The default is 'AICc'. 
  }
  \item{\dots}{
    parameters that are passed to the plot() function
  }
}

\value{
  a vector of change point indices. For example 5 means the fifth value in x
}
\references{
  On change point detection using the fused lasso method, CR Rojas, B Wahlberg - arXiv preprint arXiv:1401.5408, 2014
}
\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

\examples{
  set.seed(123456)
  n   = 100
  x   = c(rnorm(n,0,.1),rnorm(n,1,.1))
  CPs = penCPD(x = x,plot = TRUE,method = 'alasso')
  plot(x)
  abline(v=CPs)
}