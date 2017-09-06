\name{expWeight}
\alias{expWeight}
\title{
  This function computes the smooth windowing weights
}
\description{
  This function computes the exponential weights (smooth windowing function) for different shapes (k)  and bandwidth (l) and provides some plots
}
\usage{
  expWeight(t, k, l, m = 0, plot = FALSE, ...)
}
\arguments{
  \item{t}{
    Vector. a vector of positive continuous values 
  }
  \item{k}{
    A single value for sharpness
  }
  \item{l}{
    A single value for bandwidth
  }
  \item{m}{
    Vector. The location of the modes on t
  }
  \item{plot}{
    Logical flag. Set to true shows a plot of the weights.
  }
  \item{\dots}{
    The parameters that are passed to the plot() function
  }
}
\value{
  A vector of weights
}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

\seealso{
  \code{\link{SmoothWin}}
}
\examples{
  weight = expWeight(
    t = 1:100,
    k = 1,
    l = 5,
    m = c(25, 50, 75),
    plot = TRUE
  )
}

