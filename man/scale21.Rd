\name{scale21}
\alias{scale21}
\title{
  Scale data to site in (0,1) interval
}
\description{
  Scale data to site in (0,1) interval
}
\usage{
scale21(x)
}
\arguments{
  \item{x}{
  A vector of values
}
}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

\examples{
    set.seed(123456)
    x   = rnorm(100, 0, 2)
    x21 = scale21(x)
    matplot(cbind(x, x21), type = 'b', pch = 1:2)
    abline(h = c(0, 1), lty = 3)
    legend('topright',
           legend = c('Original', 'Scaled'),
           fill = 1:2)
}

