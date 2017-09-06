\name{SmoothWin}
\alias{SmoothWin}
\title{
  Implementation of the soft windowing on the linear model
}
\description{
  Implementation of the soft windowing on the linear model by imposing weights on the model. The parameters k and l control the sharpness and bandwidth of the windowing function. 
}
\usage{
  SmoothWin(
    t                                   ,
    x                                   ,
    y                                   ,
    m                                   ,
    l = SmoothWin:::lseq(
      from = 1                          ,
      to = (max(t) - min(t)) / length(m),
      length.out = 51
    )                                   ,
    k = SmoothWin:::lseq(
      from = 1                          ,
      to = 10                           ,
      length.out = 10
    )                                   ,
    min.obs = min(length(t) - length(m)                               , 
                   log(length(m) * length(t)) * length(unique(m)))    ,
    threshold = 10 ^ -8                 ,
    criteria  = 'AICc'                  ,
    method    = 'enet'                  ,
    plot = FALSE
  )
}
\arguments{
  \item{t}{
    Vector. A vector of positive continuous values for time (or any representative of time)
  }
  \item{x}{
    Matrix. The design matrix
  }
  \item{y}{
    Vector. Response 
  }
  \item{m}{
    Vector. The location of the modes on the time vector t. For example, 25 means the 25th element of t
  }
  \item{l}{
    Single value for the bandwidth of the windows
  }
  \item{k}{
    Single value for the sharpness of the windows
  }
  \item{min.obs}{
    Single value. The minimum observations that must be included in the window
  }
  \item{threshold}{
    Single Value. The values of weight below threshold would be pushed to zero.
  }
  \item{criteria}{
    Case sensitive string. Choose from 'AICc', 'BIC', 'Cp' and 'GCV' for the model selection criteria in penalized change point detection algorithm. The default is 'AICc'
  }
  \item{method}{
    Case sensitive string. Choose between 'enet' and 'alasso' for lasso and adaptive-lasso respectively
  }
  \item{plot}{
    Logical flag. Set to TRUE to see some plots
  }
}

\value{
  \item{final.k,final.l }{Final values for k and l}
  \item{finalModel,model.l,model.k }{List of models for l,k and the final model. Each models contains weights, models and details}
  \item{others}{The input parameters such as x,y,t and ... are also included in the output}
}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

\seealso{
  \code{\link{expWeight}}
}
\examples{
  ########### Generating data
  set.seed(123456)
  n       = 270
  t       = 1:n
  x       = sort(rnorm(n))
  y       = x * 0
  
  ind     = (n / 3):(2 * n / 3)
  y[ind] = x[ind] * 1 #        # regression coefficient is 1 in the window
  y       = y + rnorm(n, 0, .1)
  
  ############ Applying model & results
  plot(x, y)                   # plot data
  r = SmoothWin(               # run  windowing algorithm
    t = t,
    x = matrix(x, ncol = 1),
    y = y,
    m = n / 2,
    l = 1:60,
    k = 1:10,
    min.obs = 25,
    plot = TRUE
  )
  abline(                       # draw the true borders 
    v = c(n / 3, 2 * n / 3),
    lty = 3,
    col = 'gray',
    lwd = 2
  )
  legend(                       # legends
    x = 'topright',
    legend = c('Estimated window', 'True window'),
    fill = c('black', 'Gray')
  )
}

