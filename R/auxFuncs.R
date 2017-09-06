.onAttach <- function( lib, pkg ) {
  packageStartupMessage(
    paste0(
      '\n ========================================================================',
      '\n If you have any question about this package contact me on               ',
      '\n      hamedhaseli@gmail.com  or visit www.hamedhaseli.webs.com           ',
      '\n *** This project is supported by European Bioinformatic Institute (EBI) ',
      '\n ========================================================================'
    ),
    domain = NULL,  appendLF = TRUE )
}
###########
expF = function(x, k, l, m) {
  m  = unique(m)
  r = plogis(x,
             location = m - l,
             scale = 1 / k) *
    (1 - (plogis(
      x,
      location = m + l,
      scale = 1 / k
    )))
  return(r)
}
###########
expWeight = function(t, k, l, m = 0, plot = FALSE , ...) {
  m  = unique(m)
  lm = length(m)
  if (lm > 10)
    message('more than 10 modes detected. It can take a long time to finish!')
  r  = matrix(0, ncol = length(m), nrow = length(t))
  for (i in 1:lm) {
    r[, i] = expF (x = t,
                   k = k,
                   l = l,
                   m = m[i])
  }
  s = rowSums(r)
  if (lm > 1) {
    for (j in 2:lm) {
      cm  = combn(lm, j)
      for (i in 1:ncol(cm)) {
        s = s + (-1) ^ (j + 1) * apply(r, 1, function(x) {
          prod(x[cm[, i]])
        })
      }
    }
  }
  if (plot) {
    plot(
      t,
      s,
      xlab = 'Time',
      ylab = 'Weight',
      xaxt = 'n',
      xlim = c(min(t, na.rm = TRUE) * .85, max(t, na.rm = TRUE) * 1.5),
      ...
    )
    abline(
      v = c(m),
      lty = 1,
      col = 1:lm,
      lwd = 2
    )
    abline(v = c(m + l, m - l),
           lty = 2,
           col = 1:lm)
    axis(
      side = 1,
      at = c(m, m - l, m + l),
      labels = c(m, m - l, m + l),
      lwd.ticks = 2
    )
    legend(
      'topright',
      title = 'm/m+l/m-l',
      legend = paste(round(m, 3), round(m + l, 3), round(m - l, 3), sep = '/'),
      col = 1:lm,
      fill = 1:lm
    )
  }
  return(s)
}
###########
penCPD = function(x,
                  plot      = FALSE,
                  threshold = .01,
                  method    = 'enet',
                  criteria  = 'AICc',
                  ...) {
  #requireNamespace('msgps', 'msgps')
  x  = as.vector(x)
  lx = length(x)
  m  = matrix(1, nrow = lx, ncol = lx)
  for (i in 2:(lx)) {
    m[, i] = c(rep(0, i - 1), rep(1, (lx - i + 1)))
  }
  r = msgps::msgps(
    X = m,
    y = x,
    intercept = FALSE,
    penalty = method
  )
  if (plot)
    plot(r, ...)
  nzc = which(abs(coef(r)[, criteria]) > threshold)
  return(nzc)
}
###########
scale21 <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
    message('max and min are the same!')
    return(x)
  } else{
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
}
###########
gridSearch = function(t,
                      x                 ,
                      y                 ,
                      m = mean(y)       ,
                      l = 1             ,
                      k = 1             ,
                      plot = TRUE       ,
                      mutInd = NULL     ,
                      threshold = 10 ^ -18,
                      ...) {
  lk = length(k)
  ll = length(l)
  n  = length(y)
  m  = unique(m)
  lmodel = list()
  wmat   = matrix(0, ncol = n + 1   , nrow = ll * lk) # 1 for index
  rmat   = matrix(0, ncol = 8 + ncol(x) , nrow = ll * lk)
  colnames(rmat)          = colnames(rmat, do.NULL = FALSE)
  colnames(rmat)[1:7]     = c('Ind',
                              'Obs.in.Interval',
                              'AIC',
                              'BIC',
                              'residual.SD',
                              'k',
                              'l')
  if (!is.null(colnames(x))) {
    colnames(rmat)[-c(1:7)] = c('Int', colnames(x))
  } else{
    colnames(rmat)[-c(1:7)] = paste('x.', 0:ncol(x), sep = '')
    colnames(rmat)            = colnames(rmat, do.NULL = FALSE, prefix = 'x.')
  }
  
  counter = 1
  for (lp in l) {
    for (kp in k) {
      weight = expWeight(
        t = t,
        k = kp,
        l = lp,
        m = m,
        plot = 0
      )
      inn = abs(weight) > threshold
      if (plot) {
        plot(
          t,
          y,
          col = inn + 1,
          pch = inn + 1,
          main = paste(
            'l=',
            round(lp, 2),
            ', k=',
            round(kp, 2),
            ', #=',
            sum(inn),
            sep = ''
          ),
          xlab = 'Time',
          ylab = 'Response',
          ...
        )
        if (!is.null(mutInd))
          points(t[mutInd], y[mutInd], col = 3, pch = 19)
        lines(t,
              scale21(weight) * (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)) + min(y, na.rm = TRUE),
              lty = 3)
        abline(v = m)
      }
      slm = lm(
        y ~ . ,
        data = data.frame(y = y, x = x),
        weights  = weight / sum(weight) * (abs(weight) > threshold)
      )
      lmodel[[counter]]   = slm
      rmat  [counter, ]   = c(counter ,
                              sum(inn),
                              AIC(slm),
                              BIC(slm),
                              sd(resid(slm)),
                              kp,
                              lp,
                              coef(slm))
      wmat[counter, ]    = c(counter, weight / sum(weight))
      cat('\r', counter, '|', lk * ll)
      counter           = counter  + 1
    }
  }
  return(list(
    weights = wmat,
    output = as.data.frame(rmat),
    models = lmodel
  ))
}
###########
lseq = function (from = 1,
                 to = 5,
                 length.out = 6,
                 adj = 1)
{
  r = exp(seq(log(from) / adj, log(to), length.out = length.out))
  return(r)
}