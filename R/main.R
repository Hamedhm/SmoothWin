SmoothWin = function(t                                                ,
                     x                                                ,
                     y                                                ,
                     m                                                ,
                     l = SmoothWin:::lseq(
                       from = 1                    ,
                       to = (max(t) - min(t)) / length(m)  ,
                       length.out = 51
                     )                    ,
                     k = SmoothWin:::lseq(from = 1                    ,
                                          to = 10                     ,
                                          length.out = 10)            ,
                     min.obs   = min(length(t) - length(m), 
                                     log(length(m) * length(t)) * length(unique(m)))        ,
                     threshold = 10 ^ -8                              ,
                     criteria  = 'AICc'                               ,
                     method    = 'enet'                               ,
                     plot      = FALSE) {
  min.obs = ceiling(min.obs)
  k = sort(k[k > 0 & !is.na(k)], decreasing = TRUE)
  l = sort(l[l > 0 & !is.na(l)])
  
  if (length(k) <= 1 | length(l) <= 1)
    stop('l and k must have at least 2 values length!')
  
  if (length(unique(t[m])) > 10)
    message('More than 10 modes detected. The entire procedure can take a long time!')
  
  if (length(m) > min.obs) {
    stop('min.obs is less than the number of treatments!')
  } else{
    cat(
      '\n ** Minimum observations in the model : (#Treatment =',
      length(m),
      ', #Min observations =',
      min.obs,
      ') =',
      min.obs + length(m)
    )
  }
  ### 1. Determining l
  cat('\n 1|3 Searching for l ...\n')
  rl = gridSearch(
    t = t,
    x = x,
    y = y,
    m = t[m],
    l = l,
    k = max(k),
    mutInd    = m,
    plot      = 0,
    threshold = threshold
  )
  scp1        = scale21(rl$output$residual.SD) * 100
  if(length(unique(scp1))>1){
    cpt         = penCPD(
      scp1,
      threshold = threshold,
      criteria = criteria,
      method = method,
      plot  = plot,
      main = 'l - CPD solution path'
    )
  }else{
    cpt  = length(scp1)
  }
  finall      = rl$output$l[cpt][which(rl$output$Obs.in.Interval[cpt] >=
                                         (min.obs + length(m)))][1]
  if (is.na(finall)) {
    cat('\n An optimal l is not found. Max l will be used.')
    finall = max(l)
  }
  
  ### 2. Determining k
  cat('\n 2|3 Searching for k ... \n')
  rk    = gridSearch(
    t = t,
    x = x,
    y = y,
    m = t[m],
    l = finall,
    k = k,
    plot = 0,
    mutInd = m,
    threshold = threshold
  )
  scp2 = scale21(rk$output$residual.SD) * 100
  if(length(unique(scp2))>1){
    print(scp2)
    cptk = penCPD(
      scp2,
      threshold = threshold,
      plot  = plot,
      criteria = criteria,
      method = method,
      main = 'k - CPD solution path'
    )
  }else{
    cptk = 1
  }
  finalk = rk$output$k[cptk][which(rk$output$Obs.in.Interval[cptk] >=
                                     (min.obs + length(m)))][1]
  if (is.na(finalk)) {
    cat('\n An optimal k is not found. Max k will be used.')
    finalk = max(k)
  }
  
  
  
  ##### Plots
  if (plot & !is.na(finall)) {
    plot(
      rl$output$l,
      scale21(rl$output$residual.SD),
      ylim = c(0, 3),
      xaxt = 'n',
      ylab = 'Standardized sd',
      xlab = 'l',
      main = paste('Final l = ', round(finall, 5), sep = '')
    )
    
    abline(v = rl$output$l[cpt],
           col = 'grey',
           lty = 3)
    axis(
      side = 1,
      at = rl$output$l[cpt],
      labels = round(rl$output$l[cpt], 2),
      las = 3
    )
    abline(
      v = finall,
      col = 2,
      lwd = 2,
      lty = 3
    )
    text(
      x      = rl$output$l[cpt],
      y      = 2 + (cos(rl$output$Ind[cpt] / 2 / pi)),
      col    = 2 + ((-1) ^ rl$output$Ind[cpt]) ,
      labels = rl$output$Obs.in.Interval[cpt],
      cex    = .5
    )
    
  }
  if (plot & !is.na(finalk)) {
    plot(k,
         rk$output$BIC,
         xlab = 'k',
         ylab = 'BIC',
         main = 'k - BIC')
    abline(v = finalk)
    
    plot(
      rk$output$k,
      scale21(rk$output$residual.SD),
      ylim = c(0, 3),
      xaxt = 'n',
      ylab = 'Standardized sd',
      xlab = 'k',
      main = paste('Final k = ', round(finalk, 5), sep = '')
    )
    
    abline(v = rk$output$k[cptk],
           col = 'grey',
           lty = 3)
    axis(
      side = 1,
      at = rk$output$k[cptk],
      labels = round(rk$output$k[cptk], 2),
      las = 3
    )
    abline(
      v = finalk,
      col = 4,
      lwd = 2,
      lty = 3
    )
    text(
      x      = rk$output$k[cptk],
      y      = 2 + (cos(rk$output$Ind[cptk] / 2 / pi)),
      col    = 2 + ((-1) ^ rk$output$Ind[cptk]) ,
      labels = rk$output$Obs.in.Interval[cptk],
      cex    = .5
    )
  }
  
  
  ##### final model
  cat('\n 3|3 Forming the final model \n')
  finalr = gridSearch(
    t = t,
    x = x,
    y = y,
    m = t[m],
    l = finall,
    k = finalk ,
    mutInd = m,
    plot = plot,
    threshold = threshold
  )
  finalr$weights = finalr$weights[-1] #Remove index from the final weights
  return(
    list(
      final.k = finalk,
      final.l = finall,
      finalModel = finalr,
      model.l = rl,
      model.k = rk,
      x = x,
      y = y,
      l = l,
      k = k,
      m = m,
      min.obs   = min.obs,
      threshold = threshold,
      plot      = plot
    )
  )
  
}
