# required packages
library(boot)
library(forecast)
library(ggplot2)
library(hash)

# generate random time series data (ARMA)
arma.sim <- function(ars, mas, nobs, c = 0, sigma = 1, seed = NA){

  # set seed
  if (!is.na(seed)){
    set.seed(seed)
  }

  x <- rnorm(nobs, mean = 0, sd = sigma)
  y <- vector()

  # initial data
  for (i in 1:max(length(ars), length(mas))) {
    y[i] <- x[i] + c
  }

  for (i in (max(length(ars), length(mas))+1):nobs) {
    y[i] <- x[i] + c
    for (a in 1:length(ars)) {
      y[i] <- y[i] + y[i-a]*ars[a] # calculate ar part
    }
    for (m in 1:length(mas)) {
      y[i] <- y[i] - x[i-m]*mas[m] # calculate ma part
    }
  }

  # transform to time series data
  y <- ts(y)
  return(y)
}


# plot the time series figure
arma.plot <- function(data){

  warnings("off")

  # transform to the data frame
  tdata <- data.frame(time = c(time(data)), data = data)

  # ggplot2
  tplot <- ggplot(tdata, aes(x = time, y = data)) + geom_line(colour = "red") +
    xlab("Sequence") +
    ylab("Data")

  return(tplot)
}

# block bootstrap for time series data
arma.boot <- function(data, len, resample = TRUE, seed = NA){

  # set seed
  if (!is.na(seed)){
    set.seed(seed)
  }

  # calculate windows
  ws <- length(data) %/% len
  ws_o <- length(data) %% len

  # sample the windows
  ws_sample <- sample(ws, ws, replace = resample)

  # joint
  ws_boot <- c()
  for (w in ws_sample) {
    ws_boot <- c(ws_boot, data[((w-1)*len+1):(w*len)])
  }

  ws_o_sample <- sample((length(data)-ws_o), 1)
  ws_boot <- c(ws_boot, data[ws_o_sample:(ws_o_sample+ws_o-1)])

  # transform to time series data
  ws_boot <- ts(ws_boot)

  return(ws_boot)
}

# boot base model
lynx.fun <- function(data){
  ar.fit <- ar(data, order.max = 25)
  c(ar.fit$order, mean(data), data)
}

# estimate the model uncertainty variance
arma.muv <- function(data, max.p = 5, max.q = 5, method = "aic", stepwise = TRUE,
                     blength = 15, btype = "mbb", bsamples = 100, seed = NA){

  # check the model selection method
  if (!(method %in% c("aic", "bic", "aicc"))){
    stop("The model selection method used is not in the candidate list!")
  }

  # check the bootstrap method
  if (!(btype %in% c("sb", "mbb"))){
    stop("The bootstrap method used is not in the candidate list!")
  }

  # set seed
  if (!is.na(seed)){
    set.seed(seed)
  }

  # generate bootstrap samples
  if (btype == "mbb"){
    boot.re <- tsboot(data, lynx.fun, R = bsamples, l = blength, sim = "fixed")
  }else{
    boot.re <- tsboot(data, lynx.fun, R = bsamples, l = blength, sim = "geom")
  }
  boot.data <- boot.re$t[, 3:(length(data)+2)]

  # fit bootstrap models
  a.order <- vector(mode = "numeric", length = bsamples)
  m.order <- vector(mode = "numeric", length = bsamples)
  for (b in 1:bsamples) {
    bdata <- boot.data[b, ]
    f.model <- auto.arima(bdata, ic = method, max.p = max.p, max.q = max.q, max.P = 0, max.Q = 0,
                          max.order = (max.p + max.q), max.D = 0, max.d = 0,
                          start.p = 1, start.q = 1, stepwise = stepwise)
    a.order[b] <- f.model$arma[1]
    m.order[b] <- f.model$arma[3]
  }

  muv <- var(a.order) + var(m.order)

  return(muv)
}

# estimate the model confidence bounds
arma.mcb <- function(data, max.p = 5, max.q = 5, method = "aic", stepwise = TRUE,
                     blength = 15, btype = "mbb", bsamples = 100, seed = NA){

  # check the model selection method
  if (!(method %in% c("aic", "bic", "aicc"))){
    stop("The model selection method used is not in the candidate list!")
  }

  # check the bootstrap method
  if (!(btype %in% c("sb", "mbb"))){
    stop("The bootstrap method used is not in the candidate list!")
  }

  # set seed
  if (!is.na(seed)){
    set.seed(seed)
  }

  # generate bootstrap samples
  if (btype == "mbb"){
    boot.re <- tsboot(data, lynx.fun, R = bsamples, l = blength, sim = "fixed")
  }else{
    boot.re <- tsboot(data, lynx.fun, R = bsamples, l = blength, sim = "geom")
  }
  boot.data <- boot.re$t[, 3:(length(data)+2)]

  # fit bootstrap models
  a.order <- vector(mode = "numeric", length = bsamples)
  m.order <- vector(mode = "numeric", length = bsamples)
  for (b in 1:bsamples) {
    bdata <- boot.data[b, ]
    f.model <- auto.arima(bdata, ic = method, max.p = max.p, max.q = max.q, max.P = 0, max.Q = 0,
                          max.order = (max.p + max.q), max.D = 0, max.d = 0,
                          start.p = 2, start.q = 2, stepwise = stepwise)
    if (f.model$arma[1] == 0){
      a.order[b] <- 1
    }else{
      a.order[b] <- f.model$arma[1]
    }
    if (f.model$arma[3] == 0){
      m.order[b] <- 1
    }else{
      m.order[b] <- f.model$arma[3]
    }
  }

  # bootstrap models
  bmodels <- list()
  for (b in 1:bsamples) {
    bmodels[[b]] <- c(a.order[b], m.order[b])
  }

  # full models
  fmodels <- list()
  f <- 1
  for (i in 1:max.p) {
    for (j in 1:max.q) {
      fmodels[[f]] <- c(i, j)
      f <- f + 1
    }
  }

  # create hash map for model sequence
  pmodels <- hash()
  for (m1 in fmodels) {
    for (m2 in fmodels) {
      mwidth <- (m2[1] - m1[1] + 1) * (m2[2] - m1[2] + 1)
      if (m2[1] >= m1[1] && m2[2] >= m1[2]){
        if (as.character(mwidth) %in% keys(pmodels)){
          pmodels[[as.character(mwidth)]] <- append(pmodels[[as.character(mwidth)]], list(c(m1, m2)))
        }else{
          pmodels[mwidth] <- list(c(m1, m2))
        }
      }
    }
  }

  # create hash map for model confidence bounds
  ptmcb <- hash()
  widths <- sort(as.numeric(keys(pmodels)))
  seqS <- 0

  for (w in widths) {
    wpmodel <- pmodels[[as.character(w)]]
    # index set
    gindex = list()
    # calculate r(wp) under the same width
    for (i in 1:length(wpmodel)) {
      wp <- wpmodel[[i]]
      bn <- 0
      plower <- wp[1:2]
      pupper <- wp[3:4]
      for (bm in bmodels) {
        if (plower[1] <= bm[1] && bm[1] <= pupper[1]){
          if (plower[2] <= bm[2] && bm[2] <= pupper[2]){
            bn <- bn + 1
          }
        }
      }
      gindex <- append(gindex, bn)
    }
    if (seqS == 0){
      ptmcb[w] <- wpmodel[[which.max(gindex)]]
      ptmcb[[as.character(w)]] <- append(ptmcb[[as.character(w)]],
                                       round(gindex[[which.max(gindex)]] / bsamples, 5))
      plastl <- ptmcb[[as.character(w)]][1]
      qlastl <- ptmcb[[as.character(w)]][2]
      plastu <- ptmcb[[as.character(w)]][3]
      qlastu <- ptmcb[[as.character(w)]][4]
      bcrlast <- round(gindex[[which.max(gindex)]] / bsamples, 5)
      seqS <- 1
    }else{
      while (length(gindex) >= 1) {
        sptmcb <- wpmodel[[which.max(gindex)]]
        pnowl <- sptmcb[1]
        qnowl <- sptmcb[2]
        pnowu <- sptmcb[3]
        qnowu <- sptmcb[4]
        if (pnowl <= plastl && pnowu >= plastu){
          if (qnowl <= qlastl && qnowu >= qlastu){
            bcrnow <- round(gindex[[which.max(gindex)]] / bsamples, 5)
            if (bcrnow >= bcrlast){
              ptmcb[w] <- wpmodel[[which.max(gindex)]]
              ptmcb[[as.character(w)]] <- append(ptmcb[[as.character(w)]],
                                                 round(gindex[[which.max(gindex)]] / bsamples, 5))
              bcrlast <- bcrnow
              plastl <- pnowl
              qlastl <- qnowl
              plastu <- pnowu
              qlastu <- qnowu
              break
            }
          }
        }
        ind <- which.max(gindex)
        gindex[[ind]] <- NULL
        wpmodel[[ind]] <- NULL
      }
    }
  }

  widths <- sort(as.numeric(keys(ptmcb)))
  lbm <- vector(mode = "character", length = length(widths))
  bcr <- vector(mode = "character", length = length(widths))
  ubm <- vector(mode = "numeric", length = length(widths))

  for (i in 1:length(widths)){
    mcb <- ptmcb[[as.character(widths[i])]]
    lbm[i] <- paste("(", mcb[1], ", ", mcb[2], ")", sep = "")
    ubm[i] <- paste("(", mcb[3], ", ", mcb[4], ")", sep = "")
    bcr[i] <- mcb[5]
  }

  tmcb <- data.frame(width = widths, lbm = lbm, bcr = bcr, ubm = ubm)

  return(tmcb)
}

# plot the model uncertainty curve
arma.muc <- function(method, ...){

  mcbs <- list(...)

  cols <- c("#8ECFC9", "#FFBE7A", "#FA7F6F")

  if (length(method) == 1){
    mcb <- mcbs[[1]]
    mcb$x <- mcb$width / max(mcb$width)
    plot(mcb$bcr ~ mcb$x, type = "l", xlim = c(0,1), ylim = c(0,1),
         xlab = "msc/f", ylab = "Bootstrap Coverage Rate", lwd = 3, col = "#8ECFC9")
    legend("bottomright", legend = c(method), col = c("#8ECFC9"), lty = 1, lwd = 2)
  }else{
    mcb <- mcbs[[1]]
    mcb$x <- mcb$width / max(mcb$width)
    plot(mcb$bcr ~ mcb$x, type = "l", xlim = c(0,1), ylim = c(0,1),
         xlab = "msc/f", ylab = "Bootstrap Coverage Rate", lwd = 3, col = "#8ECFC9")
    for (m in 2:length(method)) {
      mcb <- mcbs[[m]]
      mcb$x <- mcb$width / max(mcb$width)
      lines(mcb$bcr ~ mcb$x, type = "l", xlim = c(0,1), ylim = c(0,1),
           xlab = "msc/f", ylab = "Bootstrap Coverage Rate", lwd = 3, col = cols[m])
    }

    legend("bottomright", legend = c(method), col = c(cols[1:length(method)]), lty = 1, lwd = 2)
  }

}
