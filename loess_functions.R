## this file contains loess functions that are re

loess.boot <- function(x, y, npoints=40, nreps=100, confidence=0.95, ...){
  dat <- stats::na.omit(data.frame(x=x,y=y))
  if(nrow(dat) == 0) stop ( "Error in dropping NA's")
  ndx <- order(dat$x)
  dat$x <- dat$x[ndx]
  dat$y <- dat$y[ndx]
  r <- range(dat$x, na.rm=TRUE)
  x.out <- seq(r[1], r[2], length.out=npoints)
  f <- stats::loess(y~x, data=dat, ...)
  y.fit <- stats::approx(f$x, stats::fitted(f), x.out,rule=2)$y
  len <- length(dat$x)
  mat <- matrix(0,nreps,length(x.out))
  for(i in seq(nreps)){
    ndx <- sample(len,replace=TRUE)
    x.repl <- x[ndx]
    y.repl <- y[ndx]
    f <- stats::loess(y.repl~x.repl, ...)
    mat[i,] <- stats::predict(f, newdata=x.out)
  }
  n.na <- apply(is.na(mat), 2, sum)  
  nx <- ncol(mat)
  up.lim <- rep(NA, nx)
  low.lim <- rep(NA, nx)
  stddev <- rep(NA, nx)
  for(i in 1:nx) {
    if(n.na[i] > nreps*(1.0-confidence)) {
      next
    }
    conf <- confidence*nreps/(nreps-n.na[i])
    pr <- 0.5*(1.0 - conf)
    up.lim[i] <- stats::quantile(mat[,i], 1.0-pr, na.rm=TRUE)
    low.lim[i] <- stats::quantile(mat[,i], pr, na.rm=TRUE)
    stddev[i] <- stats::sd(mat[,i], na.rm=TRUE)
  }
  ndx <- !is.na(up.lim)   
  fit <- data.frame(x=x.out[ndx], y.fit=y.fit[ndx], up.lim=up.lim[ndx],
                    low.lim=low.lim[ndx], stddev=stddev[ndx])
  fit.boot <- list(nreps=nreps, confidence=confidence,
                   span=f$pars$span, degree=f$pars$degree, 
                   normalize=f$pars$normalize, family=f$pars$family,
                   parametric=f$pars$parametric, surface=f$pars$surface,
                   data=dat, fit=fit)
  class( fit.boot ) <- "loess.boot"
  return( fit.boot )
}