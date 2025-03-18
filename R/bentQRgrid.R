bentQRgrid <- function (y, x, z, tau) 
{
  checkfun <- function(u) {
    u * (tau - ifelse(u < 0, 1, 0))
  }
  qnregprof <- function(y, x, z, tau) {
    tt <- seq(quantile(x, 0.01), quantile(x, 0.99), length = 100)
    p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
    sst <- NULL
    bet <- matrix(0, length(tt), p + 3)
    for (kk in 1:length(tt)) {
      if (p == 0) {
        datnew <- data.frame(y = y, x1 = x, x2 = pmax(x - 
                                                        tt[kk], 0))
      }
      else {
        datnew <- data.frame(y = y, x1 = x, x2 = pmax(x - 
                                                        tt[kk], 0), z)
        names(datnew) = c("y", "x1", "x2", paste0("z", 
                                                  1:p, seq = ""))
      }
      fit <- rq(y ~ ., tau = tau, data = datnew, method = "br")
      bet[kk, ] <- fit$coef
      sst[kk] <- sum(checkfun(fit$residuals))
    }
    bp.grid <- min(tt[sst == min(sst)])
    bet.grid <- bet[tt == bp.grid, ]
    return(list(bet = bet.grid, bp = bp.grid))
  }
  fit.grid <- qnregprof(y, x, z, tau) 
  bp.est <- fit.grid$bp 
  return(list(bp.est = bp.est))
}
