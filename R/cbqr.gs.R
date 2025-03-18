#' Grid search of change points for bent line quantile regression
#' @description
#' It finds the change points using grid search. The candidates are equally-spaced grid points from 1% quantile to 99% quantile of the covariate that may contain a change point. The change point is found by minimizing the quantile check loss.
#' @param y the response
#' @param x the covariate that may contain a change point
#' @param z the covariates that do not contain a change point
#' @param tau quantile level
#' @return bp.est estimates of change points

t.search <- function (y, x, z, tau)
{
  checkloss <- function(u) {
    u * (tau - ifelse(u < 0, 1, 0))
  }
  grid.search <- function(y, x, z, tau) {
    tt <- seq(quantile(x, 0.01), quantile(x, 0.99), length = 100)
    p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
    sst <- NULL
    bet <- matrix(0, length(tt), p + 3)
    for (kk in 1:length(tt)) {
      if (p == 0) {
        datnew <- data.frame(y = y, x1 = x, x2 = pmax(x - tt[kk], 0))
      }
      else {
        datnew <- data.frame(y = y, x1 = x, x2 = pmax(x - tt[kk], 0), z)
      }
      fit <- rq(y ~ ., tau = tau, data = datnew, method = "br")
      bet[kk, ] <- fit$coef
      sst[kk] <- sum(checkloss(fit$residuals))
    }
    t.hat <- min(tt[sst == min(sst)])
    return(list(t.hat = t.hat))
  }
  grid <- grid.search(y, x, z, tau)
  t.hat <- t.hat$bp
  return(list(t.hat = t.hat))
}


#' Censored bent line quantile regression
#'
#'@description
#' It implements the methodology censored bent line quantile regression Y^* = beta_1 t + beta_2 max(t - t_0, 0) + X gamma + epsilon, subject to censoring. In the first step, it estimates an informative subset, which is estimated by generalized additive model (GAM) with a probit link. In the second step, it estimates a bent line quantile regression model. The estimation of the change point is based on grid search, and the estimation of regression coefficients is based on rq in the package quantreg. The standard errors and the interval estimates are provided with three options ('nid', 'ker', 'ks'). See details in summary.rq for 'nid' (quotient) and 'ker' (kernel density estimation). 'ks' considers a bias correction using the derived Bahadur representation for the residuals in kernel density estimation. For smaller sample sizes, 'ker' is more recommended as 'ks' may produce wide intervals occasionally.
#' @param cy the response (including the censored response)
#' @param tt the covariate that may contain a change point
#' @param xx the covariates that do not contain a change point
#' @param tau quantile level
#' @param delta the censoring indicator (0 if censored)
#' @param level significance level. The default is 5%.
#' @param censoring the censoring type. 'right' means censoring from below while ’left’ means censoring from above. The default is ’right’.
#' @param method the methods for density estimation involved in asymptotic interval estimation. The options include 'nid', 'ker', and 'ks'. The default is 'ks'.
#' @return beta.est estimates of beta_1, beta_2 along with standard error estimates
#' @return t.est estimates of $t$ along with standard error estimates
#' @return ci.t confidence interval of change point
#' @return ci.beta confidence intervals of beta_1, beta_2
#' @return gamma.est estimates of gamma along with standard error estimates
#' @references Li, Chenxi, Ying Wei, Rick Chappell, and Xuming He. "Bent line quantile regression with application to an allometric study of land mammals' speed and mass." Biometrics 67, no. 1 (2011): 242-249.
#' @references Tang, Yanlin, Huixia Judy Wang, Xuming He, and Zhongyi Zhu. "An informative subset-based estimator for censored quantile regression." Test 21 (2012): 635-655.
#' @examples
#' m <- 10
#' n <- 50
#' tau <- 0.50
#' cl <- 4
#' t0 <- 1.5
#' subj.m <- seq(from = .5, to = 3, length.out = m)
#' yy <- c()
#' tt <- c()
#' cy <- c()
#' for (j in 1:m) {
#'  t <- runif(n, 0, 3)
#'  t <- t[order(t)]
#'  tt <- c(tt, t)
#'  y <- subj.m[j] + .5 * t + 2 * pmax(t - t0, 0) + rnorm(n, sd = .1)
#'
#'  # censoring
#'  c <- y
#'  if (min(which(y > cl)) < Inf) {
#'    c[min(which(y > cl)):length(c)] <- cl
#'  }
#'
#'  yy <- c(yy, y)
#'  cy <- c(cy, c)
#'}
#' # covariate
#' group <- as.factor(sapply(1:m, function(x) rep(x, n)))
#' xx <- as.matrix(model.matrix(~ group))
#' xx <- xx[, - 1]
#' ## GAM estimation
#' delta <- ifelse(yy <= cl, yes = 1, no = 0)
#' mod.gs <- cbqr.gs(cy = cy, tt = tt, xx = xx, tau = tau, delta = delta)
#' @export

# censored bent line quantile regression (grid search)
cbqr.gs <- function(cy, tt, xx, tau, delta, level = .05,
                    censoring = 'right',
                    method = 'ks') {

  ## GAM estimation
  n <- length(cy)
  mod.glm <- gam(delta ~ tt + xx, family = binomial(link = 'probit'))
  pi.hat <- predict(mod.glm, type = 'response')
  if (censoring == 'right') {
    ISUB.est <- which(pi.hat > tau + n ^ (- 1 / 4) * tau)
  } else {
    ISUB.est <- which(pi.hat > 1 - tau + n ^ (- 1 / 4) * tau)
  }
  n.isub <- length(ISUB.est)
  tt.isub <- tt[ISUB.est]
  xx.isub <- xx[ISUB.est, ]

  ## Bent QR with estimated ISUB
  fit.grid <- t.search(y = cy[ISUB.est], z = xx.isub, x = tt.isub, tau)
  t.hat <- fit.grid$t.hat
  mod.isub <- rq(cy ~ tt + pmax(tt - t.hat, 0) + xx,
                 tau = tau, subset = ISUB.est,
                 method = "br")

  # slope estimates
  beta1.hat <- coef(mod.isub)[2]
  beta2.hat <- coef(mod.isub)[3]

  # design matrix
  h <- cbind(ifelse(tt.isub <= t.hat, yes = 1, no = 0),
             tt.isub * ifelse(tt.isub <= t.hat, yes = 1, no = 0),
             ifelse(tt.isub > t.hat, yes = 1, no = 0),
             tt.isub * ifelse(tt.isub > t.hat, yes = 1, no = 0),
             xx.isub)
  g1 <- cbind(1, tt.isub - t.hat, 0, - beta1.hat, xx.isub)
  g2 <- cbind(1, 0, tt.isub - t.hat, - beta2.hat, xx.isub)
  g.tilde <- g1 * ifelse(tt.isub <= t.hat, yes = 1, no = 0) +
    g2 * ifelse(tt.isub > t.hat, yes = 1, no = 0)

  # density estimates
  resid <- residuals(mod.isub)


  if (method == 'ks') {


    # homoskedastic density estimate
    f0.hat <- akj(resid, 0)$dens
    D.hat <- - 1 / n * f0.hat * t(h) %*% g.tilde
    D.hat.inv <- solve(D.hat)
    # score
    psi.r <- function(r, tau) tau - ifelse(r < 0, 1, 0)
    br <- - (1 / n) * h %*% D.hat.inv %*% t(t(psi.r(r = resid, tau = tau)) %*% h)

    # .05 .95
    bounds <- quantile(resid + br, c(.05, .95))
    rr <- bounds[2] - bounds[1]
    # rr <- mad(as.numeric(resid + br))
    # 2
    bandwidth <- .2 * rr * n ^ (-1 / 3)
    kd <- dnorm(as.numeric(resid + br) / bandwidth) / bandwidth
    f.mat <- diag(kd)


  } else if (method == 'nid') {

    h.n <- bandwidth.rq(tau, n)
    while ((tau - h.n < 0) || (tau + h.n > 1)) h.n <- h.n/2

    # model fit for higher quantile level
    fit.grid.hi <- bentQRgrid(y = cy[ISUB.est], z = xx.isub, x = tt.isub, tau + h.n)
    bp.ini.hi <- fit.grid.hi$bp.est
    mod.hi <- rq(cy ~ tt + pmax(tt - bp.ini.hi, 0) + xx,
                 tau = tau + h.n, subset = ISUB.est,
                 method = "br")

    # model fit for lower quantile level
    fit.grid.lo <- bentQRgrid(y = cy[ISUB.est], z = xx.isub, x = tt.isub, tau - h.n)
    bp.ini.lo <- fit.grid.lo$bp.est
    mod.lo <- rq(cy ~ tt + pmax(tt - bp.ini.lo, 0) + xx,
                 tau = tau - h.n, subset = ISUB.est,
                 method = "br")

    dyhat <- fitted.values(mod.hi) - fitted.values(mod.lo)
    f.mat <- diag(pmax(0, (2 * h.n)/(dyhat - 1e-4)))



  } else if (method == 'ker') {

    h.n <- bandwidth.rq(tau, n.isub)
    while ((tau - h.n < 0) || (tau + h.n > 1)) h.n <- h.n/2
    h.n <- (qnorm(tau + h.n) - qnorm(tau - h.n)) * min(sqrt(var(resid)),
                                                       (quantile(resid, 0.75) - quantile(resid, 0.25))/1.34)
    f <- dnorm(resid / h.n) / h.n
    f.mat <- diag(f)

  }


  # covariance estimation
  D.hat <- - 1 / n * t(h) %*% f.mat %*% g.tilde
  D.hat.inv <- solve(D.hat)
  C.hat <- 1 / n * t(h) %*% h
  Sigma.hat <- tau * (1 - tau) * D.hat.inv %*% C.hat %*% t(D.hat.inv)

  ci.t <- t.hat + c(-1, 1) * qnorm(1 - level) * sqrt(Sigma.hat[4, 4] / n)
  ci.beta1 <- beta1.hat + c(-1, 1) * qnorm(1 - level) * sqrt(Sigma.hat[2, 2] / n)
  ci.beta2 <- beta1.hat + beta2.hat + c(-1, 1) * qnorm(1 - level) * sqrt(Sigma.hat[3, 3] / n)
  ci.beta <- rbind(ci.beta1, ci.beta2 - beta1.hat)

  beta.est <- cbind(coef(mod.isub)[2:3], sqrt(c(Sigma.hat[2, 2], Sigma.hat[3, 3]) / n))
  rownames(beta.est) <- c('beta1', 'beta2')
  colnames(beta.est) <- c('Estimate', 'SE')

  t.est <- c(t.hat, sqrt(Sigma.hat[4, 4] / n))
  names(t.est) <- c('Estimate of change point', 'SE')
  names(ci.t) <- c('LB', 'UB')

  gamma.est <- cbind(coef(mod.isub)[-(2:3)], sqrt(diag(Sigma.hat)[-(1:4)] / n))
  colnames(gamma.est) <- c('Estimate', 'SE')

  list(beta.est = beta.est,
       t.est = t.hat,
       ci.t = ci.t,
       ci.beta = ci.beta,
       gamma.est = gamma.est)
}


