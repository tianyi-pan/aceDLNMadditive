#' Title Generate data
#'
#' @param fEtype
#' @param wltype
#' @param Nt
#' @param kx.per500
#' @param theta
#' @param maxL
#' @param kw
#'
#' @return
#' @export
#' @importFrom mgcv s
#'
#' @examples
GenerateData <- function(fEtype.vec, wltype.vec, Nt = 1000, kx.per500 = 100,
                         interpolate = TRUE,
                         theta = 10,
                         maxL = 14,
                         exposure.list = NULL,
                         std = FALSE,
                         verbose = FALSE,
                         c1 = NULL,
                         c2 = NULL,
                         wl.set = NULL,
                         fE.set = NULL,
                         other) {

  ## load waterloo PM2.5 data (PM25.waterloo)
  if(is.null(exposure.list)) data("AirPollutionWaterloo")
  maxLreal <- maxL+1

  t <- 1:(Nt+maxL) # time 1 to (Nt+40)
  t.sim <- t[-(1:maxL)]



  if((kx.per500 > 300) || (interpolate == TRUE)) {
    # cat("Interpolate the exposure process.")
    kx <- Nt+maxL + 4 + 2 # number of knots for X(t)
    interpolate <- TRUE
  } else{
    kx <- kx.per500 * ifelse(Nt < 500, 1, round(Nt/500)) # number of knots for X(t)
  }

  ## exposure process

  ## exposure process
  if(is.null(exposure.list)) {
    PM25 <- AirPollution.Waterloo$PM25[t]
    O3 <- AirPollution.Waterloo$O3[t]
    NO2 <- AirPollution.Waterloo$NO2[t]
    exposure.list <- list(PM25, O3, NO2)
  }
  if(std) {
    PM25 <- (PM25 - mean(PM25)) / sd(PM25)
    O3 <- (O3 - mean(O3)) / sd(O3)
    NO2 <- (NO2 - mean(NO2)) / sd(NO2)
    exposure.list <- list(PM25, O3, NO2)
  }


  wl.list <- lapply(wltype.vec, function(wltype) {
    switch (wltype,
            type1 = {
              wlc <- function(l) dnorm(l, mean = 3, sd = 3.5)^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) dnorm(l, mean = 3, sd = 3.5)/wl_de
            },
            type2 = {
              wlc <- function(l) (1/(1+exp(l-8)))^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) 1/(1+exp(l-8))/wl_de
            },
            type3 = {
              wlc <- function(l) (1/(1+exp(l-1)))^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) 1/(1+exp(l-1))/wl_de
            },
            type4 = {
              wlc <- function(l) dnorm(l,8,10)^2
              wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
              wl <- function(l) dnorm(l,8,10)/wl_de
            }
            # symmetric = {
            #   wlc <- function(l) ( dnorm(l, mean = 7.5, sd = 3.2)-dnorm(0, mean = 7.5, sd = 3.2) )^2
            #   wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
            #   wl <- function(l) ( dnorm(l, mean = 7.5, sd = 3.2)-dnorm(0, mean = 7.5, sd = 3.2) )/wl_de
            # },
            # constant = {
            #   wlc <- function(l) as.numeric((c1 <= l) & (l <= c2))^2
            #   wl_de <- sqrt(integrate(wlc, lower = 0, upper = maxLreal)$value)
            #   wl <- function(l) as.numeric((c1 <= l) & (l <= c2))/wl_de
            # }

    )
    return(list(wlc = wlc,
                wl_de = wl_de,
                wl = wl))
  })




  # model exposure process
  if(!interpolate) {
    SSx <- mgcv::smoothCon(s(t, bs = "bs", k = kx),
                           absorb.cons = FALSE,
                           data = data.frame(t = t))[[1]] ## reparameterize it later
    knots_x <- SSx$knots
    X <- SSx$X
    ## sum-to-zero reparameterization for SSx
    QRx <- qr(t(X) %*% as.vector(rep(1,nrow(X))))
    Qx <- qr.Q(QRx, complete = TRUE)
    Zx <- Qx[,2:ncol(Qx)]
    ## Check whether the design matrices are identical
    # X_repa <- SSx$X %*% Zx
    # max(unname(model.matrix(xt.fit))[,-1] - SSx$X %*% Zx) # SAME
  } else {
    # knots_x <- c(rep(0,2), t, rep(Nt+maxL+1,2))
    # X <- splines::splineDesign(knots = knots_x, x = c(0, t, Nt+maxL+1), outer.ok = TRUE)
    # X <- X[-c(1,Nt+2),]
    # knots_x <- c(rep(t[1]-1-0.2,3),t[1]-1-0.2, t[1]-1-0.001, t, Nt+maxL+0.5+0.001, Nt+maxL+0.5+0.2, rep(Nt+maxL+0.5+0.2,3))
    # knots_x <- c(rep(t[1]-1-0.2,4), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,4))
    knots_x <- c(t[1]-1-0.3, rep(t[1]-1-0.2,3), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,3), t[length(t)]+1+0.3)
    X <- splines::splineDesign(knots = knots_x, x = knots_x, outer.ok = TRUE)
    # X <- X[-c(1:2,Nt+maxL+7,Nt+maxL+8),]
  }


  E.list <- lapply(1:length(exposure.list), function(kk) {
    x <- exposure.list[[kk]]

    wlc = wl.list[[kk]]$wlc
    wl_de = wl.list[[kk]]$wl_de
    wl = wl.list[[kk]]$wl

    if(!interpolate){
      xt.fit <- mgcv::gam(x~s(t, bs = "bs", k = kx), data = data.frame(x = PM25, t = t))
      ## coefficients for X(t)
      alpha_x <- xt.fit$coefficients

      # plot(t, xt.fit$fitted.values, type = "l", ylim = c(0,20))
      # lines(t, PM25, col = "red")
      Xt <- function(tnew) as.numeric(c(1, Bsplinevec2Con(tnew, knots_x, 4, Zx)) %*% alpha_x)
    } else {
      ## interpolate = TRUE
      ## set points for boundary and auxiliary boundary
      # alpha_x <- Interpolate(X, Zx, c(rep(0,4),PM25,rep(0,4)))
      Xsparse <- as(X, "dgCMatrix")
      # alpha_x <- Interpolate(Xsparse, c(rep(0,5),PM25,rep(0,5)))
      alpha_x <- Interpolate(Xsparse, c(rep(0,4), x[1], x, x[length(x)], rep(0,4)))
      xt.fit <- "interpolate"
      Xt <- function(tnew) as.numeric(Bsplinevec2(tnew, knots_x, 4) %*% alpha_x)
    }



    ## The true exposure process is termed as the fitted function

    # CHECK: sapply(t, Xt) - xt.fit$fitted.values
    # Check: interpolate
    # Xtfit <- sapply(t, Xt)
    # max(abs(Xtfit - PM25)) #<1e-12



    ## True weighted exposure E = \int X(t-l) w(l) d l
    # E <- sapply(t, function(t.) integrate(Vectorize(function(l) wl(l)*Xt(t.-l)), lower = 0, upper = maxL)$value)


    ## True weighted exposure E = \int X(t-l) w(l) d l
    # E <- sapply(t, function(t.) integrate(Vectorize(function(l) wl(l)*Xt(t.-l)), lower = 0, upper = maxL)$value)
    l.eval <- seq(0, maxLreal, length.out = max(maxLreal, 500))
    wl.true <- sapply(l.eval, wl)

    wl.true.fit <- mgcv::gam(wl~s(l, bs = "bs", k = max(maxLreal, 100)), data = data.frame(wl = wl.true, l = l.eval))
    SSw.true <- mgcv::smoothCon(s(l, bs = "bs", k = max(maxLreal, 100)), absorb.cons = FALSE, data = data.frame(l = l.eval))[[1]]
    QRw.true <- qr(t(SSw.true$X) %*% as.vector(rep(1,nrow(SSw.true$X))))
    Qw.true <- qr.Q(QRw.true, complete = TRUE)
    Zw.true <- Qw.true[,2:ncol(Qw.true)]


    if(!interpolate) {
      if(verbose) start <- Sys.time()
      integ.E <- Integral(knots_x, SSw.true$knots, kx, max(maxLreal, 100), maxLreal, Zx, Zw.true, t.sim+0.5, alpha_x, TRUE)
      if(verbose) cat("integral 1: ", Sys.time() - start, "\n")
    }else {
      integ.E <- Integral_interpolate(knots_x, SSw.true$knots, kx, max(maxLreal, 100), maxLreal, Zw.true, t.sim+0.5, alpha_x, TRUE)
    }



    E.sim <- c(integ.E$AlphaxD %*% wl.true.fit$coefficients)



    ## association of weighted exposure f(E)
    Emax <- ceiling(max(c(E.sim)))
    Emin <- floor(min(c(E.sim)))
    return(list(E.sim = E.sim,
                Emax = Emax,
                Emin = Emin))
  })


  fE.list <- lapply(1:length(exposure.list), function(kk) {
    fEtype <- fEtype.vec[[kk]]
    Emax <- E.list[[kk]]$Emax
    Emin <- E.list[[kk]]$Emin
    switch (fEtype,
            cubic = {
              fEtmp <- function(E) (E-(Emin + (Emax-Emin)*0.3))*(E-(Emin + (Emax-Emin)*0.2))*(E-(Emin + (Emax-Emin)*0.9))
              fE <- function(E) (fEtmp(E) / (fEtmp(Emax)/2.5) + 2.5) /1.7
            },
            linear = {
              fEtmp <- function(E) 0
              fE <- function(E) ((E-(Emax+Emin)/2) / ((Emax - Emin)/3.5) + 3) /1.7
            },
            # constant = {
            #   fEtmp <- function(E) 0
            #   fE <- function(E) 0
            # },
            quadratic = {
              fEtmp <- function(x){25*(dnorm(((2*((x-Emin)/(Emax-Emin)+0.18) - 1.1))))}
              fE <- function(x) (-fEtmp(x)+11) /1.7
            }
            )
    return(list(fEtmp = fEtmp,
                fE = fE,
                Emax = Emax,
                Emin = Emin))
  })

  ## linear predictor
  eta.sim <- apply(sapply(1:length(exposure.list), function(kk) {
    sapply(E.list[[kk]]$E.sim, fE.list[[kk]]$fE)
  }), 1, sum)
  # eta.sim <- sapply(E.sim1, fE) + sapply(E.sim2, fE) + sapply(E.sim3, fE)
  # eta.sim <- sapply(E.sim1, fE) + sapply(E.sim3, fE)
  Esurface.sim <- eta.sim - mean(eta.sim)
  if (!missingArg(other)) {
    eta.sim <- eta.sim + apply(other, 1, sum)
  }

  ## generate data
  y.sim <- sapply(eta.sim, function(eta.) rnbinom(n = 1, size = theta, mu = exp(eta.)))


  return(list(y.sim = y.sim,
              t.sim = t.sim,
              E.list = E.list,
              eta.sim = eta.sim,
              Esurface.sim = Esurface.sim,
              y = c(rep(0, maxL), y.sim),
              x1 = PM25,
              x2 = O3,
              x3 = NO2,
              t = t,
              # xt.fit = xt.fit,
              # alpha_x = list(alpha_x1, alpha_x2, alpha_x3), # for interpolation
              # knots_x = knots_x, # for interpolation
              true.f = list(x1 = list(fE = fE.list[[1]]$fE, Emax = fE.list[[1]]$Emax, Emin = fE.list[[1]]$Emin, fEtmp = fE.list[[1]]$fEtmp,
                                      wl = wl.list[[1]]$wl, wl_de = wl.list[[1]]$wl_de),
                            x2 = list(fE = fE.list[[2]]$fE, Emax = fE.list[[2]]$Emax, Emin = fE.list[[2]]$Emin, fEtmp = fE.list[[2]]$fEtmp,
                                      wl = wl.list[[2]]$wl, wl_de = wl.list[[2]]$wl_de),
                            x3 = list(fE = fE.list[[3]]$fE, Emax = fE.list[[3]]$Emax, Emin = fE.list[[3]]$Emin, fEtmp = fE.list[[3]]$fEtmp,
                                      wl = wl.list[[3]]$wl, wl_de = wl.list[[3]]$wl_de))
  ))
}
