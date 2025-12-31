#' Compute Rate Ratio
#'
#' @param object An object.
#' @param ... Passed to methods.
#' @export
RateRatio <- function(object, ...) {
  UseMethod("RateRatio")
}


#' Title Rate Ratio for aceDLNMadditive_fit objects
#'
#' @param object object of class \code{aceDLNMadditive_fit}.
#' @param verbose whether to print messages during the process, default \code{FALSE}.
#'
#' @return
#' @importFrom mgcv s
#' @export
RateRatio.aceDLNMadditive_fit <- function(object, x0, x1, verbose = FALSE, ...) {
  ## point estimate
  pc <- object$data$pc
  alpha_w.list <- object$point$alpha_w.list
  alpha_f.list <- object$point$alpha_f.list
  M <- length(alpha_f.list)
  kw <- object$data$kw
  kE <- object$data$kE

  Ufpen.list <- object$data$Ufpen.list
  Uwpen <- object$data$Uwpen
  Zwnew <- object$data$Zw %*% Uwpen
  maxL <- object$data$maxL

  sX.x.names.list <- object$data$sX.x.names.list
  R.CI <- nrow(object$CI.sample[[1]])



  out.list <- lapply(1:M, function(i) {
    ## point estimate
    wl.fit <- function(lnew) c(1, Bsplinevec2Con(lnew, object$data$knots_w, 4, Zwnew)) %*% alpha_w.list[[i]]
    wl.fit <- Vectorize(wl.fit)
    wl.discrete <- sapply(0:maxL, function(ll) {
      integrate(wl.fit, lower = ll, upper = ll + 1)$value
    })
    E0 <- sum(wl.discrete * x0[,i])
    E1 <- sum(wl.discrete * x1[,i])

    fE.fit <- function(Enew) mgcv::PredictMat(object$smooth$fE[[i]], data = data.frame(E = Enew)) %*% Ufpen.list[[i]] %*% alpha_f.list[[i]]

    eta.E.0 <- fE.fit(E0)
    eta.E.1 <- fE.fit(E1)

    RR.est <- data.frame(eta.E.0 = eta.E.0,
                         eta.E.1 = eta.E.1,
                         RR = exp(eta.E.1 - eta.E.0))



    ## CI
    w.id <- ((i-1)*kw+1) : (i*kw)
    f.id <- ((i-1)*(kE-1)+1) : (i*(kE-1))
    RR.sample <- lapply(1:R.CI, function(j) {
      if(verbose) {
        if(j %% 100 == 0) {
          cat("Sampling iteration:", j, " / ", R.CI, " for exposure",  i, "\n")
        }
      }
      alpha_w_sample <- object$CI.sample$alpha_w_sample[j, w.id]
      wl.fit <- function(lnew) c(1, Bsplinevec2Con(lnew, object$data$knots_w, 4, Zwnew)) %*% alpha_w_sample
      wl.fit <- Vectorize(wl.fit)


      wl.discrete <- sapply(0:maxL, function(ll) {
        integrate(wl.fit, lower = ll, upper = ll + 1)$value
      })
      E0 <- sum(wl.discrete * x0[,i])
      E1 <- sum(wl.discrete * x1[,i])

      fE.fit <- function(Enew) mgcv::PredictMat(object$smooth$fE[[i]], data = data.frame(E = Enew)) %*% Ufpen.list[[i]] %*% alpha_f.list[[i]]

      eta.E.0 <- fE.fit(E0)
      eta.E.1 <- fE.fit(E1)

      return(data.frame(eta.E.0 = eta.E.0,
                        eta.E.1 = eta.E.1,
                        RR = exp(eta.E.1 - eta.E.0)))
    })


    RR.sample <- data.table::rbindlist(RR.sample)
    return(list(RR.est = RR.est,
                RR.sample = RR.sample))
  })

  out.list
}
