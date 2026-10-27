#' Prepare an additive ACE-DLNM for Super Learner
#'
#' @param object An object of class \code{aceDLNMadditive_fit}.
#' @param verbose Whether to print progress, default \code{FALSE}.
#' @param AIC.only Used by \code{mDLNM::AICdata()} to prepare AIC quantities
#'   without leave-two-out predictions.
#' @param nthreads Number of threads for deletion calculations, default \code{1}.
#'   Uses OpenMP for parallel calculation.
#' @param dat.c Counterfactual data, including the full exposure histories, on
#'   the same predictor scale expected by \code{predict()}. It must retain the
#'   fitted reference records and model structure. This required input defines
#'   the intervention.
#' @param ... Unused; additional arguments are rejected.
#' @exportS3Method mDLNM::SLdata
SLdata.aceDLNMadditive_fit <- function(object, dat.c, verbose = FALSE,
                                     AIC.only = FALSE, nthreads = 1, ...) {
  if(identical(object$family, "gaussian")) stop("SLdata() requires a negative binomial fit.")
  if(missing(dat.c) || is.null(dat.c)) stop("dat.c is required.")
  if(length(list(...))) stop("Unused arguments.")

  dat <- object$data
  point <- object$point
  coef.samples <- NULL
  if(AIC.only) {
    samples <- object$CI.sample
    if(is.null(samples$alpha_f_sample) || nrow(samples$alpha_f_sample) < 2) {
      stop("AICdata() requires at least two stored coefficient samples.")
    }
    coef.samples <- cbind(samples$alpha_f_sample, samples$phi_sample,
                          samples$betaR_sample, samples$betaF_sample)
  }

  mod.build <- aceDLNMadditivebuild(dat$y, dat$B_inner, dat$knots_f.list,
                                  dat$SwI_large, dat$SfI.list, dat$Dw,
                                  dat$Xrand, dat$Xfix,
                                  lapply(dat$SSf.list, "[[", "Zfnew"),
                                  dat$offset$Xoffset, dat$r, dat$K, dat$a,
                                  point$alpha_f, point$phi, point$log_theta,
                                  point$log_smoothing_f, point$log_smoothing_w,
                                  point$betaR, point$betaF, point$log_smoothing)
  mod.c <- stats::predict(object, dat = dat.c, .build.only = TRUE)
  result <- SLaceDLNMadditive(mod.build$address.eigen, mod.c$address.eigen,
                             verbose, AIC.only, nthreads, coef.samples)
  theta <- exp(point$log_theta)

  if(AIC.only) {
    aic.names <- c("AIC_conventional", "AIC_cAIC")
    colnames(result$aic.loo) <- aic.names
    return(list(y = dat$y, row.id = rownames(object$modeldata), theta = theta,
                aic = stats::setNames(as.vector(result$aic), aic.names),
                aic.loo = result$aic.loo,
                delta = result$delta, delta.loo = as.vector(result$delta.loo),
                variance = stats::var(as.vector(result$delta.sample)), estimand = "ratio"))
  }

  return(list(y = dat$y,
              row.id = rownames(object$modeldata),
              theta = theta,
              mu.loo = as.vector(result$mu.loo),
              mu.lto = result$mu.lto,
              delta = result$delta,
              delta.loo = as.vector(result$delta.loo),
              estimand = "ratio"))
}
