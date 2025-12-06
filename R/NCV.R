#' Compute NCV
#'
#' @param object An object.
#' @param ... Passed to methods.
#' @export
NCV <- function(object, ...) {
  UseMethod("NCV")
}

#' Title Obtain NCV from the fitted model by aceDLNMadditive
#'
#' @param object object of class \code{aceDLNMadditive_fit}.
#' @param kNCV number of neighbors on each side, default \code{0} which is LOOCV.
#' @param NCV.nthreads number of threads for NCV calculation, default \code{1}.
#' @param verbose whether to print messages during the process, default \code{FALSE}.
#' @param ...
#'
#' @export
NCV.aceDLNMadditive_fit <- function(object, kNCV = 0, NCV.nthreads = 1, verbose = FALSE, ...) {

  dat <- object$inputdata

  

  smooth = object$formula$smooth
  formula = object$formula$formula
  fe.cont = object$formula$fe.cont
  fe.varying = object$formula$fe.varying
  
  sXformula.list <- as.list(trimws(strsplit(deparse(formula[[3]]), "\\+")[[1]]))
  
  ## change character columns to factor. support random effects defined by mgcv::s(bs = "re")
  chr_col <- which(sapply(dat, class) == "character")
  if(length(chr_col) >= 1) {
    for(col. in chr_col) dat[,col.] <- factor(dat[,col.])
  }

  kw = object$data$kw
  kE = object$data$kE
  maxL = object$data$maxL
  maxLreal = maxL+1



  sXobject1 <- eval(parse(text = sXformula.list[[1]]))

  sXdat <- dat
  colnames(sXdat)[which(colnames(sXdat) == as.character(formula[[2]]))] <- "y"
  colnames(sXdat)[which(colnames(sXdat) == sXobject1$t)] <- "t"


  byvar <- sXobject1$by
  maxLreal <- maxL+1


  sX.x.names.list <- lapply(sXformula.list, function (sXformula) {
    return(unlist(eval(parse(text = sXformula))$x))
  })


  ### CONSTRUCTIONS #######
  ## time non-varying for group
  if(length(unique(sXdat$t)) < nrow(sXdat)) {

    if(is.null(fe.cont) & is.null(byvar)) stop("The exposure process data are duplicated at some time point. Please provide fe.cout.")

    if(!is.null(fe.cont)){
      # group_name.fe.cont <- fe.cont[[2]]
      # if(length(group_name.fe.cont) >= 2) {
      #   group_name.fe.cont <- as.character(group_name.fe.cont[2:length(group_name.fe.cont)])
      # } else {
      #   group_name.fe.cont <- as.character(group_name.fe.cont)
      # }
      group_name.fe.cont <- all.vars(fe.cont)
    } else {
      group_name.fe.cont <- NULL
    }
    if(!is.null(byvar)) {
      group_name.byvar <- byvar
    } else {
      group_name.byvar <- NULL
    }

    group_name <- c(group_name.fe.cont, group_name.byvar)

    ## split data
    sXdatlist <- split(as.data.table(sXdat), by = group_name) # use split in data.table to preserve order

  } else {
    group_name <- NULL
    sXdatlist <- list(sXdat)
  }

  ## set t starting at 1 and sort
  min_t <- min(sXdat$t)
  sXdatlist <- lapply(sXdatlist, function(sXdati) {
    sXdati$t <- sXdati$t - min_t + 1
    setorder(sXdati, t)
    return(data.frame(sXdati))
  })

  ## smooths from mgcv
  SwI_large <- object$data$SwI_large
  K <- object$data$K
  a <- object$data$a

  knots_w <- object$data$knots_w
  kx.per500 <- object$data$kx.per500
  interpolate <- object$interpolate

  Zw <- object$data$Zw
  Uwpen <- object$data$Uwpen
  Zwnew <- Zw %*% Uwpen


  ### 0. distributed lag term
  ## preparations for distributed lag terms


  sX.DL <- lapply(sX.x.names.list, function(xx) {

    DLprepare <- lapply(sXdatlist, function(sXdati){

      x <- sXdati[[xx]]
      t <- sXdati$t
      removed.t <- t[1:maxL] # removed id for the original t
      t <- t - min(t) + 1 # set t starting at 1

      Nti <- length(x) - maxL
      if((kx.per500 > 300) || (interpolate == TRUE)) {
        kx <- Nti+maxL + 4 + 2 # number of knots for X(t)
        # kx <- Nti+maxL + 4 # number of knots for X(t)
        interpolate <- TRUE
      } else{
        kx <- kx.per500 * ifelse(Nti < 500, 1, round(Nti/500)) # number of knots for X(t)
      }

      ### 0.1 Model exposure process
      if(!interpolate) {
        SSx <- mgcv::smoothCon(s(t, bs = "bs", k = kx),
                               absorb.cons = FALSE,
                               data = data.frame(t = t))[[1]] ## reparameterize it later
        knots_x <- SSx$knots
        # X <- SSx$X

        ## sum-to-zero reparameterization for SSx
        QRx <- qr(t(SSx$X) %*% as.vector(rep(1,nrow(SSx$X))))
        Qx <- qr.Q(QRx, complete = TRUE)
        Zx <- Qx[,2:ncol(Qx)]
        ## Check whether the design matrices are identical
        # X_repa <- SSx$X %*% Zx
        # max(unname(model.matrix(xt.fit))[,-1] - SSx$X %*% Zx) # SAME

        xt.fit <- mgcv::gam(x~s(t, bs = "bs", k = kx), data = data.frame(x = x, t = t))
        ## coefficients for X(t)
        alpha_x <- xt.fit$coefficients
      } else {
        # knots_x <- c(rep(t[1]-1-0.2,4), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,4))
        knots_x <- c(t[1]-1-0.3, rep(t[1]-1-0.2,3), t[1]-0.5-0.001, t, t[length(t)]+0.5+0.001, rep(t[length(t)]+1+0.2,3), t[length(t)]+1+0.3)
        # knots_x <- c(rep(t[1]-1.5,3),t[1]-0.99, t, t[length(t)]+0.99, rep(t[length(t)]+1.5,3))
        X <- splines::splineDesign(knots = knots_x, x = knots_x, outer.ok = TRUE)
        ## interpolate = TRUE
        ## set points for boundary and auxiliary boundary
        Xsparse <- as(X, "dgCMatrix")
        alpha_x <- Interpolate(Xsparse, c(rep(0,4), x[1], x, x[length(x)], rep(0,4)))
        # alpha_x <- Interpolate(Xsparse, c(rep(0,4),x,rep(0,4)))
        xt.fit <- "interpolate"
      }


      ### 0.2 Integration
      t <- t[-(1:maxL)] # delete the first maxL days
      x <- x[-(1:maxL)] # delete the first maxL days

      ### integration
      if(!interpolate) {
        integral <- Integral(knots_x, knots_w, kx, kw, maxLreal, Zx, Zwnew, t+0.5, alpha_x, FALSE)
      } else {
        integral <- Integral_interpolate(knots_x, knots_w, kx, kw, maxLreal, Zwnew, t+0.5, alpha_x, FALSE)
      }
      
      ## linear predictor s(t) = f(\int w(l) X(t-l) dl) = f (B_inner(t) %*% alpha_w)
      ## where B_inner(t) = alpha_x %*% D, dim(D) = c(kx, kw), D_{p,q} = \int b_{xp}(t-l)b_{wq}(l) dl.
      ## b_{xp} is the p-th basis function for X(t)
      ## Stack B_inner(t) for t = (1:Nt)+40, we have B_inner. dim(B_inner) = c(Nt, kw)
      ## See Section X.X in Paper ...
      B_inner <- integral$AlphaxD

      ## For constraint 1: \int w(l)^2 dl = 1
      ## equivalent to w^T %*% Dw %*% w = 1
      Dw <- integral$Dw

      E.maxi <- ceiling(max(integral$Xt2))
      E.mini <- floor(min(x))

      return(list(B_inner = B_inner,
                  Dw = Dw,
                  E.max = E.maxi,
                  E.min = E.mini,
                  removed.t = removed.t
                  ))
    })


    B_inner <- do.call("rbind", lapply(DLprepare, "[[", "B_inner"))


    Dw <- DLprepare[[1]]$Dw

    E.max <- do.call("max", lapply(DLprepare, "[[", "E.max"))
    # if(missingArg(E.min))
      # E.min <- -1.0*E.max
    removed.t <- lapply(DLprepare, "[[", "removed.t")

    ## remove the starting time points
    sXdatlist.clean <- mapply(function(sXdati, removed.ti) return(sXdati[-which(sXdati$t %in% removed.ti),]),
                        sXdatlist, removed.t, SIMPLIFY = FALSE)

    return(list(B_inner = B_inner,
                Dw = Dw,
                E.max = E.max,
                removed.t = removed.t,
                sXdatlist.clean = sXdatlist.clean))
  })



  Dw <- sX.DL[[1]]$Dw
  E.max.list <- lapply(sX.DL,"[[", "E.max")
  # if(missingArg(E.min)) E.min.list <- lapply(E.max.list, function(E.max.) -1.0*E.max.)
  E.min.list <- lapply(E.max.list, function(.) 0)
  sXdat <- do.call("rbind", sX.DL[[1]]$sXdatlist.clean)

  B_inner.list <- lapply(sX.DL, "[[", "B_inner")

  # remove NA rows
  na.id <- is.na(sXdat$y)
  sXdat <- sXdat[!na.id, ]
  B_inner <- lapply(B_inner.list, function(xx) xx[!na.id, ])

  M <- length(B_inner) # the number of exposures. But they are additive in this model.
  ## For mixture exposure, see the package acmeDLNM


  ### prepare NCV. construct a list for neighborhood of i
  sXdat$ii <- 1:nrow(sXdat)
  if(!is.null(group_name)) {
    NCVsXdat <- split(as.data.table(sXdat), by = group_name)
  } else {
    NCVsXdat <- list(as.data.table(sXdat))
  }

  nei.list <- lapply(NCVsXdat, function(NCVsXdati) {
    ni <- nrow(NCVsXdati)
    lapply(1:ni, function(ii)
      NCVsXdati$ii[seq(max(ii - kNCV, 1), min(ii + kNCV, ni))]
    )
  })
  nei.list <- do.call("c", nei.list)


  ## 1. time-varying fixed effects
  if(!is.null(fe.varying)) {
    Xlin <- stats::model.matrix(fe.varying,data=sXdat)[,-1,drop=F]
  } else {
    Xlin <- matrix(1, nrow = nrow(sXdat))
  }
  # if(!is.null(unpen.smooth)) {
  #   unpen.smooth <- lapply(lapply(attr(terms(unpen.smooth),"term.labels"), function(text) parse(text = text)), eval)
  #   unpen.SS <- lapply(lapply(unpen.smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct unpenalized smooth terms
  #   unpen.Xlist <- lapply(unpen.SS,'[[','X')
  #   unpen.X <- Reduce(cbind,unpen.Xlist) ## design matrix
  #   Xlin <- cbind(Xlin, unpen.X)
  # }

  ### code following https://github.com/awstringer1/mam/blob/master/R/mam.R
  ## start following ##
  ## 2.2 smooth term
  numsmooth <- 0 # initialize
  if(!is.null(smooth)){
    smooth <- lapply(lapply(attr(terms(smooth),"term.labels"), function(text) parse(text = text)), eval)

    SS <- lapply(lapply(smooth,mgcv::smoothCon,data=sXdat,absorb.cons = TRUE),'[[',1) ## construct smooth terms
    numsmooth <- length(smooth) # Number of smooth terms

    EE <- lapply(lapply(lapply(SS,'[[','S'),'[[',1),eigen) ## eigen decomposition for penalty matrix

    p <- sapply(lapply(EE,'[[','vectors'),ncol) ## dimension of penalty matrix
    r <- sapply(lapply(EE,'[[','values'),function(x) sum(x>.Machine$double.eps)) ## rank of penalty matrix
    m <- p-r ## dim of null space (minus intercept)
    URlist <- mapply(function(x,y) x[ ,1:y],lapply(EE,'[[','vectors'),r,SIMPLIFY = FALSE)
    UFlist <- mapply(
      function(x,y,z) {
        if (y<z) return(x[ ,(1+y):z])
        newsparsemat(z,z)
      },lapply(EE,'[[','vectors'),r,p,SIMPLIFY = FALSE)
    URlist <- lapply(URlist,cbind) # Ensure they stay matrices
    UFlist <- lapply(UFlist,cbind) # Ensure they stay matrices

    UR <- Matrix::bdiag(URlist)
    UF <- Matrix::bdiag(UFlist)
    # if m=1 UF gets coerced to numeric
    if (!is.matrix(UF)) UF <- cbind(UF)

    Dpi <- Matrix::Diagonal(sum(r),1 / sqrt(Reduce(c,lapply(lapply(EE,'[[','values'),function(x) x[x>.Machine$double.eps]))))

    Xlist <- lapply(SS,'[[','X')
    X <- Reduce(cbind,Xlist) ## design matrix

    Xrand <- as.matrix(X %*% UR %*% Dpi) ## reparametrized
    Xfix <- as.matrix(X %*% UF)
    dups <- !base::duplicated(t(Xfix)) & apply(Xfix,2,function(x) !all(x==0)) # Remove the duplicated intercepts
    if (length(dups) > 1) Xfix <- Xfix[ ,which(dups)]

    model.choice = "with.smooth"
  } else{
    model.choice = "without.smooth"
  }


  # 2.3 add the intercept
  if(exists("Xfix")) {
    Xfix <- cbind(Xlin,Xfix) ## linear effect + unpenalized columns
  } else {
    Xfix <- Xlin
  }
  if (any(apply(Xfix,2,function(x) all(x==1)))) Xfix <- Xfix[,-(apply(Xfix,2,function(x) all(x==1)))]

  if (!is.null(fe.cont)){
    Xgroup <- stats::model.matrix(fe.cont,data=sXdat)[,-1,drop=FALSE]
    if (any(apply(Xgroup,2,function(x) all(x==1)))) Xgroup <- Xgroup[,-(apply(Xgroup,2,function(x) all(x==1)))] # remove the intercept
    Xfix <- cbind(1, Xgroup, Xfix) # add the intercept to the first column
  } else {
    Xfix <- cbind(1,Xfix) # add the intercept to the first column
  }
  ### END following https://github.com/awstringer1/mam/blob/master/R/mam.R #######


  Xoffset <- object$data$offset$Xoffset


  
  
  N <- nrow(sXdat)
  y <- sXdat$y
  t <- sXdat$t


  ## f functions
  SSf.list <- object$data$SSf.list

  eta1 <- apply(sapply(1:M, function(i) {
    E <- B_inner[[i]] %*% object$point$alpha_w.list[[i]]
    SSfi <- SSf.list[[i]]
    Bf <- sapply(E, function(Ei) Bsplinevec2Con(Ei, SSfi$SSf$knots, 4, SSfi$Zfnew))
    return(as.vector(t(Bf) %*% object$point$alpha_f.list[[i]]))
  }),
  1, sum
  )

  eta2 <- as.vector(Xfix %*% object$point$betaF)


  if(exists("Xrand")) eta2 <- eta2 + as.vector(Xrand %*% object$point$betaR)

  eta.est <- eta1 + eta2 + Xoffset

  mod.build <- aceDLNMadditivebuild(y, B_inner, lapply(lapply(SSf.list, "[[", "SSf"), "[[", "knots"),
                                  SwI_large, lapply(SSf.list, "[[", "SfI"), Dw,
                                  Xrand, Xfix, lapply(SSf.list, "[[", "Zfnew"), Xoffset, r,
                                  K,a,
                                  object$point$alpha_f,
                                  object$point$phi,
                                  object$point$log_theta, object$point$log_smoothing_f, object$point$log_smoothing_w,
                                  object$point$betaR, object$point$betaF, object$point$log_smoothing)

  mod.address <- mod.build$address.eigen

  NCVresults <- NCVaceDLNMadditive(mod.address, nei.list, verbose, NCV.nthreads)

  return(NCVresults)
}
