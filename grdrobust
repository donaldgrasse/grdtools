grdrobust <- function(data, run_var, out, treat, adjust, covs, theta_start, crs, p, lat_long, user_bw, input_bw){
  
  packs <- c('intamap', 'sp', 'dplyr', 'purrr', 'sf', 'rdrobust', 'fields')
  suppressMessages(lapply(packs, require, character.only = T))
  
  mKrig <- function(x, y, weights=rep(1, nrow(x)), Z = NULL,
                    cov.function="stationary.cov", 
                    cov.args = NULL, lambda = 0, m = 2, 
                    chol.args = NULL, find.trA = TRUE, NtrA = 20, 
                    iseed = 123, llambda = NULL, na.rm=FALSE, 
                    collapseFixedEffect = TRUE, ...) {
    # pull extra covariance arguments from ...  and overwrite
    # any arguments already named in cov.args
    ind<- match( names( cov.args), names(list(...) ) )
    cov.args = c(cov.args[is.na(ind)], list(...))
    #
    #If cov.args$find.trA is true, set onlyUpper to FALSE (onlyUpper doesn't
    #play nice with predict.mKrig, called by mKrig.trace)
    #
    if(find.trA == TRUE && supportsArg(cov.function, "onlyUpper"))
      cov.args$onlyUpper= FALSE
    if(find.trA == TRUE && supportsArg(cov.function, "distMat"))
      cov.args$distMat= NA
    
    if (!is.null(llambda)) {
      lambda <- exp(llambda)
    }
    # see comments in Krig.engine.fixed for algorithmic commentary
    #
    # check for duplicate x's.
    # stop if there are any
    if (any(duplicated(cat.matrix(x)))) {
      stop("locations are not unique see help(mKrig) ")
    }
    # next function also omits NAs from x,y,weights, and Z  if na.rm=TRUE.
    object<- mKrigCheckXY( x, y, weights, Z, na.rm = na.rm)
    # create fixed part of model as m-1 order polynomial
    # NOTE: if m==0 then fields.mkpoly returns a NULL to 
    # indicate no polynomial part.
    Tmatrix <- cbind(fields.mkpoly(object$x, m), object$Z)
    # set some dimensions
    np <- nrow(object$x)
    if( is.null(Tmatrix) ){
      nt<- 0
    }
    else{
      nt<- ncol(Tmatrix) 
    }
    if( is.null(object$Z)){
      nZ<- 0
    }
    else{
      nZ<- ncol(object$Z)
    }
    ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ)) 
    # as a place holder for reduced rank Kriging, distinguish between
    # observations locations and  the locations to evaluate covariance.
    # (this is will also allow predict.mKrig to handle a Krig object)
    object$knots <- object$x
    # covariance matrix at observation locations
    # NOTE: if cov.function is a sparse constuct then Mc will be sparse.
    # see e.g. wendland.cov
    Mc <- do.call(cov.function, c(cov.args, list(x1 = object$knots, x2 = object$knots)))
    #
    # decide how to handle the pivoting.
    # one wants to do pivoting if the matrix is sparse.
    # if Mc is not a matrix assume that it is in sparse format.
    #
    sparse.flag <- !is.matrix(Mc)
    #
    # set arguments that are passed to cholesky
    #
    if (is.null(chol.args)) {
      chol.args <- list(pivot = sparse.flag)
    }
    else {
      chol.args <- chol.args
    }
    # quantify sparsity of Mc for the mKrig object
    nzero <- ifelse(sparse.flag, length(Mc@entries), np^2)
    # add diagonal matrix that is the observation error Variance
    # NOTE: diag must be a overloaded function to handle sparse format.
    if (lambda != 0) {
      if(! sparse.flag)
        invisible(.Call("addToDiagC", Mc, as.double(lambda/object$weights), nrow(Mc), PACKAGE="fields")
        )
      else
        diag(Mc) = diag(Mc) + lambda/object$weights
    }
    #  MARK LINE Mc
    # At this point Mc is proportional to the covariance matrix of the
    # observation vector, y.
    #
    # cholesky decoposition of Mc
    # do.call used to supply other arguments to the function
    # especially for sparse applications.
    # If chol.args is NULL then this is the same as
    #              Mc<-chol(Mc), chol.args))
    Mc <- do.call("chol", c(list(x = Mc), chol.args))
    
    lnDetCov <- 2 * sum(log(diag(Mc)))
    
    #
    # start linear algebra to find estimates and likelihood
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    if( !is.null(Tmatrix)){
      # Efficent way to multply inverse of Mc times the Tmatrix  
      VT <- forwardsolve(Mc, x = Tmatrix, k=ncol(Mc), transpose = TRUE, upper.tri = TRUE)
      qr.VT <- qr(VT)
      
      # now do generalized least squares for d
      d.coef <- as.matrix(qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
                                                      object$y, upper.tri = TRUE)))
      
      if (collapseFixedEffect) {
        # use a common estimate of fixed effects across all replicates      
        d.coefMeans <- rowMeans(d.coef)
        d.coef <- matrix(d.coefMeans, ncol = ncol(d.coef), 
                         nrow = nrow(d.coef))
      }
      
      resid<-  object$y - Tmatrix %*% d.coef
      # GLS covariance matrix for fixed part.
      Rinv <- solve(qr.R(qr.VT))
      Omega <- Rinv %*% t(Rinv)
      #    
      #  Omega is  solve(t(Tmatrix)%*%solve( Sigma)%*%Tmatrix)
      #   where Sigma = cov.function( x,x) + lambda/object$weights    
      #   proportional to fixed effects covariance matrix.
      #    for the GLS estimates of
      #  the fixed linear part of the model. 
      #    
      #  SEdcoef = diag( Omega) * rho.MLE.FULL
      # 
      # if fixed effects are pooled across replicate fields then
      # adjust the Omega matrix to reflect a mean estimate.
      if (collapseFixedEffect) {
        Omega <- Omega/ ncol(d.coef)
      }
      
      R2diag<-  diag( qr.R(qr.VT) )^2
      lnDetOmega<- -1* sum( log(R2diag) ) 
    }
    else{
      # much is set to NULL because no fixed part of model    
      nt<- 0
      resid<- object$y
      Rinv<- NULL
      Omega<- NULL
      qr.VT<- NULL
      d.coef<- NULL
      lnDetOmega <- 0
    }
    # and now find c.
    #  the coefficents for the spatial part.
    # if linear fixed part included resid as the residuals from the 
    # GLS regression.
    c.coef <- as.matrix(forwardsolve(Mc, transpose = TRUE,
                                     resid, upper.tri = TRUE))
    # save intermediate result this is   t(y- T d.coef)( M^{-1}) ( y- T d.coef)
    quad.form <- c(colSums(as.matrix(c.coef^2)))
    # find c coefficients
    c.coef <- as.matrix(backsolve(Mc, c.coef))
    # MLE estimate of rho and sigma
    #    rhohat <- c(colSums(as.matrix(c.coef * y)))/(np - nt)
    # NOTE if y is a matrix then each of these are vectors of parameters.
    rho.MLE <- quad.form/np
    rhohat <- c(colSums(as.matrix(c.coef * object$y)))/np
    shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
    # the  log profile likehood with  rhohat  and  dhat substituted
    # leaving a profile for just lambda.
    # NOTE if y is a matrix then this is a vector of log profile
    # likelihood values.
    lnProfileLike <- (-np/2 - log(2 * pi) * (np/2) - (np/2) * 
                        log(rho.MLE) - (1/2) * lnDetCov)
    # see section 4.2 handbook of spatial statistics (Zimmermanchapter)
    lnProfileREML <-  lnProfileLike + (1/2) * lnDetOmega
    rho.MLE.FULL <- mean(rho.MLE)
    sigma.MLE.FULL <- sqrt(lambda * rho.MLE.FULL)
    # if y is a matrix then compute the combined likelihood
    # under the assumption that the columns of y are replicated
    # fields
    lnProfileLike.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
                                 (np/2) * log(rho.MLE.FULL) 
                               - (1/2) * lnDetCov)
    )
    lnProfileREML.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
                                 (np/2) * log(rho.MLE.FULL) 
                               - (1/2) * lnDetCov
                               + (1/2) * lnDetOmega  )
    )
    
    #
    # return coefficients and   include lambda as a check because
    # results are meaningless for other values of lambda
    # returned list is an 'object' of class mKrig (micro Krig)
    # also save the matrix decompositions so coefficients can be
    # recalculated for new y values.  Make sure onlyUpper and 
    # distMat are unset for compatibility with mKrig S3 functions
    if(!is.null(cov.args$onlyUpper))
      cov.args$onlyUpper = FALSE
    if(!is.null(cov.args$distMat))
      cov.args$distMat = NA
    object <- c( object, list( 
      d = d.coef, c = c.coef, nt = nt, np = np, 
      lambda.fixed = lambda, 
      cov.function.name = cov.function, 
      args = cov.args, m = m, chol.args = chol.args, call = match.call(), 
      nonzero.entries = nzero, shat.MLE = sigma.MLE, sigma.MLE = sigma.MLE, 
      rho.MLE = rho.MLE, rhohat = rho.MLE, lnProfileLike = lnProfileLike, 
      rho.MLE.FULL = rho.MLE.FULL, sigma.MLE.FULL = sigma.MLE.FULL, 
      lnProfileLike.FULL = lnProfileLike.FULL,
      lnProfileREML.FULL =  lnProfileREML.FULL,
      lnProfileREML =  lnProfileREML,
      lnDetCov = lnDetCov, lnDetOmega = lnDetOmega,
      quad.form = quad.form, Omega = Omega,lnDetOmega=lnDetOmega,
      qr.VT = qr.VT, 
      Mc = Mc,
      Tmatrix = Tmatrix, ind.drift = ind.drift, nZ = nZ,
      fixedEffectsCov = Omega * rho.MLE.FULL, 
      # dcoefSE = sqrt(diag( Omega) * rho.MLE.FULL),  
      collapseFixedEffect= collapseFixedEffect)
    )
    #
    # find the residuals directly from solution
    # to avoid a call to predict
    object$residuals <- lambda * c.coef/object$weights
    object$fitted.values <- object$y - object$residuals
    # estimate effective degrees of freedom using Monte Carlo trace method.
    if (find.trA) {
      object2 <- mKrig.trace(object, iseed, NtrA)
      object$eff.df <- object2$eff.df
      object$trA.info <- object2$trA.info
      object$GCV <- (sum(object$residuals^2)/np)/(1 - object2$eff.df/np)^2
      if (NtrA < np) {
        object$GCV.info <- (sum(object$residuals^2)/np)/(1 - object2$trA.info/np)^2
      }
      else {
        object$GCV.info <- NA
      }
    }
    else {
      object$eff.df <- NA
      object$trA.info <- NA
      object$GCV <- NA
    }
    class(object) <- "mKrig"
    return(object)
  }
  "rdist.earth" <- function(x1, x2=NULL, miles = TRUE, R = NULL) {
    if (is.null(R)) {
      if (miles) 
        R <- 3963.34
      else R <- 6378.388
    }
    coslat1 <- cos((x1[, 2] * pi)/180)
    sinlat1 <- sin((x1[, 2] * pi)/180)
    coslon1 <- cos((x1[, 1] * pi)/180)
    sinlon1 <- sin((x1[, 1] * pi)/180)
    if (is.null(x2)) {
      pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
        t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
      return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
    }
    else {
      coslat2 <- cos((x2[, 2] * pi)/180)
      sinlat2 <- sin((x2[, 2] * pi)/180)
      coslon2 <- cos((x2[, 1] * pi)/180)
      sinlon2 <- sin((x2[, 1] * pi)/180)
      pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*% 
        t(cbind(coslat2 * coslon2, coslat2 * sinlon2, sinlat2))
      return(R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp)))
    }
  }
  
  
  
  
  
  spatialProcess <- function(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                             mKrig.args = list( m=2),
                             cov.function = "stationary.cov", 
                             cov.args = list(Covariance = "Matern", smoothness = 1),
                             theta = NULL, 
                             theta.start = NULL,
                             lambda.start = .5, 
                             theta.range = NULL, 
                             abstol = 1e-4,
                             na.rm = TRUE,
                             verbose = FALSE,
                             REML = FALSE, 
                             ...) {
    
    # NOTE all ... information is assumed to be for the cov.args list
    # overwrite the default choices (some R arcania!)
    ind<- match( names( cov.args), names(list(...) ) )
    cov.args <- c(cov.args[is.na(ind)], list(...))
    if( verbose){
      cat("extra arguments from ... " , names( list(...)) , fill=TRUE)
      cat(" full list from cov.args: ", names(cov.args) )
    }
    
    # NOTE: switch to find theta MLE    is.null( theta)
    if( !is.null(theta)){
      par.grid<- list( theta = theta)
      MLEInfo <-mKrigMLEGrid(x, y,  weights = weights, Z= Z, 
                             mKrig.args = mKrig.args,
                             cov.fun = cov.function, 
                             cov.args  = cov.args,
                             par.grid = par.grid, 
                             lambda = lambda.start, 
                             lambda.profile = TRUE, 
                             na.rm = na.rm,
                             verbose = FALSE,
                             REML = REML) 
      lambda.MLE <- MLEInfo$lambda.MLE
      theta.MLE <- NA
      theta.95CI<- NA
      thetaModel <- theta
      if( verbose){
        print( MLEInfo$summary)
      }
    }
    else{
      #  
      # NOTE MLEspatialProcess omits NAs
      MLEInfo <- MLESpatialProcess(x, y, weights = weights, Z=Z, 
                                   mKrig.args = mKrig.args,
                                   cov.function = cov.function, 
                                   cov.args = cov.args,
                                   
                                   theta.start = theta.start, 
                                   theta.range = theta.range, 
                                   gridN = 20,
                                   lambda.start = lambda.start,
                                   abstol = abstol,
                                   verbose = FALSE,
                                   REML = REML
      )
      lambda.MLE <- MLEInfo$MLEJoint$pars.MLE[1] 
      theta.MLE<- MLEInfo$MLEJoint$pars.MLE[2]
      thetaModel<- theta.MLE
      # approximate confidence interval for theta 
      thetaGrid<- MLEInfo$MLEGrid$par.grid$theta
      lgProfileLike<- MLEInfo$MLEGrid$summary[,2]
      splineFit<- splint(thetaGrid, lgProfileLike, nx=500)
      cutLevel<- max(splineFit$y ) - qchisq(.95, 1) / 2
      ind<- splineFit$y> cutLevel
      lower<- min( splineFit$x[ ind] )
      upper<- max(splineFit$x[ ind])
      theta.95CI = c( lower, upper)
    }
    #  
    if( verbose){
      cat("Summary from joint optimization", fill=TRUE)
      print( MLEInfo$MLEJoint$summary )
      print( MLEInfo$MLEJoint$pars.MLE)
      print(theta.95CI )
    }
    # now fit spatial model with MLE for theta (range parameter)
    #  or the value supplied in the call
    # reestimate the other parameters for simplicity to get the complete mKrig object
    obj <- do.call( "mKrig", 
                    c( list(x=x,
                            y=y,
                            weights=weights,
                            Z=Z),
                       mKrig.args,
                       list( na.rm=na.rm),
                       list(           
                         cov.function = cov.function,
                         cov.args = cov.args, 
                         lambda = lambda.MLE,
                         theta  = thetaModel
                       )
                    )
    )
    obj <- c(obj, list(   MLEInfo = MLEInfo,
                          thetaModel = thetaModel,
                          theta.MLE = theta.MLE,
                          theta.95CI = theta.95CI,
                          lambda.MLE = lambda.MLE,
                          summary = MLEInfo$summary)
    )
    # replace call to mKrig with this top level one
    obj$call<- match.call()	
    class(obj) <- c( "spatialProcess","mKrig")
    
    return(obj)
  }
  
  "Matern" <- function(d, range = 1, alpha = 1/range, 
                       smoothness = 0.5, nu = smoothness, phi = 1.0) {
    #
    # Matern covariance function transcribed from Stein's book page 31
    # nu==smoothness, alpha ==  1/range
    #
    # GeoR parameters map to kappa==smoothness and phi == range
    # check for negative distances
    # phi is accepted as the marginal variance of the process (see below)
    # within fields, however, this parameter is "rho" and we recommend
    # not using phi.
    
    if (any(d < 0)) 
      stop("distance argument must be nonnegative")
    d <- d * alpha
    # avoid sending exact zeroes to besselK
    d[d == 0] <- 1e-10
    #
    # the hairy constant ...
    con <- (2^(nu - 1)) * gamma(nu)
    con <- 1/con
    #
    # call to  Bessel function from R base package
    #
    return(phi * con * (d^nu) * besselK(d, nu))
  }
  
  mKrigMLEGrid <- function(x, y, weights = rep(1, nrow(x)), Z = NULL,
                           mKrig.args = NULL,
                           cov.fun = "stationary.cov", 
                           cov.args = NULL,
                           na.rm = TRUE, 
                           par.grid = NULL, 
                           lambda = NULL, 
                           lambda.profile = TRUE, 
                           relative.tolerance = 1e-04,
                           REML = FALSE,
                           verbose = FALSE) {
    if( na.rm){
      obj<- mKrigCheckXY(x, y, weights, Z, na.rm)
      x<- obj$x
      y<- obj$y
      weights<- obj$weights
      Z<- obj$Z
    }
    #check which optimization options the covariance function supports
    #precompute distance matrix if possible so it only needs to be computed once
    supportsDistMat = supportsArg(cov.fun, "distMat")
    #precompute distance matrix if possible so it only needs to be computed once
    if(supportsDistMat) {
      #Get distance function and arguments if available
      #If user left all distance settings NULL, use rdist with compact option.
      #Use rdist function by default in general.
      #
      if(is.null(cov.args$Distance)) {
        cov.args$Distance  <-  "rdist"
        cov.args$Dist.args <- list(compact=TRUE)
      }
      cov.args$distMat<-do.call(cov.args$Distance, c( list(x), cov.args$Dist.args) )
      cov.args$onlyUpper<- TRUE
    }
    
    lnProfileLike.max <- -1e+20
    # find NG --  number of parameters to try
    par.grid <- data.frame(par.grid)
    if (nrow(par.grid) == 0) {
      NG<- ifelse(is.null(lambda), 1, length( lambda)) 
    }
    else {
      NG <- nrow(par.grid)
    }
    lambda.best <- NA
    # default for lambda is 1.0 for first value and exp(llambda.opt) for subsequent ones
    # this is controlled by NAs for lambda starting values.
    if (is.null(lambda)) {
      lambda <- rep(NA, NG)
    }
    # output matrix to summarize results
    summary <- matrix(NA, nrow = NG, ncol = 8)
    
    # default starting value for lambda is .5 or log lambda is 0
    lambda.opt <- .5
    optim.counts <- c(NA, NA)
    lnLike.eval <- list()
    # Define the objective function as a tricksy call to mKrig
    # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
    #
    # begin loop over covariance arguments
    lnLike.eval<- list()
    for (k in 1:NG) {
      lambda.start <- ifelse(is.na(lambda[k]), lambda.opt, (lambda[k]))
      # list of covariance arguments from par.grid with right names (some R arcania!)
      # note that this only works because 1) temp.fn will search in this frame for this object
      # par.grid has been coerced to a data frame so one has a concept of a row subscript.
      cov.args.temp <- as.list(par.grid[k, ])
      names(cov.args.temp) <- names(par.grid)
      currentCov.args<- c(cov.args.temp, cov.args) 
      # optimize over lambda if lambda.profile is TRUE
      optim.args = list(method = "BFGS", 
                        control = list(fnscale = -1, parscale = c(0.5), 
                                       ndeps = c(0.05)))
      if (lambda.profile) {
        # set up matrix to store evaluations from within optim
        MLEfit0 <- mKrigMLEJoint(x, y, weights=weights, Z=Z, 
                                 lambda.start = lambda.start, 
                                 cov.params.start = NULL, 
                                 cov.fun = cov.fun,
                                 optim.args = optim.args,
                                 cov.args = currentCov.args,
                                 na.rm = na.rm,
                                 mKrig.args = mKrig.args,
                                 REML = REML,
                                 verbose = verbose)
        lnLike.eval<- c( lnLike.eval, list(MLEfit0$lnLike.eval))
        lambda.opt<- MLEfit0$pars.MLE[1]
      }
      else {
        # no refinement for lambda so just save the the 'start' value as final one.
        lambda.opt <- lambda.start
      }
      
      # final fit at optimal value 
      #    (or starting value if not refinement/maximization for lambda)
      obj <- do.call("mKrig", c(
        list(x = x, y = y, weights = weights, Z = Z, na.rm = na.rm),
        mKrig.args,
        list(lambda=lambda.opt),
        list( cov.fun= cov.fun, cov.args = currentCov.args)
      )
      )
      nameCriterion<- ifelse( !REML,
                              "lnProfileLike.FULL",
                              "lnProfileREML.FULL" )
      if (obj[[nameCriterion]] > lnProfileLike.max) {
        lnProfileLike.max <- obj$lnProfileLike.FULL
        cov.args.MLE <- cov.args.temp
        lambda.best <- lambda.opt
      }
      
      # save results of the kth covariance model evaluation
      summary[k, 1:8] <- c(obj$eff.df, obj[[nameCriterion]], 
                           obj$GCV, obj$sigma.MLE.FULL, obj$rho.MLE.FULL, lambda.opt, 
                           optim.counts)
      dimnames(summary) <- list(NULL, c("EffDf",nameCriterion , 
                                        "GCV", "sigma.MLE", "rho.MLE", "lambda.MLE", "counts eval", 
                                        "counts grad"))
      if (verbose) {
        cat("Summary: ", k, summary[k, 1:8], fill = TRUE)
      }
    }
    return(list(summary = summary, par.grid = par.grid, cov.args.MLE = cov.args.MLE, 
                lambda.best = lambda.best, lambda.MLE = lambda.best, 
                call = match.call(), lnLike.eval = lnLike.eval)
    )
  }
  
  
  mKrig.MLE.joint <- function(x, y, weights = rep(1, nrow(x)), 
                              lambda.guess = 1, cov.params.guess=NULL, 
                              cov.fun="stationary.cov", cov.args=NULL, 
                              Z = NULL, optim.args=NULL, find.trA.MLE = FALSE, 
                              ..., verbose = FALSE) {
    
    #set default optim.args if necessary
    if(is.null(optim.args))
      optim.args = list(method = "BFGS", 
                        control=list(fnscale = -1, 
                                     ndeps = rep(log(1.1), length(cov.params.guess)+1), 
                                     reltol=1e-04, maxit=10))
    
    #check which optimization options the covariance function supports
    supportsDistMat = supportsArg(cov.fun, "distMat")
    
    #precompute distance matrix if possible so it only needs to be computed once
    if(supportsDistMat) {
      
      #Get distance function and arguments if available
      #
      Dist.fun= c(cov.args, list(...))$Distance
      Dist.args=c(cov.args, list(...))$Dist.args
      
      #If user left all distance settings NULL, use rdist with compact option.
      #Use rdist function by default in general.
      #
      if(is.null(Dist.fun)) {
        Dist.fun = "rdist"
        if(is.null(Dist.args))
          Dist.args = list(compact=TRUE)
      }
      
      distMat = do.call(Dist.fun, c(list(x), Dist.args))
    }
    
    #set cov.args for optimal performance if possible
    if(supportsDistMat)
      cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
    
    # these are all the arguments needed to call mKrig except lambda and cov.args
    mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, cov.fun=cov.fun), 
                    list(...))
    mKrig.args$find.trA = find.trA.MLE
    
    # output matrix to summarize results
    ncolSummary = 8 + length(cov.params.guess)
    summary <- matrix(NA, nrow = 1, ncol = ncolSummary)
    dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", "GCV", "sigma.MLE", 
                                      "rho.MLE", "llambda.MLE", names(cov.params.guess), 
                                      "counts eval","counts grad"))
    
    # Define the objective function as a tricksy call to mKrig
    # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
    lnProfileLike.max <- -Inf
    temp.fn <- function(parameters) {
      # Separate lambda from covariance parameters.
      # Optimization is over log-scale so exponentiate log-parameters.
      lambda = exp(parameters[1])
      if(length(parameters) > 1) {
        otherParams = as.list(exp(parameters[2:length(parameters)]))
        names(otherParams) = names(cov.params.guess)
      }
      else
        otherParams = NULL
      
      #get all this eval's covariance arguments using the input parameters
      cov.args.temp = c(cov.args, otherParams)
      
      # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
      # assign to hold the last mKrig object
      hold <- do.call("mKrig", c(mKrig.args, list(lambda = lambda),
                                 cov.args.temp))
      
      #save best mKrig object to global environment
      if(hold$lnProfileLike.FULL > lnProfileLike.max) {
        out <<- hold
        lnProfileLike.max = hold$lnProfileLike.FULL
      }
      hold = hold[c("rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
      
      # add this evalution to an object (i.e. here a matrix) in the calling frame
      temp.eval <- get("capture.evaluations")
      assign("capture.evaluations", rbind(temp.eval, c(parameters, unlist(hold))), envir = capture.env)
      return(hold$lnProfileLike.FULL)
    }
    
    #
    # optimize over covariance parameters and lambda
    
    # list of covariance arguments from par.grid with right names (some R arcania!)
    # note that this only works because 1) temp.fn will search in this frame for this object
    # par.grid has been coerced to a data frame so one has a concept of a row subscript.
    
    # set up matrix to store evaluations from within optim
    capture.evaluations <- matrix(NA, ncol = 4+length(cov.params.guess), nrow = 1,
                                  dimnames = list(NULL, c("lambda", names(cov.params.guess), "rho.MLE",
                                                          "sigma.MLE", "lnProfileLike.FULL")))
    capture.env <- environment()
    
    # call to optim with initial guess (on log-scale)
    init.guess = log(unlist(c(lambda.guess, cov.params.guess)))
    look <- do.call(optim, c(list(par=init.guess), list(temp.fn), optim.args))
    
    #get optim results
    optim.counts <- look$counts
    llambda.opt <- look$par[1]
    lambda.opt <- exp(llambda.opt)
    if(length(look$par) > 1) {
      params.opt <- exp(look$par[2:length(look$par)])
      params.opt <- as.list(params.opt)
      names(params.opt) <- names(cov.params.guess)
    }
    else
      params.opt=NULL
    
    # call to 1-d search
    #            opt.summary     <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
    #            llambda.opt <- opt.summary$maximum
    #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
    # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
    lnLik.eval <- capture.evaluations[-1,]
    
    #exponentiate lambda and covariance parameters in lnLik.eval
    lnLik.eval[, 1:length(look$par)] = exp(lnLik.eval[, 1:length(look$par)])
    
    # calculate trace of best mKrig object if necessary
    #
    find.trA = list(...)$find.trA
    if(is.null(find.trA) || find.trA) {
      
      #get arguments for mKrig.trace
      iseed = list(...)$iseed
      NtrA = list(...)$NtrA
      
      #set iseed and NtrA to default values of mKrig if NULL
      if(is.null(iseed))
        iseed = 123
      if(is.null(NtrA))
        NtrA = 20
      
      #update best mKrig object with trace results
      out2 <- mKrig.trace(out, iseed, NtrA)
      out$eff.df <- out2$eff.df
      out$trA.info <- out2$trA.info
      np <- nrow(x)
      out$GCV <- (sum(out$residuals^2)/np)/(1 - out2$eff.df/np)^2
      if (NtrA < np)
        out$GCV.info <- (sum(out$residuals^2)/np)/(1 - out2$trA.info/np)^2
      else
        out$GCV.info <- NA
    }
    
    # save results of the best covariance model evaluation in a neat table
    
    summary[1, 1:ncolSummary] <- unlist(c(out$eff.df, out$lnProfileLike.FULL, 
                                          out$GCV, out$sigma.MLE.FULL, out$rho.MLE.FULL,
                                          llambda.opt, 
                                          params.opt, optim.counts))
    
    if (verbose) {
      cat("Summary: ", 1, summary[1, 1:ncolSummary], fill = TRUE)
    }
    
    #add summary table to output mKrig object and ensure it is still 
    #of class mKrig
    out = c(out, list(summary=summary, lnLik.eval=lnLik.eval))
    class(out) = "mKrig"
    
    return(out)
  }
  
  MLESpatialProcess <- function(x, y, weights = rep(1, nrow(x)), Z = NULL,
                                mKrig.args = NULL,
                                cov.function = "stationary.cov", 
                                cov.args = list(Covariance = "Matern",
                                                smoothness = 1), 
                                lambda.start = .5,
                                theta.start = NULL, 
                                theta.range = NULL,
                                gridN = 20,
                                optim.args = NULL,
                                na.rm = TRUE,
                                verbose = FALSE,
                                abstol  = 1e-4,
                                REML = FALSE,
                                ...) {
    if( verbose){
      cat(" MLESpatialProcess extra arguments:" , full=TRUE)
      print( names( list(...)))
    }
    # combine  list(...) with cov.args and omit duplicates but favoring the ... value
    ind<- match( names( cov.args), names(list(...) ) )
    cov.args = c(cov.args[is.na(ind)], list(...))
    
    ########################################################################
    #  evaluate likelihood for a grid of theta on log scale
    # maximizing over lambda.
    #
    # if range or starting value  for range is missing use quantiles of pairwise
    # distances among data.  
    if( is.null( theta.range) ){
      if( is.null( cov.args$Distance)){
        pairwiseD<- dist(x)
      }
      else{
        pairwiseD<- do.call(cov.args$Distance, list(x))
        pairwiseD<- pairwiseD[col(pairwiseD) > row( pairwiseD) ]
      }
      theta.range<- quantile( pairwiseD, c(.02,.97))
    }
    thetaGrid<- 10**seq( log10(theta.range[1]), log10(theta.range[2]), length.out=gridN )
    # 
    par.grid<- list( theta= thetaGrid)
    timeGrid<- system.time(
      MLEGrid<- mKrigMLEGrid(x, y,  weights = weights, Z= Z, 
                             mKrig.args = mKrig.args,
                             cov.fun = cov.function, 
                             cov.args  = cov.args,
                             par.grid = par.grid, 
                             lambda = lambda.start, 
                             lambda.profile = TRUE, 
                             na.rm = na.rm,
                             verbose = verbose,
                             REML = REML)
    )
    ##################################################################################
    #refine MLE for lambda and theta use the best value of theta from grid search if
    # starting value not passed. 
    if ( is.null(theta.start) ) {
      ind<- which.max( MLEGrid$summary[,2] )
      theta.start <-  par.grid$theta[ind]
      lambda.start<- MLEGrid$lambda.best
    }
    timeOptim<- system.time(
      MLEJoint <- mKrigMLEJoint(x, y, weights = weights, Z = Z,
                                mKrig.args = mKrig.args,
                                cov.fun = cov.function,
                                cov.args = cov.args, 
                                lambda.start = lambda.start, 
                                cov.params.start = list(theta=theta.start), 
                                optim.args = optim.args,
                                abstol = abstol,
                                na.rm = na.rm,
                                verbose = verbose,
                                REML = REML)
    )
    #####################################################################################
    # evaluate likelihood on grid of log lambda with MLE for theta
    #NOTE lambda.profile = FALSE makes this work.
    lambdaGrid<-   (10^(seq( -4,4,,gridN)  ))*MLEJoint$pars.MLE[1]
    par.grid<- list( theta= rep(MLEJoint$pars.MLE[2], gridN) )
    if( verbose){ print( par.grid) }
    timeProfile<- system.time(
      MLEProfileLambda <- mKrigMLEGrid(x, y,  weights = weights, Z= Z,
                                       cov.fun = cov.function, 
                                       cov.args  = cov.args,
                                       mKrig.args = mKrig.args,
                                       par.grid = par.grid, 
                                       lambda = lambdaGrid, 
                                       lambda.profile = FALSE, 
                                       na.rm = na.rm,
                                       verbose = verbose,
                                       REML = REML) 
    )
    timingTable<- cbind( timeGrid, timeOptim, timeProfile)
    return(
      list( summary= MLEJoint$summary, MLEGrid= MLEGrid, MLEJoint=MLEJoint, 
            MLEProfileLambda=MLEProfileLambda, call=match.call(), 
            timingTable= timingTable)
    )
  }
  
  mKrigMLEJoint <- function(x, y, weights = rep(1, nrow(x)),  Z = NULL,
                            mKrig.args = NULL,
                            na.rm = TRUE,
                            cov.fun = "stationary.cov", cov.args=NULL, 
                            lambda.start = .5,
                            cov.params.start = NULL,
                            optim.args = NULL,
                            abstol = 1e-4,
                            parTransform = NULL, 
                            REML = FALSE, 
                            verbose = FALSE) {
    # overwrite basic data to remove NAs this has be done in case distance 
    # matrices are precomputed (see below)
    if( na.rm){
      obj<- mKrigCheckXY(x, y, weights, Z, na.rm)
      x<- obj$x
      y<- obj$y
      weights<- obj$weights
      Z<- obj$Z
    }
    #set default optim.args if necessary
    # abstol is anticipating this is a likelihood so differencs of 1e-4 are not appreciable
    # 
    if(is.null(optim.args)){
      optim.args = list(method = "BFGS", 
                        control=list(fnscale = -1,
                                     ndeps = rep(log(1.1),length(cov.params.start)+1), 
                                     abstol = abstol,
                                     maxit = 20)
      )
    }
    # main way to keep track of parameters to optimize -- lambda always included  
    parNames<- c( "lambda", names(cov.params.start))
    if( is.null(parTransform)){
      # parTransform: log/exp
      parTransform<- function( ptemp, inv=FALSE){
        if( !inv){ log( ptemp)}
        else{
          exp(ptemp)
        }
      }
    }
    ########bug
    if(verbose){
      cat("parameters to optimze: ", parNames, fill=TRUE)
    }
    #check which optimization options the covariance function supports
    supportsDistMat = supportsArg(cov.fun, "distMat")
    #precompute distance matrix if possible so it only needs to be computed once
    if(supportsDistMat & is.null( cov.args$distMat)) {
      #Get distance function and arguments if available
      #
      Dist.fun= c(cov.args)$Distance
      Dist.args=c(cov.args)$Dist.args
      
      #If user left all distance settings NULL, use rdist with compact option.
      #Use rdist function by default in general.
      #
      if(is.null(Dist.fun)) {
        Dist.fun = "rdist"
        if(is.null(Dist.args))
          Dist.args = list(compact=TRUE)
      }
      distMat = do.call(Dist.fun, c(list(x), Dist.args))
      #set cov.args for optimal performance
      cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
    }
    # these are all the arguments needed to call mKrig except lambda and cov.args
    mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z),
                    mKrig.args,
                    list(cov.fun=cov.fun) 
    )
    # reset switch so trace is not found for each evaluation of the likelihood.   
    mKrig.args$find.trA = FALSE
    # output matrix to summarize results
    ncolSummary = 7 + length(parNames)
    summary <- matrix(NA, nrow = 1, ncol = ncolSummary)
    dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", "GCV", "sigma.MLE", 
                                      "rho.MLE", parNames, 
                                      "counts eval","counts grad"))
    
    lnProfileLike.max <- -Inf
    
    #
    # optimize over (some) covariance parameters and lambda
    capture.evaluations <- matrix(NA, ncol = length(parNames) + 4 , nrow = 1,
                                  dimnames = list(NULL,
                                                  c(   parNames,
                                                       "rho.MLE",
                                                       "sigma.MLE", 
                                                       "lnProfileLike.FULL",
                                                       "lnProfileREML.FULL")
                                  )
    ) 
    capture.env <- environment()
    # call to optim with initial start (default is log scaling )
    init.start <- parTransform( unlist(c(lambda.start, cov.params.start)), inv=FALSE)
    #  cat("init.start",  init.start, fill=TRUE)
    optimResults <- do.call(optim, c(    list(par=init.start),
                                         list(mKrigJointTemp.fn),
                                         optim.args,
                                         list(  parNames = parNames,
                                                parTransform = parTransform,
                                                mKrig.args = mKrig.args,
                                                cov.args = cov.args, 
                                                capture.env = capture.env,
                                                REML = REML)
    )
    )
    #get optim results
    optim.counts <- optimResults$counts
    parOptimum<- parTransform(optimResults$par, inv=TRUE)
    # first row is just NAs  
    lnLike.eval <- capture.evaluations[-1,]
    
    nameCriterion<- ifelse( !REML,
                            "lnProfileLike.FULL",
                            "lnProfileREML.FULL" )
    ind<- which( lnLike.eval[ , nameCriterion]
                 == optimResults$value )
    ind<- max( ind)
    # below is an aspect from optim I dont understand and thought to flag
    #if( length(ind)!=1 ){
    #     cat( "Weirdness in optimization. See lnLike.eval rows: ", ind,
    #         fill=TRUE )
    # ind<- max( ind)
    #}
    # save results of the best covariance model evaluation in a neat table
    summary <- c(          optimResults$value, 
                           parOptimum,
                           lnLike.eval[ind,"sigma.MLE"],
                           lnLike.eval[ind,"rho.MLE"],
                           optim.counts)
    names(summary) <- c(nameCriterion, parNames, 
                        "sigmaMLE", "rhoMLE", "funEval", "gradEval")
    out = c( list(summary=summary, lnLike.eval = lnLike.eval, optimResults=optimResults,
                  pars.MLE=parOptimum, parTransform = parTransform))
    return(out)
  }
  
  mKrigJointTemp.fn <- function(parameters,
                                mKrig.args, cov.args, parTransform, parNames,
                                REML=FALSE,
                                capture.env) {
    # optimization is over a transformed scale ( so need to back transform for mKrig)
    tPars<- parTransform( parameters, inv=TRUE)
    names( tPars)<- parNames
    #get all this eval's covariance arguments using the input parameters
    cov.args.temp = c(cov.args, tPars)
    # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
    # assign to hold the last mKrig object
    hold <- do.call("mKrig", c(mKrig.args,
                               cov.args.temp))
    
    hold = hold[c("rho.MLE.FULL",
                  "sigma.MLE.FULL",
                  "lnProfileLike.FULL",
                  "lnProfileREML.FULL"
    )]
    
    # add this evalution to an object (i.e. here a matrix) in the calling frame
    temp.eval <- get("capture.evaluations", envir=capture.env)
    assign("capture.evaluations", rbind(temp.eval,
                                        c(parTransform(parameters, inv=TRUE),
                                          unlist(hold))), 
           envir = capture.env)
    if( !REML){
      return(hold$lnProfileLike.FULL)
    }
    else{
      return(hold$lnProfileREML.FULL)
    }  
  }
  
  
  mKrigCheckXY <- function(x, y,  weights, Z, na.rm) 
  {
    #
    # check for missing values in y or X.
    #
    # save logical indicating where there are NA's
    # and check for NA's
    #
    ind <- is.na(y)
    if (any(ind) & !na.rm) {
      stop("Need to remove missing values or use: na.rm=TRUE in the call")
    }
    #
    # coerce x to be a matrix
    x <- as.matrix(x)
    #
    # coerce y to be a matrix
    #
    y <- as.matrix(y)
    #
    #default weights ( reciprocal variance of errors).
    #
    if (is.null(weights)) 
      weights <- rep(1, nrow(y))
    #
    # check that dimensions agree
    #
    if (nrow(y) != nrow(x)) {
      stop(" length of y and number of rows of x differ")
    }
    if (nrow(y) != length(weights)) {
      stop(" length of y and weights differ")
    }
    #  if Z is not NULL coerce to be  a matrix
    # and check  # of rows
    if (!is.null(Z)) {
      if (!is.matrix(Z)) {
        Z <- as.matrix(Z)
      }
      if (nrow(y) != nrow(Z)) {
        stop(" number of rows of y and number of rows of Z differ")
      }
    }
    # if NAs can be removed then remove them and warn the user
    if (na.rm) {
      ind <- is.na(y)
      if(all(ind)){
        stop("Oops! All y values are missing!")
      }
      if (any(ind)) {
        y <- y[!ind]
        x <- as.matrix(x[!ind, ])
        if (!is.null(Z)) {
          Z <- as.matrix(Z[!ind, ])
        }
        weights <- weights[!ind]
      }
    }
    #
    # check for NA's in x matrix -- there should not be any !
    if (any(c(is.na(x)))) {
      stop(" NA's in x matrix")
    }
    #
    # check for NA's in Z matrix
    if (!is.null(Z)) {
      if (any(c(is.na(Z)))) {
        stop(" NA's in Z matrix")
      }
    }
    
    # save x, weights  and y w/o NAs
    N <- length(y)
    return(list(N = N, y = y, x = x, weights = weights, Z = Z, 
                NA.ind = ind) )
  }
  
  supportsArg = function(fun=stationary.cov, arg) {
    
    if(is.null(fun)) {
      #set fun to the default covariance function if not specified
      fun = stationary.cov
    }
    
    argNames = names(as.list(args(fun)))
    return(any(argNames == arg))
  }
  "stationary.cov" <- function(x1, x2=NULL, Covariance = "Exponential", Distance = "rdist", 
                               Dist.args = NULL, theta = 1, V = NULL, C = NA, marginal = FALSE, 
                               derivative = 0, distMat = NA, onlyUpper = FALSE, ...) {
    
    # get covariance function arguments from call
    cov.args <- list(...)
    # coerce x1 and x2 to matrices
    if (is.data.frame(x1)) 
      x1 <- as.matrix(x1)
    if (!is.matrix(x1)) 
      x1 <- matrix(c(x1), ncol = 1)
    
    if(is.null(x2))
      x2 = x1
    if (is.data.frame(x2)) 
      x2 <- as.matrix(x2)
    if (!is.matrix(x2)& !is.null(x2)) 
      x2 <- matrix(c(x2), ncol = 1)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    #
    # separate out a single scalar transformation and a
    # more complicated scaling and rotation.
    # this is done partly to have the use of great circle distance make sense
    # by applying the scaling  _after_ finding the distance.
    #
    if (length(theta) > 1) {
      stop("theta as a vector matrix has been depreciated use the V argument")
    }
    #
    # following now treats V as a full matrix for scaling and rotation.
    #
    # try to catch incorrect conbination  of great circle distance and V
    if (Distance == "rdist.earth" & !is.null(V)) {
      stop("V not supported with great circle distance")
    }
    # use V the anisotropic scaling and rotation from covariance arguments 
    # if exists
    if( !is.null(cov.args$V) ){
      V<- cov.args$V
    }
    if (!is.null(V)) {
      if (theta != 1) {
        stop("can't specify both theta and V!")
      }
      x1 <- x1 %*% t(solve(V))
      x2 <- x2 %*% t(solve(V))
    }
    #
    # locations are now scaled and rotated correctly
    # now apply covariance function to pairwise distance matrix, or multiply
    # by C vector or just find marginal variance
    # this if block finds the cross covariance matrix
    if (is.na(C[1]) && !marginal) {
      #
      # if distMat is supplied, evaluate covariance for upper triangular part only
      #
      if(is.na(distMat[1])) {
        # distMat not supplied so must compute it along with covariance matrix
        # note overall scaling by theta (which is just theta under isotropic case)
        if(is.null(x2))
          distMat <- do.call(Distance, c(list(x1), Dist.args))
        else
          distMat <- do.call(Distance, c(list(x1=x1, x2=x2), Dist.args))
        
      }
      
      #
      # now convert distance matrix to covariance matrix
      #
      if(inherits(distMat, "dist")) {
        #
        # distMat is in compact form, so evaluate covariance over all distMat and convert to matrix form
        
        diagVal = do.call(Covariance, c(list(d=0), cov.args))
        
        if(onlyUpper)
          return(compactToMat(do.call(Covariance, c(list(d=distMat*(1/theta)), cov.args)), diagVal))
        else
          # if onlyUpper==FALSE, also set lower triangle of covariance matrix
          return(compactToMat(do.call(Covariance, c(list(d=distMat*(1/theta)), cov.args)), diagVal, lower.tri=TRUE))
      }
      else {
        # distMat is a full matrix
        return(do.call(Covariance, c(list(d = distMat/theta), cov.args)))
      }
    }
    # or multiply cross covariance by C
    # as coded below this is not particularly efficient of memory
    else if(!is.na(C[1])) {
      if(onlyUpper) {
        #the onlyUpper option is not compatible with the C option
        onlyUpper = FALSE
      }
      
      if(is.null(x2))
        bigD <- do.call(Distance, c(list(x1, x1), Dist.args))
      else
        bigD <- do.call(Distance, c(list(x1=x1, x2=x2), Dist.args))
      
      if (derivative == 0) {
        return(do.call(Covariance, c(list(d = bigD*(1/theta)), cov.args)) %*% C)
      }
      else {
        # find partial derivatives
        tempW <- do.call(Covariance, c(list(d = bigD*(1/theta)), 
                                       cov.args, derivative = derivative))
        # loop over dimensions and knock out each partial accumulate these in
        # in temp
        temp <- matrix(NA, ncol = d, nrow = n1)
        for (kd in 1:d) {
          # Be careful if the distance (tempD) is close to zero.
          # Note that the x1 and x2 are in transformed ( V inverse) scale
          sM <- ifelse(bigD == 0, 0, tempW * outer(x1[, kd], x2[, kd], "-")/(theta * bigD))
          # accumlate the new partial
          temp[, kd] <- sM %*% C
        }
        # transform back to original coordinates.
        if (!is.null(V)) {
          temp <- temp %*% t(solve(V))
        }
        return(temp)
      }
    }
    # or find marginal variance and return  a vector.
    else if (marginal) {
      sigma2 <- do.call(Covariance, c(list(d = 0), cov.args))
      return(rep(sigma2, nrow(x1)))
    }
    
    # should not get here based on sequence of conditional if statements above.
  }
  
  "cat.matrix" <- function(mat, digits = 8) {
    nc <- ncol(mat)
    temp <- matrix(match(c(signif(mat, digits)), unique(c(signif(mat, 
                                                                 digits)))), ncol = nc)
    temp2 <- format(temp[, 1])
    if (nc > 1) {
      for (k in 2:nc) {
        temp2 <- paste(temp2, temp[, k], sep = "X")
      }
    }
    match(temp2, unique(temp2))
  }
  
  "fields.mkpoly" <- function(x, m = 2) {
    if (m < 0) 
      stop("'m' has to be zero or larger.")
    if( m==0){
      #      warning("Note: There is no polynomial fixed component")
      return( NULL)
    }
    if (!is.matrix(x)) 
      x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
    nterms <- choose((m + d - 1), d)
    temp <- .Fortran("dmaket",PACKAGE="fields", m = as.integer(m), n = as.integer(n), 
                     dim = as.integer(d), des = as.double(x), lddes = as.integer(n), 
                     npoly = as.integer(nterms), tmatrix = as.double(rep(0, 
                                                                         n * (nterms))), ldt = as.integer(n), wptr = as.integer(rep(0, 
                                                                                                                                    d * m)), info = as.integer(0), ptab = as.integer(rep(0, 
                                                                                                                                                                                         nterms * d)), ldptab = as.integer(nterms))
    temp2 <- matrix(temp$tmatrix, nrow = n)
    attr(temp2, "ptab") <- matrix(temp$ptab, nrow = nterms, ncol = d)
    temp2
  }
  
  mKrig.trace <- function(object, iseed, NtrA) {
    set.seed(iseed)
    # if more MonteCarlo samples > number of data points just
    # find A exactly using  np  calls to predict.
    np<- object$np
    if (NtrA >= object$np) {
      Ey <- diag(1, np)
      NtrA <- np
      hold <- diag(predict.mKrig(object, ynew = Ey, collapseFixedEffect=FALSE))
      trA.info<- NA
      trA.est <- sum(hold)
    }
    else {
      # if fewer tests then use random trace method
      # find fitted.values  for iid N(0,1) 'data' to calculate the
      # the Monte Carlo estimate of tr A(lambda)
      # basically repeat the steps above but take some
      # short cuts because we only need fitted.values
      # create random normal 'data'
      Ey <- matrix(rnorm(np * NtrA), nrow = np, 
                   ncol = NtrA)
      trA.info <- colSums(Ey * (predict.mKrig(object, ynew = Ey,
                                              collapseFixedEffect=FALSE)))
      trA.est <- mean(trA.info)
    }
    if (NtrA < np) {
      MSE<-(sum(object$residuals^2)/np) 
      GCV <-       MSE/(1 - trA.est /np)^2
      GCV.info <- MSE/( 1 - trA.info/np)^2
    }
    else{
      GCV<- NA
      GCV.info <- NA
    }	
    return(
      list(trA.info = trA.info, eff.df = trA.est,
           GCV= GCV, GCV.info=GCV.info)
    )
  }
  
  mKrig.coef <- function(object, y, collapseFixedEffect=TRUE) {
    # given new data y and the matrix decompositions in the
    # mKrig object find coefficients d and c.
    # d are the coefficients for the fixed part
    # in this case hard coded for a low order polynomial
    # c are coefficients for the basis functions derived from the
    # covariance function.
    #
    # see mKrig itself for more comments on the linear algebra
    #
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    # generalized least squares for d
    if( any(is.na(y))){
      stop("mKrig can not omit missing values in observation vecotor")
    }
    if( object$nt>0){
      d.coef <- as.matrix(qr.coef(object$qr.VT, forwardsolve(object$Mc, 
                                                             transpose = TRUE, y, upper.tri = TRUE)))
      d.coefMeans<- rowMeans( d.coef)
      if( collapseFixedEffect){ 
        d.coef<- matrix( d.coefMeans, ncol=ncol(d.coef), nrow= nrow( d.coef))
      }
      #  residuals from subtracting off fixed part
      #  of model as m-1 order polynomial
      resid <- y - object$Tmatrix %*% d.coef
    }
    else{
      d.coef<- NULL
      resid <- y
    }
    # and now find c.
    c.coef <- forwardsolve(object$Mc, transpose = TRUE, resid, 
                           upper.tri = TRUE)
    c.coef <- as.matrix(backsolve(object$Mc, c.coef))
    out <- list(d = (d.coef), c = (c.coef))
    return(out)
  }
  
  summary.mKrig <- function(object, ...) {
    
    outObject<- list()  
    digits<- 4
    
    if (is.matrix(object$residuals)) {
      n <- nrow(object$residuals)
      nData <- ncol(object$residuals)
    }
    else {
      n <- length(object$residuals)
      nData <- 1
    }
    outObject$nData<-nData
    
    c1 <- "Number of Locations:"
    c2 <- n
    
    if (nData > 1) {
      c1 <- c(c1, "Number of data sets fit:")
      c2 <- c(c2, nData)
      
    }
    c1 <- c(c1, "Degree of polynomial null space ( base model):")
    if(object$m !=0 ){
      c2 <- c(c2, object$m - 1)
    }
    else{
      c2 <- c(c2, NA)
    }
    c1 <- c(c1, "Total number of parameters in base model")
    c2 <- c(c2, object$nt)
    if (object$nZ > 0) {
      c1 <- c(c1, "Number of additional covariates (Z)")
      c2 <- c(c2, object$nZ)
    }
    if (!is.na(object$eff.df)) {
      c1 <- c(c1, " Eff. degrees of freedom")
      c2 <- c(c2, signif(object$eff.df, digits))
      if (length(object$trA.info) < object$np) {
        c1 <- c(c1, "   Standard Error of estimate: ")
        c2 <- c(c2, signif(sd(object$trA.info)/sqrt(length(object$trA.info)), 
                           digits))
      }
    }
    c1 <- c(c1, "Smoothing parameter")
    c2 <- c(c2, signif(object$lambda.fixed, digits))
    
    if (nData == 1) {
      c1 <- c(c1, "MLE sigma ")
      c2 <- c(c2, signif(object$shat.MLE, digits))
      c1 <- c(c1, "MLE rho")
      c2 <- c(c2, signif(object$rho.MLE, digits))
    }
    
    c1 <- c(c1, "Nonzero entries in covariance")
    c2 <- c(c2, object$nonzero.entries)
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    ########### save table of information and the call
    outObject$summaryTable<- sum
    outObject$call<- object$call
    outObject$collapseFixedEffect<- object$collapseFixedEffect
    
    ########### information for SE for fixed effects
    if( outObject$collapseFixedEffect | (nData==1)){
      outObject$fixedEffectsCov<- object$fixedEffectsCov
      SE<- sqrt(diag(outObject$fixedEffectsCov))
      d.coef<-  object$d[,1]
      print( d.coef)
      pValue<- pnorm(abs(d.coef/SE), lower.tail = FALSE)*2
      outObject$fixedEffectsTable<- cbind( signif(d.coef, digits), 
                                           signif(SE, digits),
                                           signif(pValue, digits)
      )
      if( is.null( object$fixedEffectNames )){
        outObject$fixedEffectNames<- paste0("d",1:(object$nt) )
      }
      else{
        outObject$fixedEffectNames<- object$fixedEffectNames
      }
      dimnames( outObject$fixedEffectsTable) <- list( outObject$fixedEffectNames,
                                                      c("estimate", "SE", "pValue") )
      ########### save covariance information
    } 
    outObject$cov.function<- object$cov.function
    outObject$args<- object$args
    class( outObject)<-"mKrigSummary"
    return( outObject)
  }
  
  print.mKrig<- function (x, digits = 4,...){
    object<- summary.mKrig( x,...)
    print.mKrigSummary(object, digits = digits)
  }
  
  print.mKrigSummary <- function (x, digits = 4, ...) 
  {
    cat("Call:\n")
    dput(x$call)
    print(x$summaryTable, quote = FALSE)
    
    # fixed effects are reported differently when fields are replicated.    
    nData<- x$nData
    cat(" ", fill = TRUE)
    if( nData == 1 | x$collapseFixedEffect ){
      cat(" ", fill = TRUE)
      cat("Summary of fixed effects", fill = TRUE)
      print( x$fixedEffectsTable)
    }
    else {
      cat("Estimated fixed effects found separately for each replicate field", 
          fill = TRUE)
    }
    cat(" ", fill = TRUE)
    cat("Covariance Model:", x$cov.function, fill = TRUE)
    if (x$cov.function == "stationary.cov") {
      cat("   Covariance function:  ", ifelse(is.null(x$args$Covariance), 
                                              "Exponential", x$args$Covariance), fill = TRUE)
    }
    if (!is.null(x$args)) {
      cat("   Non-default covariance arguments and their values ", 
          fill = TRUE)
      nlist <- as.character(names(x$args))
      NL <- length(nlist)
      for (k in 1:NL) {
        cat("   Argument:", nlist[k], " ")
        if (object.size(x$args[[k]]) <= 1024) {
          cat("has the value(s): ", fill = TRUE)
          print(x$args[[k]])
        }
        else {
          cat("too large to print value, size > 1K ...", 
              fill = TRUE)
        }
      }
    }
    invisible(x)
  }
  
  predict.mKrig <- function(object, xnew = NULL, ynew = NULL, grid.list=NULL,
                            derivative = 0, Z = NULL, drop.Z = FALSE, just.fixed = FALSE,
                            collapseFixedEffect = object$collapseFixedEffect, 
                            ...) {
    # the main reason to pass new args to the covariance is to increase
    # the temp space size for sparse multiplications
    # other optional arguments that typically describe the covariance function 
    # from mKrig are passed along in the list object$args
    cov.args <- list(...)
    # predict at observation locations by default
    if( !is.null(grid.list)){
      xnew<- make.surface.grid(grid.list)
    }
    if (is.null(xnew)) {
      xnew <- object$x
    }
    if (is.null(Z) & (length(object$ind.drift) >0 )) {
      Z <- object$Tmatrix[, !object$ind.drift]
    }
    if (!is.null(ynew)) {
      coef.hold <- mKrig.coef(object, ynew,
                              collapseFixedEffect=collapseFixedEffect)
      c.coef <- coef.hold$c
      d.coef <- coef.hold$d
    }
    else {
      c.coef <- object$c
      d.coef <- object$d
    }
    # fixed part of the model this a polynomial of degree m-1
    # Tmatrix <- fields.mkpoly(xnew, m=object$m)
    # only do this if nt>0, i.e. there is a fixed part.
    #
    if( object$nt>0){
      if (derivative == 0) {
        if (drop.Z | object$nZ == 0) {
          # just evaluate polynomial and not the Z covariate
          temp1 <- fields.mkpoly(xnew, m = object$m) %*% 
            d.coef[object$ind.drift, ]
        }
        else {
          if( nrow( xnew) != nrow(as.matrix(Z)) ){
            stop(paste("number of rows of covariate Z",
                       nrow(as.matrix(Z)), 
                       " is not the same as the number of locations",
                       nrow( xnew) )
            )
          }
          temp0 <-  cbind(fields.mkpoly(xnew, m = object$m),as.matrix(Z)) 
          temp1 <- temp0 %*% d.coef
        }
      }
      else {
        if (!drop.Z & object$nZ > 0) {
          stop("derivative not supported with Z covariate included")
        }
        temp1 <- fields.derivative.poly(xnew, m = object$m, d.coef[object$ind.drift, 
        ])
      }
      if (just.fixed) {
        return(temp1)
      }
    }  
    # add nonparametric part. Covariance basis functions
    # times coefficients.
    # syntax is the name of the function and then a list with
    # all the arguments. This allows for different covariance functions
    # that have been passed as their name.
    if (derivative == 0) {
      # argument list are the parameters and other options from mKrig
      #  locations and coefficients,
      temp2 <- do.call(object$cov.function.name, c(object$args, 
                                                   list(x1 = xnew, x2 = object$knots, C = c.coef), cov.args))
    }
    else {
      temp2 <- do.call(object$cov.function.name, c(object$args, 
                                                   list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
                                                   cov.args))
    }
    # add two parts together and coerce to vector
    if( object$nt>0){
      return((temp1 + temp2))
    }
    else{
      return(  temp2)
    }
  }
  
  
  "splint" <- function(x, y, xgrid, wt = NULL, derivative = 0, 
                       lam = 0, df = NA, lambda = NULL, nx=NULL) {
    #
    # reform calling args if passed as a matrix or list
    
    if (is.matrix(x)) {
      if (ncol(x) > 1) {
        xgrid <- y
        y <- x[, 2]
        x <- x[, 1]
      }
    }
    if (is.list(x)) {
      xgrid <- y
      y <- x$y
      x <- x$x
    }
    if (any(duplicated(x))) {
      stop("duplicated x values, use sreg")
    }
    if ((derivative > 2) | (derivative < 0)) 
      stop("derivative must be 0,1,2")
    if (length(x) != length(y)) 
      stop("Lengths of x and y must match")
    n <- length(x)
    #default values for weights
    # NOTE: weights do not matter when interpolating (lam==0)
    if (is.null(wt)) {
      wt <- rep(1, n)
    }
    # find lambda from eff degrees of freedom if it is passed
    if (!is.na(df)) {
      if ((df < 2) | (df > n)) {
        stop("df out of range")
      }
      lam <- sreg.df.to.lambda(df, x, wt)
    }
    # use lambda is it is passed
    if (!is.null(lambda)) {
      lam <- lambda
    }
    igcv <- ifelse(lam == 0, 2, 0)
    # call to FORTRAN -- only return the evaluated poiints (ygrid).
    if( !is.null(nx)){
      xgrid<- seq( min( x), max(x),,nx)
    }
    
    ygrid<- .Fortran("css",PACKAGE="fields",
                     h = as.double(ifelse(igcv == 2, 1, log(lam))),
                     as.integer(n),
                     as.double(x),
                     as.double(y), 
                     wt = as.double(1/sqrt(wt)),
                     sy = as.double(rep(0, n)), 
                     as.double(1),
                     as.double(1),
                     as.double(1),
                     as.integer(length(xgrid)), 
                     as.double(xgrid),
                     ygrid = as.double(rep(0, length(xgrid))), 
                     job = as.integer(c(igcv, 3, 0)),
                     as.integer(derivative), 
                     as.integer(0)
    )$ygrid
    if(!is.null(nx) ){ 
      return(list( x=xgrid, y=ygrid))
    }
    else{ 
      return( ygrid)
    }
  }
  
  data <- data %>% 
    st_transform(st_crs(crs))
  
  data$Y <- as.numeric(map_chr(data$geometry, 2))
  data$X <- as.numeric(map_chr(data$geometry, 1))
  
  crs <- crs 
  
  summary_of_df <- cbind.data.frame('out' =  out, 'treat' = treat)
  mean_dv <- mean(summary_of_df$out[summary_of_df$treat==0])
  sd_dv <- sd(out)
  
  table_data <- function(ols,Vcov){
    X = qr.X(ols$qr)
    N = nrow(X)
    V = t(X) %*% (Vcov) %*% X
    r.inv = solve(qr.R(ols$qr))
    xx.inv = r.inv %*% t(r.inv)
    Cov = xx.inv %*% V %*% xx.inv
    se = sqrt(diag(Cov))
    t.stat = summary(ols)$coef[, 1] / se
    p.value = 2 * (1 - pt(abs(t.stat), df = ols$df.residual))
    coef = data.frame(
      Estimate = summary(ols)$coef[, 1],
      'Std Error' = se,
      't value' = t.stat,
      'Pr t' = p.value,
      'N' = size_of_df, 
      'bw' = round(bw,2), 
      'mean_dv' = round(mean_dv, 2), 
      'sd_dv' = round(sd_dv, 2), 
      'Range' = round(Range,2), 
      'Structure' = round(Structure,2) 
    )
    return(coef)
  }
  
  if(user_bw == FALSE){
    
    
    if(lat_long == TRUE){
      
      if(adjust == TRUE){
        
        bw = rdbwselect(out, run_var, covs = covs, kernel = 'uni', p = p)$bws[1]
        size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
        
        ols <- lm(out ~ treat + poly(X, p)*poly(Y, p) + covs,  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
        
        Residuals <- ols$residuals
        data$Residuals <- ols$residuals
        
        data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
        data$geometry <- NULL 
        coords = data %>% dplyr::select(X,Y)
        data=SpatialPointsDataFrame(coords=coords,data=data)
        proj4string(data)=CRS(crs)
        
        kappa=0.5
        
        fit = spatialProcess(
          x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                    estimateAnisotropy(data,Residuals,Residuals~1))@coords,
          y = data$Residuals,
          Distance = "rdist.earth",
          cov.args = list(Covariance = "Matern", smoothness = kappa),
          mKrig.args = list(m = 1),
          theta.start = theta_start
          
        )
        cov.params=as.vector(fit$summary)
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        VCov = cov.params[5] * 
          Matern(
            rdist.earth(x1 = fit$x),
            range = cov.params[3],
            smoothness = kappa) +
          diag(fit$N) * cov.params[4] ^ 2
        
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        
        table_data(ols, VCov)
      }
      
      else{
        bw = rdbwselect(out, run_var, kernel = 'uni', p = p)$bws[1]
        size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
        
        ols <- lm(out ~ treat + poly(X, p)*poly(Y, p),  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
        
        Residuals <- ols$residuals
        data$Residuals <- ols$residuals
        
        data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
        data$geometry <- NULL 
        coords = data %>% dplyr::select(X,Y)
        data=SpatialPointsDataFrame(coords=coords,data=data)
        proj4string(data)=CRS(crs)
        
        kappa=0.5
        
        fit = spatialProcess(
          x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                    estimateAnisotropy(data,Residuals,Residuals~1))@coords,
          y = data$Residuals,
          Distance = "rdist.earth",
          cov.args = list(Covariance = "Matern", smoothness = kappa),
          mKrig.args = list(m = 1),
          theta.start = theta_start
          
        )
        cov.params=as.vector(fit$summary)
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        VCov = cov.params[5] * 
          Matern(
            rdist.earth(x1 = fit$x),
            range = cov.params[3],
            smoothness = kappa) +
          diag(fit$N) * cov.params[4] ^ 2
        
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        
        table_data(ols, VCov)
      }
    }
    
    else{ 
      if(adjust == TRUE){
        if(p == 1){
          
          bw = rdbwselect(out, run_var, covs = covs, 
                          kernel = 'uni', p = p)$bws[1]
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var + covs,  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }
        else{
          
          bw = rdbwselect(out, run_var, covs = covs, kernel = 'uni', p = p)$bws[1]
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var + covs + treat*I(run_var^p),  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }}
      
      else{
        
        if(p == 1){
          bw = rdbwselect(out, run_var, kernel = 'uni', p = p)$bws[1]
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var,  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }
        else{
          bw = rdbwselect(out, run_var, kernel = 'uni', p = p)$bws[1]  
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var + treat*I(run_var^p),  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }
      }
    }
  }
  else{
    if(lat_long == TRUE){
      
      if(adjust == TRUE){
        
        bw = input_bw
        size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
        
        ols <- lm(out ~ treat + poly(X, p)*poly(Y, p) + covs,  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
        
        Residuals <- ols$residuals
        data$Residuals <- ols$residuals
        
        data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
        data$geometry <- NULL 
        coords = data %>% dplyr::select(X,Y)
        data=SpatialPointsDataFrame(coords=coords,data=data)
        proj4string(data)=CRS(crs)
        
        kappa=0.5
        
        fit = spatialProcess(
          x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                    estimateAnisotropy(data,Residuals,Residuals~1))@coords,
          y = data$Residuals,
          Distance = "rdist.earth",
          cov.args = list(Covariance = "Matern", smoothness = kappa),
          mKrig.args = list(m = 1),
          theta.start = theta_start
          
        )
        cov.params=as.vector(fit$summary)
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        VCov = cov.params[5] * 
          Matern(
            rdist.earth(x1 = fit$x),
            range = cov.params[3],
            smoothness = kappa) +
          diag(fit$N) * cov.params[4] ^ 2
        
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        
        table_data(ols, VCov)
      }
      
      else{
        bw = input_bw
        size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
        
        ols <- lm(out ~ treat + poly(X, p)*poly(Y, p),  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
        
        Residuals <- ols$residuals
        data$Residuals <- ols$residuals
        
        data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
        data$geometry <- NULL 
        coords = data %>% dplyr::select(X,Y)
        data=SpatialPointsDataFrame(coords=coords,data=data)
        proj4string(data)=CRS(crs)
        
        kappa=0.5
        
        fit = spatialProcess(
          x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                    estimateAnisotropy(data,Residuals,Residuals~1))@coords,
          y = data$Residuals,
          Distance = "rdist.earth",
          cov.args = list(Covariance = "Matern", smoothness = kappa),
          mKrig.args = list(m = 1),
          theta.start = theta_start
          
        )
        cov.params=as.vector(fit$summary)
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        VCov = cov.params[5] * 
          Matern(
            rdist.earth(x1 = fit$x),
            range = cov.params[3],
            smoothness = kappa) +
          diag(fit$N) * cov.params[4] ^ 2
        
        Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
        Range=as.vector(sqrt(8*kappa)*fit$summary[3])
        
        table_data(ols, VCov)
      }
    }
    
    else{ 
      if(adjust == TRUE){
        if(p == 1){
          
          bw = input_bw
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var + covs,  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }
        else{
          
          bw = input_bw
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var + covs + treat*I(run_var^p),  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }}
      
      else{
        
        if(p == 1){
          bw = input_bw
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var,  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }
        else{
          bw = input_bw
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          ols <- lm(out ~ treat*run_var + treat*I(run_var^p),  data = data, weights = ifelse(abs(run_var) <= bw, 1,0))
          
          size_of_df <- dim(subset(data, abs(run_var) <= bw))[1]
          
          Residuals <- ols$residuals
          data$Residuals <- ols$residuals
          
          data <- as.data.frame(subset(data, abs(run_var) <= bw)) 
          data$geometry <- NULL 
          coords = data %>% dplyr::select(X,Y)
          data=SpatialPointsDataFrame(coords=coords,data=data)
          proj4string(data)=CRS(crs)
          
          kappa=0.5
          
          fit = spatialProcess(
            x = rotateAnisotropicData(SpatialPoints(coords=coords),
                                      estimateAnisotropy(data,Residuals,Residuals~1))@coords,
            y = data$Residuals,
            Distance = "rdist.earth",
            cov.args = list(Covariance = "Matern", smoothness = kappa),
            mKrig.args = list(m = 1),
            theta.start = theta_start
            
          )
          cov.params=as.vector(fit$summary)
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Effective.Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          VCov = cov.params[5] * 
            Matern(
              rdist.earth(x1 = fit$x),
              range = cov.params[3],
              smoothness = kappa) +
            diag(fit$N) * cov.params[4] ^ 2
          
          Structure=cov.params[5]/(cov.params[5]+cov.params[4]^2)
          Range=as.vector(sqrt(8*kappa)*fit$summary[3])
          
          table_data(ols, VCov)
        }
      }
    }
  }
}
