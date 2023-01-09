# rdpg_snapshot_bs: large visible function

# helper: rdpg Dirichlet spline coordinates
coord_rdpg <- function(n,q,d,
                       B_mat,
                       alpha,density){
  W <- array(0,c(n,q,d))
  for(ii in 1:n){
    temp <- MCMCpack::rdirichlet(q,rep(alpha,d))
    normalizer <- max(sqrt(rowSums((B_mat %*% temp)^2)))
    if(normalizer > 1){
      tempnorm <- temp/normalizer
    }
    else{
      tempnorm <- temp
    }
    W[ii,,] <- sqrt(d*density)*tempnorm
  }
  return(W)
}


rdpg_snapshot_bs <- function(n,d,
                             m=NULL,
                             self_loops=TRUE,
                             spline_design=list(),
                             process_options=list()){
  # error checking
  if(is.null(m) & is.null(spline_design$x_vec)){
    stop('must specify either m or an index vector')
  }
  if(is.null(spline_design$q)){
    stop('must specify B-spline dimension')
  }

  # spline design parameter checking
  if(is.null(m)){
    m <- length(spline_design$x_vec)
  }
  if(is.null(spline_design$x_vec)){
    spline_design$x_vec <- seq(0,1,length.out=m)
  }
  if(is.null(spline_design$x_min)){
    spline_design$x_min <- min(spline_design$x_vec)
  }
  if(is.null(spline_design$x_max)){
    spline_design$x_max <- max(spline_design$x_vec)
  }
  # consistency of parameters
  if(spline_design$x_min > min(spline_design$x_vec) || spline_design$x_max > max(spline_design$x_vec)){
    stop('inconsistent index space')
  }
  if(spline_design$q > m || spline_design$q < 4){
    stop('invalid choice of q')
  }
  # spline design populate
  # type
  spline_design$type <- 'bs'
  # B-spline design
  spline_design$spline_mat <- B_func(spline_design$q,
                                     spline_design$x_min,
                                     spline_design$x_max)(spline_design$x_vec)
  # process options checking
  if(is.null(process_options$alpha_coord)){
    process_options$alpha_coord <- .1
  }
  if(is.null(process_options$density)){
    process_options$density <- 1/d
  }
  # parameter checking
  if(process_options$density > (1/d)){
    cat(paste0('requested density is too high, networks will have edge density approximately ',round(1/d,2)))
    process_options$density <- 1/d
  }
  # populate Theta, A
  W <- coord_rdpg(n,spline_design$q,d,
                  spline_design$spline_mat,
                  process_options$alpha_coord,
                  process_options$density)
  # start with mean structure
  Theta <- WB_to_Theta(W,spline_design$spline_mat,self_loops)
  A <- array(0,dim(Theta))
  for(kk in 1:m){
    # truncation step for P
    Theta.tri <- pmin(pmax(c(Theta[,,kk][lower.tri(Theta[,,kk],diag=T)]),0),1)
    A.temp <- matrix(NA,n,n)
    A.temp[lower.tri(A.temp,diag=T)] <- as.numeric(stats::rbinom(n=length(Theta.tri),size=1,prob=Theta.tri))
    A.temp[upper.tri(A.temp)] <- t(A.temp)[upper.tri(A.temp)]
    A[,,kk] <- A.temp
  }
  # populate output
  out <- list(A=A,W=W,spline_design=spline_design)
  return(out)
}

