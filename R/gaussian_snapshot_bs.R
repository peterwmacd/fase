# gaussian_snapshot_bs: large visible function

gaussian_snapshot_bs <- function(n,d,
                                 m=NULL,
                                 self_loops=TRUE,
                                 spline_design=list(),
                                 sigma_edge=1,
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
  if(is.null(process_options$sigma_coord)){
    process_options$sigma_coord <- 1
  }
  # populate Theta, A
  W <- array(stats::rnorm(n*spline_design$q*d,sd=process_options$sigma_coord),c(n,spline_design$q,d))
  # start with mean structure
  A <- WB_to_Theta(W,spline_design$spline_mat,self_loops)
  for(kk in 1:m){
    # error matrix
    E <- matrix(stats::rnorm(n^2,sd=sigma_edge),n,n)*(1 - as.integer(!self_loops)*diag(n))
    E <- upper.tri(E)*E
    E <- E + t(E)
    if(self_loops){
      diag(E) <- stats::rnorm(n,sd=sigma_edge)
    }
    A[,,kk] <- A[,,kk]+E
  }

  # populate output
  out <- list(A=A,W=W,spline_design=spline_design)
  return(out)
}

