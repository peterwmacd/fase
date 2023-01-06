# fase: main visible function

fase <- function(A,d,self_loops=TRUE,
                 spline_design=list(),
                 lambda=0,
                 optim_options=list(),
                 output_options=list()){
  # dimensions
  n <- dim(A)[1]
  m <- dim(A)[3]
  # check spline_design contents and construct objects
  if(is.null(spline_design$type)){
    spline_design$type <- 'bs'
  }
  # for B-splines
  if(spline_design$type=='bs'){
    # spline parameter checking
    if(is.null(spline_design$q)){
      stop('must specify B-spline dimension')
    }
    if(is.null(spline_design$x_min)){
      spline_design$x_min <- 0
    }
    if(is.null(spline_design$x_max)){
      spline_design$x_max <- 1
    }
    if(is.null(spline_design$x_vec)){
      spline_design$x_vec <- seq(spline_design$x_min,
                                 spline_design$x_max,
                                 length.out=m)
    }
    # consistency of parameters
    if(spline_design$x_min > min(spline_design$x_vec) || spline_design$x_max > max(spline_design$x_vec)){
      stop('inconsistent index space')
    }
    if(spline_design$q > m || spline_design$q < 4){
      stop('invalid choice of q')
    }
    # populate
    # B-spline design
    spline_design$spline_mat <- B_func(spline_design$q,
                                       spline_design$x_min,
                                       spline_design$x_max)(spline_design$x_vec)
    # ridge matrix
    if(lambda > 0){
      spline_design$ridge_mat <- diag(m)
    }
  }
  else{
    if(spline_design$type=='ss'){
      # spline paramter checking
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
      if(spline_design$x_min != min(spline_design$x_vec) || spline_design$x_max != max(spline_design$x_vec)){
        stop('inconsistent index space')
      }
      # populating
      # S-spline design
      spline_design$spline_mat <- splines2::naturalSpline(spline_design$x_vec,
                                                          intercept=TRUE,
                                                          knots=spline_design$x_vec[2:(m-1)],
                                                          Boundary.knots=c(spline_design$x_vec[1],spline_design$x_vec[m]))
      # ridge matrix
      if(lambda > 0){
        spline_design$ridge_mat <- ss_ridge(spline_design)
      }
    }
    else{
      stop('FASE supports B-splines and S-splines')
    }
  }

  # optimization parameter checking
  if(is.null(optim_options$eta)){
    optim_options$eta <- 1/(n*m)
  }
  if(is.null(optim_options$eps)){
    optim_options$eps <- 1e-5
  }
  if(is.null(optim_options$K_max)){
    optim_options$K_max <- 2e3
  }
  if(is.null(optim_options$verbose)){
    optim_options$verbose <- FALSE
  }
  if(is.null(optim_options$init_band)){
    if(spline_design$type=='bs'){
      optim_options$init_band <- spline_design$q/(m*(spline_design$x_max - spline_design$x_min))
    }
    else{
      optim_options$init_band <- 10/(m*(spline_design$x_max - spline_design$x_min))
    }
  }
  if(is.null(optim_options$init_q)){
    if(spline_design$type=='bs'){
      optim_options$init_q <- spline_design$q
    }
    else{
      optim_options$init_q <- 10
    }
  }
  # consistency of parameters
  if(optim_options$init_q > m){
    stop('invalid choice of initial q (too large)')
  }
  if((spline_design$type=='bs' & optim_options$init_q < spline_design$q) || optim_options$init_q < 4){
    stop('invalid choice of initial q (too small)')
  }
  # initialization
  W_init <- kern_orth_embed(A,d,
                            spline_design,
                            init_q=optim_options$init_q,
                            band=optim_options$init_band)
  # gradient descent
  if(lambda==0){
    gd_out <- gd_fit(A,self_loops,
                     W_init,
                     spline_design$spline_mat,
                     optim_options$eta,
                     optim_options$eps,
                     optim_options$K_max,
                     verbose=optim_options$verbose)
  }
  else{
    gd_out <- gd_fit(A,self_loops,
                     W_init,
                     spline_design$spline_mat,
                     optim_options$eta,
                     optim_options$eps,
                     optim_options$K_max,
                     ridge=TRUE,
                     ridge_mat=spline_design$ridge_mat,
                     ridge_pen=lambda,
                     verbose=optim_options$verbose)
  }

  # output parameter checking
  if(is.null(output_options$align_output)){
    output_options$align_output <- TRUE
  }
  if(is.null(output_options$return_fn)){
    output_options$return_fn <- FALSE
  }
  if(is.null(output_options$return_coords)){
    output_options$return_coords <- FALSE
  }
  if(spline_design$type=='ss'){
    output_options$return_ngcv <- FALSE
  }
  if(spline_design$type=='bs'){
    if(is.null(output_options$return_ngcv)){
      output_options$return_ngcv <- TRUE
    }
  }
  # populate output
  out <- list()
  # estimated processes
  # process snapshots
  if(!output_options$return_fn){
    # align first
    if(output_options$align_output){
      out$Z <- Z_align_proc(coord_to_snap(gd_out$W_hat,spline_design$spline_mat))
    }
    # use coordinates directly
    else{
      out$Z <- coord_to_snap(gd_out$W_hat,spline_design$spline_mat)
    }
  }
  # process function
  else{
    # align first
    if(output_options$align_output){
      temp_Z <- Z_align_proc(coord_to_snap(gd_out$W_hat,spline_design$spline_mat))
      temp_spline_design <- list(type='ss',
                                 x_vec=spline_design$x_vec)
      out$Z <- coord_to_func(aperm(temp_Z,c(1,3,2)),temp_spline_design)
    }
    # use coordinates directly
    else{
      out$Z <- coord_to_func(gd_out$W_hat,spline_design)
    }
  }
  # coordinates
  if(output_options$return_coords){
    out$W <- gd_out$W_hat
  }
  else{
    out$W <- NA
  }
  # spline design
  out$spline_design <- spline_design
  # NGCV criterion
  if(output_options$return_ngcv){
    n_node <- ifelse(self_loops,n*m,(n-1)*m)
    out$ngcv <- log(sqrt(gd_out$obj/(n*n_node))) - log(1 - ((2*spline_design$q*d)/n_node))
  }
  else{
    out$ngcv <- NA
  }
  # convergence information
  out$K <- gd_out$K
  out$converged <- gd_out$converged
  # return
  return(out)
}
