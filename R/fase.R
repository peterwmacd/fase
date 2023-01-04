# fase main visible function

fase <- function(A,d,self_loops=TRUE,
                 spline_design=list(type='bs',
                                    x_vec=NULL,
                                    x_min=0,
                                    x_max=1,
                                    q=(dim(A)[3]/10)),
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
    if(is.null(spline_design$q)){
      stop('must specify B-spline dimension')
    }
    if(is.null(spline_design$x_min)){
      spline_design$x_min <- 0
    }
    # etc checking

    # add B matrix
    spline_design$B_mat <- B_func(spline_design$q,
                                  spline_design$x_min,
                                  spline_design$x_max)(spline_design$x_vec)
    # add ridge matrix
    if(lambda > 0){
      spline_design$ridge_mat <- diag(m)
    }
  }
  else{
    if(spline_design$type=='ss'){
      # checking

      # populating
    }
    else{
      stop('FASE supports B-splines and S-splines')
    }
  }

  # optimization parameter checking

  # initialization

  # gradient descent

  # output parameter checking

  out <- list()
  # estimated processes


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
    ngcv_out <- log(sqrt(gd_out$obj/prod(dA))) - log(1 - ((2*q_curr*d_curr)/(n*m)))
  }
  else{
    ngcv_out <- NA
  }
  # convergence information
  out$K <- gd_out$K
  out$convergence <- gd_out$convergence
  # return
  return(out)
}
