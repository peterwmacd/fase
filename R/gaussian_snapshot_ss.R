# gaussian_snapshot_ss: large visible function

# helper: generate nonparametric (scaled trig) processes
# given dimensions and snapshot times
process_nonpar <- function(n,d,x_vec,
                           amplitude,
                           frequency){
  # dimensions
  m <- length(x_vec)
  # scaled trig functions:
  # fixed parameters
  # minimum amplitude, intercept sd
  amp_min <- int_sd <- .5
  # populate snapshots
  out <- array(NA,c(n,d,m))
  for(ii in 1:n){
    for(rr in 1:d){
      # intercept
      G <- rnorm(1,sd=int_sd)
      # increase/decrease amplitude
      B <- rbinom(1,1,.5)
      # sin x-intercepts
      U <- runif(1)
      # function
      temp <- function(x){
        ((amplitude*(sin(6.28*(frequency*x - U))))/(1 + ((amplitude/amp_min) - 1)*(x + B*(1-2*x)))) + G
      }
      # populate
      out[ii,rr,] <- temp(x_vec)
    }
  }
  return(out)
}

gaussian_snapshot_ss <- function(n,d,
                                 m=NULL,x_vec=NULL,
                                 self_loops=TRUE,
                                 sigma_edge=1,
                                 process_options=list()){
  # error checking
  if(is.null(m) & is.null(x_vec)){
    stop('must specify either m or an index vector')
  }
  # spline design parameter checking
  if(is.null(m)){
    m <- length(x_vec)
  }
  if(is.null(x_vec)){
    x_vec <- seq(0,1,length.out=m)
  }
  # process options checking
  if(is.null(process_options$amplitude)){
    process_options$amplitude <- 3
  }
  if(is.null(process_options$frequency)){
    process_options$frequency <- 2
  }
  if(is.null(process_options$return_fn)){
    process_options$return_fn <- FALSE
  }
  # Z populate
  Z_array <- process_nonpar(n,d,x_vec,
                            process_options$amplitude,
                            process_options$frequency)
  # start with mean structure
  A <- Z_to_Theta(Z_array,self_loops)
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
  spline_design <- list(type='ss',x_vec=x_vec)
  # if specified, interpolate to get functional output
  if(process_options$return_fn){
    Z <- coord_to_func(aperm(Z_array,c(1,3,2)),spline_design)
  }
  else{
    Z <- Z_array
  }
  out <- list(A=A,Z=Z,spline_design=list(type='ss',x_vec=x_vec))
  return(out)
}

