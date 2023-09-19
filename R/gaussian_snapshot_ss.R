# gaussian_snapshot_ss: large visible function

# helper: generate nonparametric (scaled trig) processes
# given dimensions and snapshot times
process_nonpar <- function(n,d,x_vec,
                           amplitude,
                           frequency,
                           int_sd){
  # dimensions
  m <- length(x_vec)
  x_min <- min(x_vec)
  x_max <- max(x_vec)
  # scaled trig functions:
  # fixed parameters
  # minimum amplitude, intercept sd
  amp_min <- .5
  # populate snapshots
  out <- array(NA,c(n,d,m))
  for(ii in 1:n){
    for(rr in 1:d){
      # intercept
      G <- stats::rnorm(1,sd=int_sd[rr])
      # increase/decrease amplitude
      B <- stats::rbinom(1,1,.5)
      # sin x-intercepts
      U <- stats::runif(1,min=x_min,max=x_max)
      # function
      temp <- function(x){
        ((amplitude*(sin(6.28*(frequency/(x_max-x_min))*(x - U))))/(1 + ((amplitude/amp_min) - 1)*(x + B*(x_max-2*x)))) + G
      }
      # populate
      out[ii,rr,] <- temp(x_vec)
    }
  }
  return(out)
}

#' Simulate Gaussian edge networks with nonparametric latent processes
#'
#' \code{gaussian_snapshot_ss} simulates a realization of a functional network
#' with Gaussian edges, according to an inner product latent process model.
#' The latent processes are randomly generated sinusoidal functions.
#'
#' The the latent process for node \eqn{i} in latent dimension \eqn{r} is given independently by
#' \deqn{z_{i,r}(x) = \frac{a \sin [2f\pi(x - U) / (x_{max} - x_{min})]}{1 + (2a-1)[x + B(x_{max} - 2x)]} + G}
#' Where \eqn{G} is Gaussian with mean \code{0} and standard deviation
#' \eqn{\sigma_{int,r}}, \eqn{B} is Bernoulli with mean \code{1/2}, and \eqn{U} is uniform
#' with minimum \code{spline_design$x_min} and maximum \code{spline_design$x_max}.
#' \eqn{f} is a frequency parameter specified with
#' \code{process_options$frequency}, and \eqn{a} is a maximum amplitude parameter
#' specified with \code{process_options$amplitude}.
#' Roughly, each process is a randomly shifted sine function which goes through
#' \code{f} cycles on the index set, with amplitude either increasing or
#' decreasing between \eqn{1/2} and \eqn{a}.
#'
#' Then, the \eqn{n \times n} symmetric adjacency matrix for
#' snapshot \eqn{k=1,...,m} has independent Gaussian entries
#' with standard deviation \code{sigma_edge} and mean
#' \deqn{E([A_k]_{ij}) = z_i(x_k)^{T}z_j(x_k)}
#' for \eqn{i \leq j} (or \eqn{i < j} with no self loops).
#'
#' This function may return the latent processes as an \eqn{n \times d \times m}
#' array evaluated at the prespecified snapshot indices, or as a function which
#' takes a vector of indices and returns the corresponding evaluations of
#' the latent process matrices.
#' It also returns the spline design information required to
#' fit a FASE embedding to this data with a natural cubic spline.
#'
#' @usage
#' gaussian_snapshot_ss(n,d,m,x_vec,self_loops=TRUE,
#'                      sigma_edge=1,process_options)
#'
#' @param n A positive integer, the number of nodes.
#' @param d A positive integer, the number of latent space dimensions.
#' @param m A positive integer, the number of snapshots.
#' If this argument is not specified, it
#' is determined from the snapshot index vector \code{x_vec}.
#' @param x_vec A vector, the snapshot evaluation indices for the data.
#' Defaults to an equally spaced sequence of length
#' \code{m} from \code{0} to \code{1}.
#' @param self_loops A Boolean, if \code{FALSE}, all diagonal adjacency matrix
#' entries are set to zero. Defaults to \code{TRUE}.
#' @param sigma_edge A positive scalar,
#' the entry-wise standard deviation for the Gaussian edge variables.
#' Defaults to \code{1}.
#' @param process_options A list, containing additional optional arguments:
#' \describe{
#'     \item{amplitude}{A positive scalar, the maximum amplitude of the
#'     randomly generated latent processes. Defaults to \code{3}.}
#'     \item{frequency}{A positive scalar, frequency of the randomly
#'     generated latent processes. Defaults to \code{2}.}
#'     \item{sigma_int}{A positive scalar, or a vector of length \eqn{d}.
#'     If it is a vector, the entries correspond to the standard deviation
#'     of the random intercepts of the node processes for each latent dimension.
#'     If is is a scalar, it corresponds to the standard deviation of the
#'     random intercepts in all dimensions. Defaults to \code{0.5}.}
#'     \item{return_fn}{A Boolean, if \code{TRUE}, then the latent processes
#'     are returned as a function which
#'     takes a vector of indices and returns the corresponding evaluations of
#'     the latent process matrices. Otherwise, the latent processes are returned
#'     as an \eqn{n \times d \times m} array evaluated at the prespecified
#'     snapshot indices. Defaults to \code{FALSE}.}
#' }
#'
#'
#'
#' @return A list is returned with the realizations of the basis coordinates,
#' spline design, and the multiplex network snapshots:
#' \item{A}{An array of dimension \eqn{n \times n \times m}, the realized
#' functional network data.}
#' \item{Z}{If \code{process_options$return_fn} is \code{TRUE}, a function,
#' which takes a vector of indices and returns the corresponding evaluations of
#' the latent process matrices. Otherwise, an array of dimension \eqn{n \times d \times m},
#' the latent processes evaluated at the prespecified
#' snapshot indices.}
#' \item{spline_design}{A list, describing the \eqn{B}-spline design:
#' \describe{
#'     \item{type}{The string \code{'ss'}.}
#'     \item{x_vec}{A vector, the snapshot evaluation indices for the data.}
#' }}
#'
#' @examples
#'
#' # Gaussian edge data with sinusoidal latent processes
#' # NOTE: latent processes are returned as a function
#'
#' data <- gaussian_snapshot_ss(n=100,d=2,
#'                              x_vec=seq(0,3,length.out=80),
#'                              self_loops=TRUE,
#'                              sigma_edge=4,
#'                              process_options=list(amplitude=4,
#'                                                   frequency=3,
#'                                                   return_fn=TRUE))
#'
#' @export
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
  if(is.null(process_options$sigma_int)){
    process_options$sigma_int <- rep(0.5,d)
  }
  else{
    if(length(process_options$sigma_int)==1){
      process_options$sigma_int <- rep(process_options$sigma_int,d)
    }
    else{
      if(length(process_options$sigma_int) != d){
        stop('invalid dimension of intercept standard deviations, must be 1 or d')
      }
    }
  }
  if(is.null(process_options$return_fn)){
    process_options$return_fn <- FALSE
  }
  # Z populate
  Z_array <- process_nonpar(n,d,x_vec,
                            process_options$amplitude,
                            process_options$frequency,
                            process_options$sigma_int)
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
  spline_design <- list(type='ss',x_vec=x_vec,x_min=min(x_vec),x_max=max(x_vec))
  # if specified, interpolate to get functional output
  if(process_options$return_fn){
    Z <- coord_to_func(aperm(Z_array,c(1,3,2)),spline_design)
  }
  else{
    Z <- Z_array
  }
  out <- list(A=A,Z=Z,spline_design=spline_design)
  return(out)
}

