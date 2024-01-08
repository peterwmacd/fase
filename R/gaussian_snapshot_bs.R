#' Simulate Gaussian edge networks with B-spline latent processes
#'
#' \code{gaussian_snapshot_bs} simulates a realization of a functional network
#' with Gaussian edges, according to an inner product latent process model.
#' The latent processes are generated from a \eqn{B}-spline basis with equally
#' spaced knots.
#'
#' The spline design of the functional network data (snapshot indices,
#' basis dimension) is generated using the information provided in
#' \code{spline_design}, producing a \eqn{q}-dimensional cubic
#' \eqn{B}-spline basis with equally spaced knots.
#'
#' The latent process basis coordinates are generated as iid
#' Gaussian random variables with standard deviation
#' \code{process_options$sigma_coord}. Each latent process is given by
#' \deqn{z_{i,r}(x) = w_{i,r}^{T}B(x).}
#'
#' Then, the \eqn{n \times n} symmetric adjacency matrix for
#' snapshot \eqn{k=1,...,m} has independent Gaussian entries
#' with standard deviation \code{sigma_edge} and mean
#' \deqn{E([A_k]_{ij}) = z_i(x_k)^{T}z_j(x_k)}
#' for \eqn{i \leq j} (or \eqn{i < j} with no self loops).
#'
#' @usage
#' gaussian_snapshot_bs(n,d,m,self_loops=TRUE,
#'                      spline_design,sigma_edge=1,
#'                      process_options)
#'
#' @param n A positive integer, the number of nodes.
#' @param d A positive integer, the number of latent space dimensions.
#' @param m A positive integer, the number of snapshots.
#' If this argument is not specified, it
#' is determined from the snapshot index vector \code{spline_design$x_vec}.
#' @param self_loops A Boolean, if \code{FALSE}, all diagonal adjacency matrix
#' entries are set to zero. Defaults to \code{TRUE}.
#' @param spline_design A list, describing the \eqn{B}-spline design:
#' \describe{
#'     \item{q}{A positive integer, the dimension of the \eqn{B}-spline basis.
#'     Must be at least \code{4} and at most \code{m}.}
#'     \item{x_vec}{A vector, the snapshot evaluation indices for the data.
#'     Defaults to an equally spaced sequence of length
#'     \code{m} from \code{0} to \code{1}.}
#'     \item{x_max}{A scalar, the maximum of the index space.
#'     Defaults to \code{max(spline_design$x_vec)}.}
#'     \item{x_min}{A scalar, the minimum of the index space.
#'     Defaults to \code{min(spline_design$x_vec)}.}
#' }
#' @param sigma_edge A positive scalar,
#' the entry-wise standard deviation for the Gaussian edge variables.
#' Defaults to \code{1}.
#' @param process_options A list, containing additional optional arguments:
#' \describe{
#'     \item{sigma_coord}{A positive scalar, or a vector of length \eqn{d}.
#'     If it is a vector, the entries correspond to the standard deviation
#'     of the randomly generated basis coordinates for each latent dimension.
#'     If is is a scalar, it corresponds to the standard deviation of the
#'     basis coordinates in all dimensions. Defaults to \code{1}.}
#' }
#'
#'
#'
#' @return A list is returned with the realizations of the basis coordinates,
#' spline design, and the multiplex network snapshots:
#' \item{A}{An array of dimension \eqn{n \times n \times m}, the realized
#' functional network data.}
#' \item{W}{An array of dimension \eqn{n \times q \times d},
#' the realized basis coordinates.}
#' \item{spline_design}{A list, describing the \eqn{B}-spline design:
#' \describe{
#'     \item{type}{The string \code{'bs'}.}
#'     \item{q}{A positive integer, the dimension of the \eqn{B}-spline basis.}
#'     \item{x_vec}{A vector, the snapshot evaluation indices for the data.}
#'     \item{x_max}{A scalar, the maximum of the index space.}
#'     \item{x_min}{A scalar, the minimum of the index space.}
#'     \item{spline_matrix}{An \eqn{m \times q} matrix, the B-spline basis
#'     evaluated at the snapshot indices.}
#' }}
#'
#' @examples
#'
#' # Gaussian edge data with B-spline latent processes, Gaussian coordinates
#' # NOTE: x_vec is automatically populated given m
#'
#' data <- gaussian_snapshot_bs(n=100,d=4,m=100,
#'                              self_loops=FALSE,
#'                              spline_design=list(q=12),
#'                              sigma_edge=3,
#'                              process_options=list(sigma_coord=.75))
#'
#' @export
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
  if(spline_design$x_min > min(spline_design$x_vec) || spline_design$x_max < max(spline_design$x_vec)){
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
                                     spline_design$x_max,
                                     spline_design$x_vec)(spline_design$x_vec)
  # process options checking
  if(is.null(process_options$sigma_coord)){
    process_options$sigma_coord <- rep(1,d)
  }
  else{
    if(length(process_options$sigma_coord)==1){
      process_options$sigma_coord <- rep(process_options$sigma_coord,d)
    }
    else{
      if(length(process_options$sigma_coord) != d){
        stop('invalid dimension of coordinate standard deviations, must be 1 or d')
      }
    }
  }
  # populate Theta, A
  W <- aperm(array(stats::rnorm(d*n*spline_design$q,sd=process_options$sigma_coord),c(d,n,spline_design$q)),c(2,3,1))
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

