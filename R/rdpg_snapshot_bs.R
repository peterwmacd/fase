# helper: rdpg Dirichlet spline coordinates
coord_rdpg <- function(n,q,d,
                       B_mat,
                       alpha,density){
  W <- array(0,c(n,q,d))
  for(ii in 1:n){
    # generate Dirichlet variables
    # matrix of independent chi-square ( == gamma(alpha,2) ) rvs
    vv <- t(matrix(stats::rchisq(d*q,2*alpha),nrow=d))
    # rescale rows to the d-dimensional simplex
    temp <- vv / rowSums(vv)
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

#' Simulate binary edge networks with B-spline latent processes
#'
#' \code{rdpg_snapshot_bs} simulates a realization of a functional network
#' with Bernoulli edges, according to an inner product latent process model.
#' The latent processes are generated from a \eqn{B}-spline basis with equally
#' spaced knots.
#'
#' The spline design of the functional network data (snapshot indices,
#' basis dimension) is generated using the information provided in
#' \code{spline_design}, producing a \eqn{q}-dimensional cubic
#' \eqn{B}-spline basis with equally spaced knots.
#'
#' The (\eqn{q \times d}) latent process basis coordinates \eqn{W_i}
#' for each node are generated as \eqn{q} iid Dirichlet
#' random variables with \eqn{d}-dimensional parameter
#' \code{process_options$alpha_coord} or
#' \code{rep(process_options$alpha_coord,d)} depending on the dimension
#' of \code{process_options$alpha_coord}.
#' Roughly, smaller values of \code{process_options$alpha_coord} will
#' tend to generate latent positions closer to the corners of the simplex.
#'
#' \eqn{W_i} is then rescaled so the overall network density is approximately
#' \code{process_options$density}, and the Euclidean norm of \eqn{z_i(x)}
#' never exceeds \code{1}.
#' If the density requested is too high, it will revert to the maximum density
#' under this model (\eqn{1/d}).
#' Then each latent process is given by
#' \deqn{z_{i}(x) = W_i^{T}B(x).}
#'
#' The \eqn{n \times n} symmetric adjacency matrix for
#' snapshot \eqn{k=1,...,m} has independent Bernoulli entries
#' with mean
#' \deqn{E([A_k]_{ij}) = z_i(x_k)^{T}z_j(x_k)}
#' for \eqn{i \leq j} (or \eqn{i < j} with no self loops).
#'
#'
#' @usage
#' rdpg_snapshot_bs(n,d,m,self_loops=TRUE,
#'                  spline_design,process_options)
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
#' @param process_options A list, containing additional optional arguments:
#' \describe{
#'     \item{alpha_coord}{A positive scalar, or a vector of length \eqn{d}.
#'     If it is a vector, it corresponds to the Dirichlet parameter of the
#'     basis coordinates.
#'     If is is a scalar, the basis coordinates have Dirichlet parameter
#'     \code{rep(alpha_coord,d)}. Defaults to \code{0.1}.}
#'     \item{density}{A scalar between \code{0} and \code{1}, which controls the
#'     approximate overall edge density of the resulting multiplex matrix.
#'     Defaults to \code{1/d}. If specified larger than \code{1/d}, this
#'     argument is reset to \code{1/d} and a warning is given.}
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
#' # Bernoulli edge data with B-spline latent processes, Dirichlet coordinates
#' # NOTE: for B-splines, x_max and x_min do not need to coincide with the
#' # max and min snapshot times.
#'
#' data <- rdpg_snapshot_bs(n=100,d=10,
#'                          self_loops=FALSE,
#'                          spline_design=list(q=8,
#'                                             x_vec=seq(-1,1,length.out=50),
#'                                             x_min=-1.1,x_max=1.1),
#'                          process_options=list(alpha_coord=.2,
#'                          density=1/10))
#'
#' @export
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
  if(is.null(process_options$alpha_coord)){
    process_options$alpha_coord <- .1
  }
  if(is.null(process_options$alpha_coord)){
    process_options$sigma_int <- rep(.1,d)
  }
  else{
    if(length(process_options$alpha_coord)==1){
      process_options$alpha_coord <- rep(process_options$alpha_coord,d)
    }
    else{
      if(length(process_options$alpha_coord) != d){
        stop('invalid dimension of intercept standard deviations, must be 1 or d')
      }
    }
  }
  if(is.null(process_options$density)){
    process_options$density <- 1/d
  }
  # parameter checking
  if(process_options$density > (1/d)){
    warning(paste0('requested density is too high, networks will have edge density approximately ',round(1/d,2),'\n'))
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

