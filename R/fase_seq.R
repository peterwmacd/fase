#' Functional adjacency spectral embedding (sequential algorithm)
#'
#' \code{fase_seq} fits a functional adjacency spectral embedding to snapshots
#' of (undirected) functional network data, with each
#' of the \eqn{d} latent dimensions fit sequentially. The latent processes are fit
#' in a spline basis specified by the user, with additional options for
#' ridge penalization.
#'
#' Note that \code{fase_seq} is a wrapper for \code{fase}. When \eqn{d=1},
#' \code{fase_seq} coincides with \code{fase}.
#'
#' \code{fase_seq} finds a functional adjacency spectral embedding of an
#' \eqn{n \times n \times m} array \eqn{A} of
#' symmetric adjacency matrices on a common set of nodes, where
#' each \eqn{n \times n} slice is associated to a scalar index \eqn{x_k}
#' for \eqn{k=1,...,m}.
#' Embedding requires the specification of a latent space dimension
#' \eqn{d} and spline design information (with the argument
#' \code{spline_design}).
#'
#' \code{fase_seq} can fit latent processes using either a cubic \eqn{B}-spline
#' basis with
#' equally spaced knots, or a natural cubic spline basis with a second
#' derivative (generalized ridge) smoothing penalty: a smoothing spline.
#' To fit with a \eqn{B}-spline design (\code{spline_design$type = 'bs'}),
#' one must minimally provide a basis
#' dimension \eqn{q} of at least \code{4} and at most \eqn{m}.
#'
#' When fitting with a smoothing spline design, the generalized ridge
#' penalty is scaled by
#' \eqn{\lambda/n}, where \eqn{\lambda} is specified by the argument \code{lambda}.
#' see \href{https://arxiv.org/abs/2210.07491}{MacDonald et al., (2022+)},
#' Appendix E for more details.
#' \code{lambda} can also be used to introduce a ridge penalty on the
#' basis coordinates when fitting with \eqn{B}-splines.
#'
#' Fitting minimizes a least squares loss,
#' using gradient descent (Algorithm 1) on the basis coordinates \eqn{w_{i,r}}
#' of each component process
#' \deqn{z_{i,r}(x) = w_{i,r}^{T}B(x).}
#' Additional options for the fitting algorithm, including initialization,
#' can be specified by the argument \code{optim_options}.
#' For more details on the fitting and initialization algorithms, see
#' \href{https://arxiv.org/abs/2210.07491}{MacDonald et al., (2022+)},
#' Section 3.
#'
#' By default, \code{fase_seq} will return estimates of the latent processes
#' evaluated at the snapshot indices as an \eqn{n \times d \times m} array, after
#' performing a Procrustes alignment of the consecutive snapshots.
#' This extra alignment step can be skipped.
#' \code{fase_seq} will also return the spline design information used to fit the
#' embedding, convergence information for gradient descent, and (if specified)
#' the basis coordinates.
#'
#' When fitting with \eqn{B}-splines, \code{fase_seq} can return a
#' network generalized cross validation criterion, described in
#' \href{https://arxiv.org/abs/2210.07491}{MacDonald et al., (2022+)},
#' Section 3.3. This criterion can be minimized to choose appropriate values
#' for \eqn{q} and \eqn{d}.
#'
#' @usage
#' fase_seq(A,d,self_loops,spline_design,lambda,optim_options,output_options)
#'
#' @param A An \eqn{n \times n \times m} array containing
#' the snapshots of the functional network.
#' @param d A positive integer, the number of latent space dimensions of the
#' functional embedding.
#' @param self_loops A Boolean, if \code{FALSE}, all diagonal entries are ignored in
#' optimization. Defaults to \code{TRUE}.
#' @param spline_design A list, containing the spline design information.
#' For fitting with a \eqn{B}-spline design (the default):
#' \describe{
#'     \item{type}{The string \code{'bs'}.}
#'     \item{q}{A positive integer, the dimension of the \eqn{B}-spline basis.}
#'     \item{x_vec}{A vector, the snapshot evaluation indices for the data.
#'     Defaults to an equally spaced vector of length \eqn{m} from \code{0}
#'     to \code{1}.}
#'     \item{x_max}{A scalar, the maximum of the index space. Defaults to
#'     \code{max(spline_design$x_vec)}.}
#'     \item{x_min}{A scalar, the minimum of the index space. Defaults to
#'     \code{min(spline_design$x_vec)}.}
#'     \item{spline_matrix}{An \eqn{m \times q} matrix, the B-spline basis
#'     evaluated at the snapshot indices. If not specified, it will be
#'     calculated internally.}
#'     \item{ridge_mat}{The \eqn{m \times m} matrix for the generalized
#'     ridge penalty. If \code{lambda}\eqn{> 0},
#'     defaults to \code{diag(m)}.}
#' }
#' For fitting with a smoothing spline design:
#' \describe{
#'     \item{type}{The string \code{'ss'}.}
#'     \item{x_vec}{A vector, the snapshot evaluation indices for the data.
#'     Defaults to an equally spaced vector of length \eqn{m} from \code{0}
#'     to \code{1}.}
#'     \item{x_max}{A scalar, the maximum of the index space. Defaults to
#'     \code{max(spline_design$x_vec)}.}
#'     \item{x_min}{A scalar, the minimum of the index space. Defaults to
#'     \code{min(spline_design$x_vec)}.}
#'     \item{spline_matrix}{An \eqn{m \times m} matrix, the
#'     natural cubic spline basis
#'     evaluated at the snapshot indices. If not specified, it will be
#'     calculated internally.}
#'     \item{ridge_mat}{The \eqn{m \times m} matrix for the generalized
#'     ridge penalty. Defaults to the second derivatives of the natural cubic spline
#'     basis evaluated at the snapshot indices.}
#' }
#' @param lambda A positive scalar, the scale factor for the generalized ridge
#' penalty (see Details). Defaults to \code{0}.
#' @param optim_options A list, containing additional optional arguments controlling
#' the gradient descent algorithm.
#' \describe{
#'     \item{eps}{A positive scalar, the convergence threshold for gradient
#'     descent in terms of relative change in objective value.
#'     Defaults to \code{1e-5}.}
#'     \item{eta}{A positive scalar, the step size for gradient descent.
#'     Defaults to \code{1/(n*m)}.}
#'     \item{K_max}{A positive integer, the maximum iterations for gradient
#'     descent. Defaults to \code{2e3}.}
#'     \item{verbose}{A Boolean, if \code{TRUE}, console output will provide
#'     updates on the progress of gradient descent. Defaults to
#'     \code{FALSE}.}
#'     \item{init_W}{A 3-dimensional array
#'     containing initial basis coordinates for gradient descent. Dimension should be
#'     \eqn{n \times}\code{spline_design$q}\eqn{ \times d} for \eqn{B}-spline designs,
#'     and \eqn{n \times m \times d} for smoothing spline designs. If included,
#'     \code{init_M}, \code{init_L} and \code{init_sigma} are ignored.}
#'     \item{init_sigma}{A positive scalar, the estimated edge dispersion parameter to calibrate
#'     initialization. If not provided, it is either estimated using the robust method proposed by
#'     Gavish and Donoho (2014) for weighted edge networks, or set to a default value \code{0.5}
#'     for binary edge networks.}
#'     \item{init_L}{A positive integer, the number of contiguous groups used for initialization.
#'     Defaults to the floor of \eqn{(2nm/\texttt{init\_sigma}^2)^{1/3}}.}
#'     \item{init_M}{A positive integer, the number of snapshots averaged in each group for
#'     initialization. Defaults use all snapshots.}
#' }
#' @param output_options A list, containing additional optional arguments controlling
#' the output of \code{fase}.
#' \describe{
#'     \item{return_coords}{A Boolean, if \code{TRUE}, the basis coordinates for
#'     each latent process component are also returned as an array.
#'     Defaults to \code{FALSE}.}
#'     \item{return_ngcv}{A Boolean, if \code{TRUE} and \code{spline_design$type=='bs'},
#'      the network generalized cross validation criterion is returned.
#'     Defaults to \code{TRUE}.}
#' }
#'
#' @return A list is returned with the functional adjacency spectral embedding,
#' the spline design information, and some additional optimization
#' output:
#' \item{Z}{An \eqn{n \times d \times m} array containing the latent process embedding
#' evaluated at the indices in \code{spline_design$x_vec}.}
#' \item{W}{For \eqn{B}-spline designs, an \eqn{n \times q \times d} array; or for
#' smoothing spline designs, an \eqn{n \times m \times d} array of estimated basis
#' coordinates. If \code{output_options$return_coords} is \code{FALSE},
#' this is not returned.}
#' \item{spline_design}{A list, describing the spline design:
#' \describe{
#'     \item{type}{A string, either \code{'bs'} or \code{'ss'}.}
#'     \item{q}{A positive integer, the dimension of the \eqn{B}-spline basis.
#'     Only returned for \eqn{B}-spline designs.}
#'     \item{x_vec}{A vector, the snapshot evaluation indices for the data.}
#'     \item{x_max}{A scalar, the maximum of the index space.}
#'     \item{x_min}{A scalar, the minimum of the index space.}
#'     \item{spline_matrix}{For \eqn{B}-spline designs, an \eqn{m \times q} matrix;
#'     or for smoothing spline designs, an \eqn{m \times m} matrix, the basis
#'     evaluated at the snapshot indices.}
#'     \item{ridge_matrix}{An \eqn{m \times m} matrix used in the generalized
#'     ridge penalty. Only returned for \code{lambda > 0}.}
#' }
#' }
#' \item{ngcv}{A scalar, the network generalized cross validation criterion
#' (see Details). Only returned for \eqn{B}-spline designs and when
#' \code{output_options$return_ngcv} is \code{TRUE}.}
#' \item{K}{A positive integer, the number of iterations run in
#' gradient descent.}
#' \item{converged}{An integer convergence code, \code{1} if
#' gradient descent converged in fewer than \code{optim_options$K_max} iterations,
#' \code{0} otherwise.}
#'
#' @examples
#' # Gaussian edge data with sinusoidal latent processes
#' set.seed(1)
#' data <- gaussian_snapshot_ss(n=50,d=2,
#'                              x_vec=seq(0,1,length.out=50),
#'                              self_loops=FALSE,sigma_edge=4)
#'
#'
#' # fase fit with B-spline design
#' fit_bs <- fase_seq(data$A,d=2,self_loops=FALSE,
#'                    spline_design=list(type='bs',q=9,x_vec=data$spline_design$x_vec),
#'                    optim_options=list(eps=1e-4,K_max=40),
#'                    output_options=list(return_coords=TRUE))
#'
#' # fase fit with smoothing spline design
#' fit_ss <- fase_seq(data$A,d=2,self_loops=FALSE,
#'                    spline_design=list(type='ss',x_vec=data$spline_design$x_vec),
#'                    lambda=.5,
#'                    optim_options=list(eta=1e-4,K_max=40,verbose=FALSE))
#'
#' #NOTE: both models fit with small optim_options$K_max=40 for demonstration
#'
#' @export
fase_seq <- function(A,d,self_loops=TRUE,
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
    # populate
    # B-spline design
    if(is.null(spline_design$spline_mat)){
      spline_design$spline_mat <- B_func(spline_design$q,
                                         spline_design$x_min,
                                         spline_design$x_max,
                                         spline_design$x_vec)(spline_design$x_vec)
    }
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
      if(is.null(spline_design$spline_mat)){
        spline_design$spline_mat <- splines2::naturalSpline(spline_design$x_vec,
                                                            intercept=TRUE,
                                                            knots=spline_design$x_vec[2:(m-1)],
                                                            Boundary.knots=c(spline_design$x_vec[1],spline_design$x_vec[m]))
      }
      # ridge matrix (see utils.R)
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
  # initialization parameters
  if(is.null(optim_options$init_W)){
    if(is.null(optim_options$init_M)){
      if(is.null(optim_options$init_sigma)){
        # check if binary (checks if upper triangle of first snapshot takes 2 unique values)
        if(length(unique(A[,,1][upper.tri(A[,,1])]))==2){
          optim_options$init_sigma <- 0.5
        }
        else{
          optim_options$init_sigma <- mean(apply(A,3,estim_sigma_mad))
        }
      }
      optim_options$init_M <- Inf
    }
    if(is.null(optim_options$init_L)){
      if(is.null(optim_options$init_sigma)){
        # check if binary (checks if upper triangle of first snapshot takes 2 unique values)
        if(length(unique(A[,,1][upper.tri(A[,,1])]))==2){
          optim_options$init_sigma <- 0.5
        }
        else{
          optim_options$init_sigma <- mean(apply(A,3,estim_sigma_mad))
        }
      }
      optim_options$init_L <- min(max(floor((2*n*m / (optim_options$init_sigma^2))^(1/3)),1),m)
    }
    # consistency of parameters
    if(optim_options$init_L > m){
      stop('invalid choice of initial L (too large)')
    }
    init_given <- FALSE
  }
  else{
    # check dimensions
    if(dim(optim_options$init_W)[1] != n){
      stop('invalid dimensions of initial coordinates (number of nodes)')
    }
    if(dim(optim_options$init_W)[2] != ncol(spline_design$spline_mat)){
      stop('invalid dimensions of initial coordinates (basis dimension)')
    }
    # store and reset
    init_given <- TRUE
    init_W <- optim_options$init_W
    optim_options$init_W <- NA
  }
  # output parameter checking
  if(is.null(output_options$return_coords)){
    output_options$return_coords <- TRUE
    final_coords <- FALSE
  }
  else{
    final_coords <- output_options$return_coords
    output_options$return_coords <- TRUE
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
  # initialize some parameters
  out$Z <- NA
  out$W <- NA
  out$spline_design <- spline_design
  out$ngcv <- NA
  out$K <- rep(NA,d)
  out$converged <- rep(NA,d)
  # sequential fitting
  A_seq <- A
  Z_seq <- array(NA,c(n,d,m))
  W_seq <- array(NA,c(n,ncol(spline_design$spline_mat),d))
  # loop over dimensions
  for(rr in 1:d){
    # check initializer and maybe pass it into fase
    # initialization
    if(init_given){
      optim_options$W_init <- array(init_W[,,rr],
                                    c(n,ncol(spline_design$spline_mat),1))
    }
    # run fase
    temp_fase <- fase(A_seq,1,self_loops,
                      spline_design = spline_design,
                      lambda = lambda,
                      optim_options = optim_options,
                      output_options = output_options)
    # update A
    A_seq <- A_seq - Z_to_Theta(temp_fase$Z,self_loops=self_loops)
    # update Z
    Z_seq[,rr,] <- temp_fase$Z[,1,]
    # update W
    W_seq[,,rr] <- temp_fase$W[,,1]
    # update ngcv
    out$ngcv <- temp_fase$ngcv
    # update K
    out$K[rr] <- temp_fase$K
    # update convergence
    out$converged[rr] <- temp_fase$converged
  }
  # estimated processes
  # process snapshots
  out$Z <- Z_seq
  # check if
  if(final_coords){
    out$W <- W_seq
  }
  # NGCV criterion adjustment
  if(output_options$return_ngcv){
    n_node <- ifelse(self_loops,n*m,(n-1)*m)
    out$ngcv <- out$ngcv + log(1 - ((2*spline_design$q)/n_node)) - log(1 - ((2*spline_design$q*d)/n_node))
  }
  # return
  return(out)
}
