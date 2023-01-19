#' Functional adjacency spectral embedding
#'
#' \code{fase} fits a functional adjacency spectral embedding to snapshots
#' of (undirected) functional network data. The latent processes are fit
#' in a spline basis specified by the user, with additional options for
#' ridge penalization.
#'
#' \code{fase} finds a functional adjacency spectral embedding of an
#' \eqn{n \times n \times m} array \eqn{A} of
#' symmetric adjacency matrices on a common set of nodes, where
#' each \eqn{n \times n} slice is associated to a scalar index \eqn{x_k}
#' for \eqn{k=1,...,m}.
#' Embedding requires the specification of a latent space dimension
#' \eqn{d} and spline design information (with the argument
#' \code{spline_design}).
#'
#' \code{fase} can fit latent processes using either a cubic \eqn{B}-spline
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
#' using gradient descent on the basis coordinates \eqn{w_{i,r}}
#' of each component process
#' \deqn{z_{i,r}(x) = w_{i,r}^{T}B(x).}
#' Additional options for the fitting algorithm, including initialization,
#' can be specified by the argument \code{optim_options}.
#' For more details on the fitting and initialization algorithms, see
#' \href{https://arxiv.org/abs/2210.07491}{MacDonald et al., (2022+)},
#' Section 3.
#'
#' By default, \code{fase} will return estimates of the latent processes
#' evaluated at the snapshot indices as an \eqn{n \times d \times m} array, after
#' performing a Procrustes alignment of the consecutive snapshots.
#' This extra alignment step can be skipped.
#' \code{output_options$return_fn} provides an option to return
#' the latent processes as a function which
#' takes a vector of indices and returns the corresponding evaluations of
#' the latent process matrices.
#' \code{fase} will also return the spline design information used to fit the
#' embedding, convergence information for gradient descent, and (if specified)
#' the basis coordinates.
#'
#' When fitting with \eqn{B}-splines, \code{fase} can return a
#' network generalized cross validation criterion, described in
#' \href{https://arxiv.org/abs/2210.07491}{MacDonald et al., (2022+)},
#' Section 3.3. This criterion can be minimized to choose appropriate values
#' for \eqn{q} and \eqn{d}.
#'
#' @usage
#' fase(A,d,self_loops,spline_design,lambda,optim_options,output_options)
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
#' penalty (see Details).
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
#'     \item{init_band}{A positive scalar, the bandwidth for local averaging in
#'     the initialization algorithm. Defaults to
#'     \code{spline_design$q/(m*(spline_design$x_max - spline_design$x_min))}
#'     for \eqn{B}-spline designs
#'     and \code{10/(m*(spline_design$x_max - spline_design$x_min))}
#'     for smoothing spline designs.}
#'     \item{init_q}{A positive integer, the size of the reduced grid of embedding
#'     indices used in the initialization algorithm. Defaults to \code{spline_design$q}
#'     for \eqn{B}-spline designs
#'     and \code{10} for smoothing spline designs.}
#' }
#' @param output_options A list, containing additional optional arguments controlling
#' the output of \code{fase}.
#' \describe{
#'     \item{align_output}{A Boolean, if \code{TRUE}, the returned latent processes
#'     have been aligned according to a Procrustes alignment which minimizes
#'     (in terms of Frobenius norm) the overall discrepancies between consecutive
#'     snapshots. Defaults to \code{TRUE}.}
#'     \item{return_fn}{A Boolean, if \code{TRUE}, the latent processes are
#'     returned as a function which
#'     takes a vector of indices and returns the corresponding evaluations of
#'     the latent process matrices. Otherwise, the latent processes are returned as
#'     an \eqn{n \times d \times m} array. Defaults to \code{FALSE}.}
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
#' \item{Z}{Either an \eqn{n \times d \times m} array containing the latent process embedding
#' evaluated at the indices in \code{spline_design$x_vec}, or a function which
#' takes a vector of indices and returns the corresponding evaluations of
#' the latent process matrices, depending on \code{output_options$return_fn}.}
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
#' \item{convergence}{An integer convergence code, \code{0} if
#' gradient descent converged in fewer than \code{optim_options$K_max} iterations,
#' \code{1} otherwise.}
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
#' fit_bs <- fase(data$A,d=2,self_loops=FALSE,
#'                spline_design=list(type='bs',q=9,x_vec=data$spline_design$x_vec),
#'                optim_options=list(eps=1e-4,K_max=40,init_q=12),
#'                output_options=list(return_coords=TRUE))
#'
#' # fase fit with smoothing spline design
#' fit_ss <- fase(data$A,d=2,self_loops=FALSE,
#'                spline_design=list(type='ss',x_vec=data$spline_design$x_vec),
#'                lambda=.5,
#'                optim_options=list(eta=1e-4,K_max=40,verbose=FALSE),
#'                output_options=list(align_output=FALSE,return_fn=TRUE))
#'
#' #NOTE: both models fit with small optim_options$K_max=40 for demonstration
#'
#' @export
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
  if(spline_design$type=='bs'){
    if((optim_options$init_q < spline_design$q) || (optim_options$init_q < 4)){
      stop('invalid choice of initial q (too small)')
    }
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
                     ridge_pen=lambda/n,
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
