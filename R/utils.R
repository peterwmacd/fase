# utility functions

# hollowization for square matrices
hollowize <- function(M){
  M - diag(diag(M))
}

# hollowization for 3D arrays
hollowize3 <- function(A){
  array(apply(A,3,hollowize),dim(A))
}

# helper which converts W and B_mat to an array of Theta = EA
WB_to_Theta <- function(W,B_mat,self_loops){
  dW <- dim(W)
  m <- nrow(B_mat)
  Theta <- array(apply(array(apply(W,3,function(x){x %*% t(B_mat)}),c(dW[1],m,dW[3])),2,tcrossprod),
                 c(dim(W)[1],dim(W)[1],m))
  if(self_loops){
    return(Theta)
  }
  else{
    return(hollowize3(Theta))
  }
}

# helper which converts Z to an array of Theta = EA
Z_to_Theta <- function(Z,self_loops){
  Theta <- array(apply(Z,3,tcrossprod),
                 c(dim(Z)[1],dim(Z)[1],dim(Z)[3]))
  if(self_loops){
    return(Theta)
  }
  else{
    return(hollowize3(Theta))
  }
}

# cubic B-spline evaluation function
B_func <- function(q,x_min,x_max,x_vec=NULL){
  # calculate evenly spaced knots
  nknots <- q-4
  if(nknots < 0){
    return(NULL)
  }
  if(nknots == 0){
    knot_seq <- NULL
  }
  if(nknots > 0){
    if(is.null(x_vec)){
      # pass x_vec=NULL for evenly spaced knots
      knot_seq <- seq(from=x_min,to=x_max,length.out=nknots+2)[-c(1,(nknots+2))]
    }
    else{
      # alternative with quantiles of x_vec
      knot_seq <- x_vec[floor(seq(from=1,to=length(x_vec),length.out=nknots+2))[-c(1,(nknots+2))]]
    }
  }
  B <- function(x){
    in_int <- (x >= x_min & x <= x_max)
    out <- matrix(0,length(x),q)
    out[in_int,] <- splines2::bSpline(x[in_int],
                                      knots=knot_seq,
                                      degree=3,
                                      intercept=TRUE,
                                      Boundary.knots=c(x_min,x_max))
    return(out)
  }
  return(B)
}

# adjacency spectral embedding
ase <- function(M,d,positive=TRUE){
  if(d==0){
    X.hat <- matrix(0,nrow(M),1)
  }
  else{
    temp <- RSpectra::eigs_sym(M,d,which=ifelse(positive,'LA','LM'))
    perm <- order(abs(temp$values),decreasing=TRUE)
    if(d>1){
      X.hat <- temp$vectors[,perm] %*% diag(sqrt(abs(temp$values[perm])))
    }
    else{
      X.hat <- temp$vectors[,perm]*sqrt(abs(temp$values[perm]))
    }
  }
  attr(X.hat,'signs') <- sign(temp$values[perm])
  return(X.hat)
}

# converting (n x q x d) coordinates to (n x d x m) snapshots
coord_to_snap <- function(W,B_mat){
  aperm(rTensor::ttm(rTensor::as.tensor(W),B_mat,m=2)@data,c(1,3,2))
}

#' Procrustes alignment
#'
#' \code{proc_align} orthogonally transforms the columns of a matrix \eqn{A} to
#' find the best approximation (in terms of Frobenius norm) to a
#' second matrix \eqn{B}. Optionally, it may also return the optimal transformation
#' matrix.
#'
#' @usage
#' proc_align(A,B,return_orth=FALSE)
#'
#' @param A An \eqn{n \times d} matrix.
#' @param B An \eqn{n \times d} matrix.
#' @param return_orth A Boolean which specifies whether to return the
#' orthogonal transformation.
#' Defaults to \code{FALSE}.
#'
#' @return If \code{return_orth} is \code{FALSE}, returns the \eqn{n \times d}
#' matrix resulting from applying the optimal aligning transformation to
#' the columns of \code{A}.
#' Otherwise, returns a list with two entries:
#' \item{Ao}{The \eqn{n \times d}
#' matrix resulting from applying the optimal aligning transformation to the
#' columns of \code{A}.}
#' \item{orth}{The \eqn{d \times d} optimal aligning orthogonal transformation
#' matrix.}
#'
#'
#' @export
proc_align <- function(A,B,return_orth=FALSE){
  # factorize crossproduct
  temp <- svd(crossprod(A,B))
  # get transformation
  orth <- tcrossprod(temp$u,temp$v)
  A_orth <- A %*% orth
  if(return_orth){
    return(list(Ao=A_orth,orth=orth))
  }
  else{
    return(A_orth)
  }
}

#' Procrustes alignment for 3-mode tensors
#'
#' \code{proc_align3} applies one orthogonal transformation
#' to the columns of each of the \eqn{n \times d} slices of an
#' \eqn{n \times d \times m} array \eqn{A} to
#' find the best approximation (in terms of matrix Frobenius norm, averaged
#' over the  \eqn{n \times d} slices) to a
#' second \eqn{n \times d \times m} array \eqn{B}.
#' Optionally, it may also return the optimal transformation
#' matrix.
#'
#' @usage
#' proc_align3(A,B,return_orth=FALSE)
#'
#' @param A An \eqn{n \times d \times m} array.
#' @param B An \eqn{n \times d \times m} array.
#' @param return_orth A Boolean which specifies whether to return the
#' orthogonal transformation.
#' Defaults to \code{FALSE}.
#'
#' @return If \code{return_orth} is \code{FALSE}, returns the \eqn{n \times d \times m}
#' array resulting from applying the optimal aligning transformation to
#' the columns of the \eqn{n \times d} slices of \code{A}.
#' Otherwise, returns a list with two entries:
#' \item{Ao}{The \eqn{n \times d}
#' matrix resulting from applying the optimal aligning transformation to the
#' columns of the \eqn{n \times d} slices of \code{A}.}
#' \item{orth}{The \eqn{d \times d} optimal aligning orthogonal transformation matrix.}
#'
#'
#' @export
proc_align3 <- function(A,B,return_orth=FALSE){
  matricize_align <- proc_align(t(rTensor::k_unfold(rTensor::as.tensor(B),m=2)@data),
                                t(rTensor::k_unfold(rTensor::as.tensor(A),m=2)@data),
                                return_orth=TRUE)
  orth <- matricize_align$orth
  A_orth <- rTensor::ttm(rTensor::as.tensor(A),orth,m=2)@data
  if(return_orth){
    return(list(Ao=A_orth,orth=orth))
  }
  else{
    return(A_orth)
  }
}

#' Slicewise Procrustes alignment for 3-mode tensors
#'
#' \code{proc_align_slicewise3} applies an orthogonal transformation
#' to the columns of each of the \eqn{n \times d} slices of an
#' \eqn{n \times d \times m} array \eqn{A} to
#' find the best approximation (in terms of matrix Frobenius norm) to
#' the corresponding \eqn{n \times d} slice of a
#' second \eqn{n \times d \times m} array \eqn{B}.
#'
#' @usage
#' proc_align_slicewise3(A,B)
#'
#' @param A An \eqn{n \times d \times m} array.
#' @param B An \eqn{n \times d \times m} array.
#'
#' @return Returns the \eqn{n \times d \times m}
#' array resulting from applying the optimal aligning transformations to
#' the columns of the \eqn{n \times d} slices of \code{A}.
#'
#'
#' @export
proc_align_slicewise3 <- function(A,B){
  A_rot <- array(NA,dim(A))
  for(ii in 1:dim(A)[3]){
    A_rot[,,ii] <- proc_align(A[,,ii],B[,,ii])
  }
  return(A_rot)
}


# helper to convert (n x d x m) snapshots to (n x q x d) smoothed coordinates
B_smooth <- function(Z,B_mat){
  aperm(rTensor::ttm(rTensor::as.tensor(Z),solve(crossprod(B_mat),t(B_mat)),m=3)@data,c(1,3,2))
}

# translating coordinates between bases
# assume the (n x q1 x d) coordinates are in B1 and we want them (n x q2 x d) in B2
B_translate <- function(W,B1,B2){
  rTensor::ttm(rTensor::as.tensor(W),solve(crossprod(B2),crossprod(B2,B1)),m=2)@data
}

# Procrustes aligning consecutive process evaluations
Z_align_proc <- function(Z,orth1=FALSE){
  # dimensions
  m <- dim(Z)[3]
  temp <- array(NA,dim(Z))
  if(orth1){
    orth <- svd(Z[,,1])$v
    temp[,,1] <- Z[,,1] %*% orth
  }
  else{
    temp[,,1] <- Z[,,1]
  }
  if(m >1){
    for(kk in 2:m){
      temp[,,kk] <- proc_align(Z[,,kk],temp[,,kk-1])
    }
  }
  return(temp)
}

# ridge penalty matrix for S-splines
ss_ridge <- function(spline_design){
  # extended natural spline matrix for approximating
  m <- length(spline_design$x_vec)
  m_ext <- max(1e3+1,m)
  x_vec_ext <- seq(spline_design$x_min,spline_design$x_max,length.out=m_ext)
  # approximation of ridge matrix
  Spp_mat <- splines2::naturalSpline(x_vec_ext,
                                     intercept=TRUE,
                                     knots=spline_design$x_vec[2:(m-1)],
                                     Boundary.knots=c(spline_design$x_min,spline_design$x_max),
                                     derivs=2)
  # numerical integration to approximate Omega_N
  Spp_cross_1 <- crossprod(Spp_mat[-1,])*mean(diff(x_vec_ext))
  Spp_cross_m <- crossprod(Spp_mat[-m_ext,])*mean(diff(x_vec_ext))
  Omega_S <- .5*(Spp_cross_1 + Spp_cross_m)
  return(Omega_S)
}

# convert to functional output
coord_to_func <- function(W,spline_design){
  # define a function of a possible vector x
  Z <- function(x){
    Z_temp <- array(0,c(dim(W)[1],dim(W)[3],length(x)))
    x_sub <- ((x >= spline_design$x_min) & (x <= spline_design$x_max))
    if(sum(x_sub) > 0){
      # B-spline new design matrix
      if(spline_design$type=='bs'){
        spline_mat_temp <- B_func(spline_design$q,
                                  spline_design$x_min,
                                  spline_design$x_max,
                                  spline_design$x_vec)(x[x_sub])
      }
      # S-spline new design matrix
      else{
        m <- length(spline_design$x_vec)
        spline_mat_temp <- splines2::naturalSpline(x[x_sub],
                                                   intercept=TRUE,
                                                   knots=spline_design$x_vec[2:(m-1)],
                                                   Boundary.knots=c(spline_design$x_vec[1],spline_design$x_vec[m]))
      }
      Z_temp[,,x_sub] <- coord_to_snap(W,spline_mat_temp)
    }
    return(Z_temp)
  }
  return(Z)
}

# helper for subsetting the middle elements of a vector
midM <- function(v,M){
  L <- length(v)
  diff <- L - M
  if(diff <= 0){
    return(v)
  }
  else{
    if((diff %% 2)==0){
      d1 <- 1:(diff/2)
      d2 <- (L-(diff/2)+1):L
    }
    else{
      d1 <- 1:floor(diff/2)
      d2 <- (L-floor(diff/2)):L
    }
    return(v[-c(d1,d2)])
  }
}

# (ported over from multiness package)
# robust noise estimation for a low-rank matrix with the MAD estimator
# of Gavish and Donoho
# This is a simpler adaptation of the function denoiseR::estim_sigma
# as that package has non-trivial dependencies, and
# its maintenance status is unclear
estim_sigma_mad <- function(M){
  # center columns
  Mc <- scale(M, scale = F)
  # dimensions
  n = nrow(Mc)
  p = ncol(Mc)
  # auxiliary quantities
  beta <- min(n, p)/max(n, p)
  lambdastar <- sqrt(2 * (beta + 1) + 8 * beta/((beta +
                                                   1 + (sqrt(beta^2 + 14 * beta + 1)))))
  wbstar <- 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta +
    1.43
  # sigma estimate
  sigma <- stats::median(svd(Mc)$d)/(sqrt(max(n, p)) * (lambdastar/wbstar))
  return(sigma)
}
