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
B_func <- function(q,x_min,x_max){
  # calculate evenly spaced knots
  nknots <- q-4
  if(nknots < 0){
    return(NULL)
  }
  if(nknots == 0){
    knot_seq <- NULL
  }
  if(nknots > 0){
    knot_seq <- seq(from=x_min,to=x_max,length.out=nknots+2)[-c(1,(nknots+2))]
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

# MAKE VISIBLE
# Procrustes alignment routine
proc_align <- function(A,B,return_rot=FALSE){
  # factorize crossproduct
  temp <- svd(crossprod(A,B))
  # get rotation
  rot <- tcrossprod(temp$u,temp$v)
  A_rot <- A %*% rot
  if(return_rot){
    return(list(Ar=A_rot,rot=rot))
  }
  else{
    return(A_rot)
  }
}

# MAKE VISIBLE
# proc_align after matricizing in 3rd mode
proc_align3 <- function(A,B,return_rot=FALSE){
  matricize_align <- proc_align(t(rTensor::k_unfold(rTensor::as.tensor(B),m=3)@data),
                                t(rTensor::k_unfold(rTensor::as.tensor(A),m=3)@data),
                                return_rot=TRUE)
  rot <- matricize_align$rot
  A_rot <- rTensor::ttm(rTensor::as.tensor(A),rot,m=3)@data
  if(return_rot){
    return(list(Ar=A_rot,rot=rot))
  }
  else{
    return(A_rot)
  }
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
    orthrot <- svd(Z[,,1])$v
    temp[,,1] <- Z[,,1] %*% orthrot
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
    # B-spline new design matrix
    if(spline_design$type=='bs'){
      spline_mat_temp <- B_func(spline_design$q,
                                spline_design$x_min,
                                spline_design$x_max)(x)
    }
    # S-spline new design matrix
    else{
      m <- length(spline_design$x_vec)
      spline_mat_temp <- splines2::naturalSpline(x,
                                                 intercept=TRUE,
                                                 knots=spline_design$x_vec[2:(m-1)],
                                                 Boundary.knots=c(spline_design$x_vec[1],spline_design$x_vec[m]))
    }
    Z_temp <- coord_to_snap(W,spline_mat_temp)
    return(Z_temp)
  }
  return(Z)
}

