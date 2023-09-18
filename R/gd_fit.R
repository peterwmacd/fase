# gd_fit: large internal function

# residual matrix
resid <- function(A,W,B_mat,self_loops){
  A - WB_to_Theta(W,B_mat,self_loops)
}

# objective and gradient calculations
objgrad <- function(A,W,B_mat,B_pos,self_loops){
  # dims
  dW <- dim(W)
  m <- nrow(B_mat)
  # residual
  R <- resid(A,W,B_mat,self_loops)
  # objective
  obj <- sum(R^2)
  # gradient
  # allocate space
  temp <- array(0,dW)
  # sum over slices
  for(kk in 1:m){
    qpos <- B_pos[kk,]
    nqpos <- sum(qpos)
    if(nqpos > 1){
      if(dW[3] > 1){
        temp[,qpos,] <- temp[,qpos,] + array(apply(W,3,function(x){colSums(t(R[,,kk] %*% x[,qpos])*B_mat[kk,qpos]) %o% B_mat[kk,qpos]}),c(dW[1],nqpos,dW[3]))
      }
      else{
        temp[,qpos,] <- temp[,qpos,] + matrix(apply(W,3,function(x){colSums(t(R[,,kk] %*% x[,qpos])*B_mat[kk,qpos]) %o% B_mat[kk,qpos]}),dW[1],nqpos)
      }
    }
    else{
      if(dW[3] > 1){
        temp[,qpos,] <- temp[,qpos,] + matrix(apply(W,3,function(x){(R[,,kk] %*% x[,qpos])}),dW[1],dW[3])
      }
      else{
        temp[,qpos,] <- temp[,qpos,] + c(apply(W,3,function(x){(R[,,kk] %*% x[,qpos])}))
      }
    }
  }
  return(list(obj=obj,grad=(-4)*temp))
}

# objective and gradient for generalized ridge component
objgrad_ridge <- function(W,G,lambda){
  # matrix calculations
  wg2 <- apply(W,c(1,3),function(v){
    temp <- G %*% v
    return(c(sum(v*temp),temp))
  })
  # objective
  obj <- lambda*sum(wg2[1,,])
  # gradient
  grad <- 2*lambda*aperm(wg2[-1,,,drop=FALSE],c(2,1,3))
  return(list(obj=obj,grad=grad))
}

# main gradient descent fitting function
gd_fit <- function(A,self_loops,
                   W_init,B_mat,
                   eta,eps,K_max,
                   ridge=FALSE,ridge_pen=NULL,ridge_mat=NULL,
                   verbose=FALSE){
  # precalculate B_pos, dimension
  B_pos <- apply(B_mat,2,function(v){v > 0})
  nnm <- ifelse(self_loops,length(A),dim(A)[1]*(dim(A)[1]-1)*dim(A)[3])
  # current step size
  eta_curr <- eta
  # initialize
  W_old <- W_init
  # starting objective
  if(ridge){
    objgrad_old <- objgrad(A,W_old,B_mat,B_pos,self_loops)
    objgrad_ridge_old <- objgrad_ridge(W_old,ridge_mat,ridge_pen)
    obj_old <- objgrad_old$obj + objgrad_ridge_old$obj
  }
  else{
    objgrad_old <- objgrad(A,W_old,B_mat,B_pos,self_loops)
    obj_old <- objgrad_old$obj
  }
  converged <- 0
  # gradient descent routine
  K <- 0
  while(K <= K_max){
    K <- K+1
    if(verbose){
      cat(paste0('iter ',K,'\n'))
    }
    # leave loop if the step size gets too small
    if(eta_curr < 1e-64){
      break
    }
    # update iterate and objective
    if(ridge){
      W_new <- W_old - eta_curr*(objgrad_old$grad + objgrad_ridge_old$grad)
    }
    else{
      W_new <- W_old - eta_curr*objgrad_old$grad
    }
    # check convergence
    if(ridge){
      objgrad_new <- objgrad(A,W_new,B_mat,B_pos,self_loops)
      objgrad_ridge_new <- objgrad_ridge(W_new,ridge_mat,ridge_pen)
      obj_new <- objgrad_new$obj + objgrad_ridge_new$obj
    }
    else{
      objgrad_new <- objgrad(A,W_new,B_mat,B_pos,self_loops)
      obj_new <- objgrad_new$obj
    }
    if(verbose){
      cat(paste0('RMSE: ',round(sqrt(obj_new/nnm),4),'\n'))
    }
    # check convergence
    crit <- (obj_old - obj_new)/obj_old
    if(verbose){
      cat(paste0('convergence: ',crit,'\n'))
    }
    if(crit >= 0){
      if(crit < eps){
        converged <- 1
        break
      }
      else{
        W_old <- W_new
        objgrad_old <- objgrad_new
        if(ridge){
          objgrad_ridge_old <- objgrad_ridge_new
        }
        obj_old <- obj_new
        eta_curr <- min(eta,2*eta_curr)
      }
    }
    else{
      eta_curr <- .5*eta_curr
    }
  }
  return(list(W_hat=W_new,K=K,converged=converged,obj=obj_new))
}



