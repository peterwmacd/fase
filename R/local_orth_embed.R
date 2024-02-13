# local_orth_embed: large internal function

# local average spectral embedding for initialization
local_orth_embed <- function(A,d,
                            spline_design,
                            L,M){
  # design/dimensions
  n <- dim(A)[1]
  m <- length(spline_design$x_vec)
  # partition indices
  groups <- ceiling(L*1:m / m)
  # dictionary of local embeddings
  Z_tilde <- array(0,c(n,d,L))
  for(ll in 1:L){
    iloc <- midM(which(groups==ll),M)
    temp <- ase(apply(A[,,iloc],c(1,2),mean),d,positive=TRUE)
    if(ll > 1){
      if(is.matrix(temp)){
        Z_tilde[,,ll] <- proc_align(temp,Z_tilde[,,ll-1])
      }
      else{
        Z_tilde[,,ll] <- temp * sign(sum(temp*Z_tilde[,,ll-1]))
      }
    }
    else{
      Z_tilde[,,ll] <- temp
    }
  }
  # populate full processes
  Z_init <- array(0,c(n,d,m))
  for(kk in 1:m){
    Z_init[,,kk] <- Z_tilde[,,groups[kk]]
  }
  # get initial coordinates
  W_init <- B_smooth(Z_init,spline_design$spline_mat)
  return(W_init)
}
