# kern_orth_embed: large internal function

# kernel smoothed embedding for initialization
kern_orth_embed <- function(A,d,
                            spline_design,
                            init_q,band){
  # design/dimensions
  if(spline_design$type=='ss'){
    B_mat <- B_func(init_q,
                    spline_design$x_min,
                    spline_design$x_max)(spline_design$x_vec)
    xq_vec <- spline_design$x_vec[apply(B_mat,2,which.max)]
    Bq_mat <- B_func(init_q,
                     spline_design$x_min,
                     spline_design$x_max)(xq_vec)

  }
  else{
    if(init_q==spline_design$q){
      B_mat <- spline_design$spline_mat
      xq_vec <- spline_design$x_vec[apply(spline_design$spline_mat,2,which.max)]
    }
    else{
      B_mat <- B_func(init_q,
                      spline_design$x_min,
                      spline_design$x_max)(spline_design$x_vec)
      xq_vec <- seq(spline_design$x_min,
                    spline_design$x_max,
                    length.out=init_q)
    }
    Bq_mat <- B_func(init_q,
                     spline_design$x_min,
                     spline_design$x_max)(xq_vec)
  }
  # increase d for kernel
  # embed each slice and align it to the true value
  Zq_hat <- array(0,c(dim(A)[1],d,init_q))
  # embed the first slice
  close <- which(abs(xq_vec[1]-spline_design$x_vec)<band)
  # align to place in the positive orthant
  Zq_hat[,,1] <- ase(apply(A[,,close],c(1,2),mean),d,positive=TRUE)
  # backwards from center
  for(ll in 2:init_q){
    close <- which(abs(xq_vec[ll]-spline_design$x_vec)<band)
    # get ASE
    temp <- ase(apply(A[,,close],c(1,2),mean),d,positive=TRUE)
    # get permutation
    if(is.matrix(temp)){
      Zq_hat[,,ll] <- proc_align(temp,Zq_hat[,,ll-1])
    }
    else{
      Zq_hat[,,ll] <- temp * sign(sum(temp*Zq_hat[,,ll-1]))
    }
  }
  # smooth the slices against B_mat
  Wq_hat <- B_smooth(Zq_hat,Bq_mat)
  # return translated coordinate estimate
  if(spline_design$type=='ss'){
    W_hat <- B_translate(Wq_hat,B_mat,spline_design$spline_mat)
  }
  else{
    if(init_q==spline_design$q){
      W_hat <- Wq_hat
    }
    else{
      W_hat <- B_translate(Wq_hat,B_mat,spline_design$spline_mat)
    }
  }
  return(W_hat)
}
