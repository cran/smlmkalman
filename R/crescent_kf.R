#####
#Code written by Andrew Buist, updated 16/05/22
#A kalman filter for the tracking of 2-d coordinates without a supplied time dimension,
#inferring point likelihood by local, crescent-shaped regions for the z-update step.
#Use cases: tracking dense, time-absent data; image segmentation of filamentous structures
#####

crescent_kf = function(x, p_length = 1000, x_hat_s = c(0,0), sigma_s = c(10,0,0,10),
                       B_s = 0, d = 1,
                       Q = c(1,0,0,1), R = c(10,0,0,10),
                       alpha = 2,  beta = 30, gamma = 20, overwrite = F){
  data_mult = x
  rownames(data_mult) = NULL

  p_length = p_length
  data_allocation = data_mult
  data_allocation[,3] = rep(0, nrow(data_allocation))

  #Initial conditions
  x_hat = matrix(x_hat_s, ncol = 1, nrow = 2)
  P = matrix(sigma_s, ncol = 2 , nrow = 2)
  B_s = tan((B_s*(pi/180)))

  predictions = as.data.frame(matrix(NA, nrow = p_length, ncol = 2))
  error = as.data.frame(matrix(NA, nrow = p_length, ncol = 2))
  search = as.data.frame(matrix(NA, nrow = p_length, ncol = 2))

  #Matrices
  Fk = matrix(c(1,0,0,1), ncol = 2, nrow = 2)
  B = matrix(c(1,1), ncol = 1, nrow = 2)
  d = d
  Q = matrix(Q, ncol = 2, nrow = 2)
  H = matrix(c(1,0,0,1), ncol = 2, nrow = 2)
  R = matrix(R, ncol = 2, nrow = 2)
  I = matrix(c(1,0,0,1), ncol = 2, nrow = 2)

  K = matrix(c(1,0,0,1), ncol = 2, nrow = 2)

  #Value for P expansion if Pk-1 = Pk
  alpha = alpha
  #Value of minimum number of datapoints to constitute a valid crescent
  #Interestingly, higher value = more general, "sharp" trend
  #Lower value = smoother trend, but requires more predictions to reach
  #filament end
  justify = beta
  #Value of max number of rounds of P expansion before algorithm "gives up"
  propagation_limit = gamma

  #Kalman Filter Loop
  for(i in 1:p_length){
    cat("\n")
    #Prediction
    if(i > 2){
      B[1] = predictions[(i-1),1] - predictions[(i-2),1]
      B[2] = predictions[(i-1),2] - predictions[(i-2),2]
    } else {
      B[1] = d
      B[2] = B_s
      }
    #X_hat_k = Fk*x_hat_k-1
    x_hat = (Fk %*% x_hat) + (B %*% d)
    #Pk = Fk*Ftk*P_k-1 + Qk
    P = (Fk %*% t(Fk) %*% P) + Q

    #Kk= Hk*Pk-1/(H*Pk-1*Htk + Rk)
    K = (P %*% t(H)) %*% solve(((H %*% P %*% t(H)) + R))

    #Propagating Crescent
    #collect k z with standard radius P
    z_head = data_mult[which((((data_mult[,1] - x_hat[1,])^2)/(P[1,1]^2)) + (((data_mult[,2] - x_hat[2,])^2)/(P[2,2]^2)) <= 1),]
    index_head = as.numeric(rownames(z_head))
    if(i == 1){z_crescent = z_head} else {z_crescent = z_head[which(!(index_head %in% index_tail)),]}

    track_expansion = 0
    W = P
    #increase P if Pk-1 = Pk
    while(nrow(z_crescent) < justify & track_expansion < propagation_limit){
      W = W * alpha
      z_head =  data_mult[which((((data_mult[,1] - x_hat[1,])^2)/(W[1,1]^2)) + (((data_mult[,2] - x_hat[2,])^2)/(W[2,2]^2)) <= 1),]
      index_head = as.numeric(rownames(z_head))
      if(i == 1){z_crescent = z_head} else {z_crescent = z_head[which(!(index_head %in% index_tail)),]}
      track_expansion = track_expansion + 1
    }
    #If desired, the covariance radius can be overwitten by the value of the search radius. This is false by default.
    if(overwrite == TRUE){P = W}

    search[i,] = c(W[1],W[4])

    index_crescent = as.numeric(row.names(z_crescent))
    data_allocation[index_crescent,3] = i

    z_m = t(t(colMeans(z_crescent)))

    #Update
    #x_hat_k = x_hat_k + K*(zk - H*x_hat_k)
    x_hat = x_hat + (K %*% (z_m - (H %*% x_hat)))
    #Pk = (Ik - Hk*K)*Pk*(Ik - Hk*Kk)t + (Kk*Rk*Kkt)
    P = ((I - (K %*% H)) %*% P %*% t(I - (K %*% H))) + (K %*% R %*% t(K))

    #Save steps
    predictions[i,] = c(x_hat[1], x_hat[2])
    error[i,] = c((P[1] + P[2]),(P[3] + P[4]))

    #k-1 radius of z data
    if(i == 1){z_tail = z_crescent} else {z_tail = rbind(z_tail, z_crescent)}
    index_tail = as.numeric(rownames(z_tail))

    #Print final value (k, number of expanses, elements in z_crescent)
    cat("\r", "[", i, "]","{", track_expansion, "}","(", nrow(z_crescent),")", sep = "")

    if(track_expansion >= propagation_limit){break}
  }

  predictions_out = predictions[complete.cases(predictions),]
  error_out = error[complete.cases(error),]
  search_out = search[complete.cases(search),]
  output = list(x_hat = predictions_out, sigma = error_out, search = search_out, allocated = data_allocation)
  output
}
