library(parallel)
library(tm)
library(numDeriv)
library(truncnorm)
library(tmvtnorm)
library(msm)
library(mvtnorm)
library(parallel)
library(condMVNorm)
library(MASS)

rbf<- function(zeta, delta) function(t_slice){
  zeta^2 * exp(-as.matrix(dist(t_slice))^2/delta^2)
}

init_var_params<- function(N, K, T, M, m_len){
  var_params<- list(
    m_eta = array(0, dim = c(N, K, T)),
    # m_mu = aperm(array(
    #   unlist(
    #     lapply(1:T, function(t) log(nmf_res@fit@W))
    #   ), dim = c(N, K, T)
    # ), c(1, 2, 3)),
    
    m_mu = array(0, dim = c(N, K, T)),
    S_eta = lapply(1:N, function(n) diag(K * T)),
    
    m_w = lapply(1:T, function(t) matrix(0, sum(m_len[[t]]), K)),
    exp_w = lapply(1:T, function(t) matrix(0, sum(m_len[[t]]), K)),
    S_w = lapply(1:T, function(t) lapply(1:sum(m_len[[t]]), function(d) diag(K))),
    
    # m_eta = aperm(array(
    #   unlist(
    #     lapply(1:T, function(t) log(nmf_res@fit@W))
    #   ), dim = c(N, K, T)
    # ), c(1, 2, 3)),
    # S_mu = diag(K * T),
    
    m_l = matrix(-2, K, T),
    S_l = diag(T),
    
    m_r = rep(.0001, K*T),
    diag_S_r = rep(.0001, K*T),
    S_r<- diag(K*T) * 1,
    r_u = rep(.0001, T),
    m_u_r = matrix(.0001, N, T),
    
    m_s = rep(.0001, M*T),
    diag_S_s = rep(.0001, M*T),
    S_s<- diag(M*T) * 1,
    s_u = rep(0.0001, T),
    m_u_s = matrix(.0001, K, T),
    
    Sigma_eta = diag(K * T) * .1,
    Sigma_alpha = diag(M * T) * .1,
    
    m_alpha = array(-2, dim = c(K, M, T)),
    
    lambda_eta = 0.1,
    lambda_alpha = 0.1
  )
  
  lambda_eta<- var_params$lambda_eta
  m_r<- var_params$m_r
  r_u<- var_params$r_u
  
  D_eta_u<- matrix(0, length(m_r), length(r_u))
  for(t in 1:T){
    decay<- exp(-lambda_eta * (t - 1))
    D_eta_u[(((t-1) * K) + 1) : (t * K) , t]<- matrix(decay, K, 1)
  }
  
  D_eta_eta<- matrix(0, length(m_r), length(m_r))
  for(t_i in 1:T){
    for(t_j in 1:T){
      decay<- exp(-lambda_eta * abs(t_i - t_j))
      D_eta_eta[(((t_i-1) * K) + 1) : (t_i * K) , (((t_j-1) * K) + 1) : (t_j * K)]<- matrix(decay, K, K)
    }
  }
  
  var_params$Kappa_eta<- D_eta_eta * m_r %*% t(m_r)
  var_params$Kappa_eta_u<- D_eta_u * m_r %*% t(r_u)
  
  lambda_alpha<- var_params$lambda_alpha
  m_s<- var_params$m_s
  s_u<- var_params$s_u
  # 
  # 
  D_alpha_u<- matrix(0, length(m_s), length(s_u))
  for(t in 1:T){
    decay<- exp(-lambda_alpha * (t - 1))
    D_alpha_u[(((t-1) * M) + 1) : (t * M) , t]<- matrix(decay, M, 1)
  }
  
  D_alpha_alpha<- matrix(0, length(m_s), length(m_s))
  for(t_i in 1:T){
    for(t_j in 1:T){
      decay<- exp(-lambda_alpha * abs(t_i - t_j))
      D_alpha_alpha[(((t_i-1) * M) + 1) : (t_i * M) , (((t_j-1) * M) + 1) : (t_j * M)]<- matrix(decay, M, M)
    }
  }
  
  var_params$Kappa_alpha<- D_alpha_alpha * m_s %*% t(m_s)
  var_params$Kappa_alpha_u<- D_alpha_u * m_s %*% t(s_u)
  
  var_params
}

update_w<- function(var_params, hyper_params, mc.cores = 10){
  m_eta<- var_params$m_eta
  S_w<- var_params$S_w
  m_alpha<- var_params$m_alpha
  
  m_len<- hyper_params$m_len
  sigma_d<- hyper_params$sigma_d
  N<- hyper_params$N
  M<- hyper_params$M
  T<- hyper_params$T
  
  dis_list<- lapply(1:T, function(t){
    exp_z_t<- apply(m_eta[,,t], 2, function(x) exp(x) / sum(exp(x)))
    Sigma_ard_t<- do.call(cbind, lapply(1:M, function(m) sapply(1:m_len[[t]][m], function(d) exp(2 * m_alpha[,m,t]))))
      
    m_len_t<- m_len[[t]]
    X_t<- X[[t]]
    
    list(exp_z_t = exp_z_t, Sigma_ard_t = Sigma_ard_t, m_len_t = m_len_t, X_t = X_t)
  })
  
  res<- pbmclapply(dis_list, function(l){
    m_len_t<- l$m_len_t
    exp_z_t<- l$exp_z_t
    Sigma_ard_t<- l$Sigma_ard_t
    X_t<- l$X_t
    
    S_W_t<- lapply(1:sum(m_len_t), function(d){
      ginv(t(exp_z_t) %*% exp_z_t + diag(1/Sigma_ard_t[,d]))
    })
    
    m_w_t<- t(sapply(1:sum(m_len_t), function(d) S_W_t[[d]] %*% t(exp_z_t) %*% X_t[,d]))
    
    exp_w_t<- t(sapply(1:sum(m_len_t), function(d){
      out_val<- etruncnorm(a = 0, mean = m_w_t[d,], sd = sigma_d*sqrt(diag(S_W_t[[d]])))
      out_val[is.na(out_val)]<- 0
      out_val
    }))
    
    list(m_w_t = m_w_t, S_W_t = S_W_t, exp_w_t = exp_w_t)
  }, mc.cores = mc.cores)
  
  var_params$S_w<- lapply(res, function(x) x$S_W_t)
  var_params$m_w<- lapply(res, function(x) x$m_w_t)
  var_params$exp_w<- lapply(res, function(x) x$exp_w_t)
  
  var_params
}

update_w_interp<- function(var_params, hyper_params, inter_t, mc.cores = 10){
  m_eta<- var_params$m_eta
  S_w<- var_params$S_w
  m_alpha<- var_params$m_alpha
  
  m_len<- hyper_params$m_len
  sigma_d<- hyper_params$sigma_d
  N<- hyper_params$N
  M<- hyper_params$M
  T<- hyper_params$T
  
  dis_list<- lapply(1:T, function(t){
    exp_z_t<- apply(m_eta[,,t], 2, function(x) exp(x) / sum(exp(x)))
    Sigma_ard_t<- do.call(cbind, lapply(1:M, function(m) sapply(1:m_len[[t]][m], function(d) exp(2 * m_alpha[,m,t]))))
    
    m_len_t<- m_len[[t]]
    X_t<- X[[t]]
    
    list(exp_z_t = exp_z_t, Sigma_ard_t = Sigma_ard_t, m_len_t = m_len_t, X_t = X_t, t_stamp = t)
  })
  
  res<- pbmclapply(dis_list, function(l){
    m_len_t<- l$m_len_t
    exp_z_t<- l$exp_z_t
    Sigma_ard_t<- l$Sigma_ard_t
    X_t<- l$X_t
    t_stamp<- l$t_stamp
    
    if(t_stamp == inter_t){
      ## prior information only
      S_W_t<- lapply(1:sum(m_len_t), function(d) diag(Sigma_ard_t[,d]))
      m_w_t<- matrix(0, sum(m_len_t), nrow(Sigma_ard_t))
    }else{
      S_W_t<- lapply(1:sum(m_len_t), function(d){
        ginv(t(exp_z_t) %*% exp_z_t + diag(1/Sigma_ard_t[,d]))
      })
      
      m_w_t<- t(sapply(1:sum(m_len_t), function(d) S_W_t[[d]] %*% t(exp_z_t) %*% X_t[,d]))
    }
    
    exp_w_t<- t(sapply(1:sum(m_len_t), function(d){
      out_val<- etruncnorm(a = 0, mean = m_w_t[d,], sd = sigma_d*sqrt(diag(S_W_t[[d]])))
      out_val[is.na(out_val)]<- 0
      out_val
    }))
    
    list(m_w_t = m_w_t, S_W_t = S_W_t, exp_w_t = exp_w_t)
  }, mc.cores = mc.cores)
  
  var_params$S_w<- lapply(res, function(x) x$S_W_t)
  var_params$m_w<- lapply(res, function(x) x$m_w_t)
  var_params$exp_w<- lapply(res, function(x) x$exp_w_t)
  
  var_params
}

update_alpha<- function(var_params, hyper_params, mc.cores = 10){
  Kappa_alpha<- var_params$Kappa_alpha
  lambda_alpha<- var_params$lambda_alpha
  exp_w<- var_params$exp_w
  
  m_alpha<- var_params$m_alpha

  K<- hyper_params$K
  T<- hyper_params$T
  M<- hyper_params$M
  epsilon_alpha<- hyper_params$epsilon_alpha
  
  prod_const<- (as.vector(t(matrix(var_params$m_s, M, T))) %*% 
                  t(as.vector(t(matrix(var_params$m_s, M, T))))) + hyper_params$epsilon_alpha * diag(T * M) 
  
  time_decay<- kronecker(
    matrix(1, M, M), sapply(1:T, function(t_i){
      sapply(1:T, function(t_j){
        exp(-var_params$lambda_alpha * (t_i - t_j)^2)
      })
    })
  ) 
  
  Sigma_alpha<- prod_const * time_decay
  Sigma_l<- kronecker(diag(1, M), hyper_params$Kappa_l + diag(T) * hyper_params$epsilon_l)
  
  cost_fun<- function(k) function(alpha){
    alpha<- matrix(alpha, M, T)
    
    -dmvnorm(as.vector(t(alpha)), rep(0, M * T), Sigma_l + Sigma_alpha, log = T) -
      -sum(
        sapply(1:T, function(t){
          ard_sd<- unlist(mapply(rep, x = exp(alpha[,t]), times = m_len[[t]], SIMPLIFY = F))
          msm::dtnorm(exp_w[[t]][,k], mean = 0, sd = ard_sd, lower = 0, upper = Inf, log = T)
        })
      )
  }
  
  dis_list<- lapply(1:K, function(k){
    init<- m_alpha[k,,]
    cost_fun_k<- cost_fun(k)
    
    list(init = init, cost_fun_k = cost_fun_k)
  })
  
  opt_res<- pbmclapply(dis_list, function(l){
    init<- l$init
    cost_fun_k<- l$cost_fun_k
    optim(init, cost_fun_k, method = "SANN", control = list(maxit = 20))$par
  })
  
  for(k in 1:K){
    m_alpha[k,,]<- opt_res[[k]]
  }
  
  var_params$m_alpha<- m_alpha
  var_params
}

update_eta<- function(var_params, hyper_params, mc.cores = 10){
  m_eta<- var_params$m_eta
  m_w<- var_params$m_w
  r_u<- var_params$r_u
  m_u_r<- var_params$m_u_r
  exp_w<- var_params$exp_w
  m_mu<- var_params$m_mu
  lambda_eta<- var_params$lambda_eta
  Kappa_eta<- var_params$Kappa_eta
  Kappa_eta_u<- var_params$Kappa_eta_u
  
  sigma_d<- hyper_params$sigma_d
  epsilon_eta<- hyper_params$epsilon_eta
  K<- hyper_params$K
  N<- hyper_params$N
  T<- hyper_params$T
  
  prod_const<- (as.vector(t(matrix(var_params$m_r, K, T))) %*% 
      t(as.vector(t(matrix(var_params$m_r, K, T)))) + hyper_params$epsilon_eta * diag(T * K)) 
  
  time_decay<- kronecker(
    matrix(1, K, K), sapply(1:T, function(t_i){
      sapply(1:T, function(t_j){
        exp(-var_params$lambda_eta * (t_i - t_j)^2)
      })
    })
  ) 
  
  Sigma_eta<- prod_const * time_decay
  Sigma_mu<- kronecker(diag(1, K), hyper_params$Kappa_mu + diag(T) * hyper_params$epsilon_mu)
  
  cost_fun<- function(n) function(eta_n){
    # eta_n is a K*T matrix
    m_eta_sub<- m_eta
    eta_n<- matrix(eta_n, K, T)
    
    z_mode_t<- lapply(1:T, function(t){
      m_eta_sub[n,,t]<- eta_n[,t]
      eta_t_old<- m_eta[n,,t]
      eta_t_n<- eta_n[,t]
      diff_exp<- exp(eta_t_n) - exp(eta_t_old)
      denom<- apply(m_eta[,,t], 2, function(x) sum(exp(x)))
      t(t(exp(m_eta_sub[,,t])) / (denom + diff_exp))
    })
    
    # a_n<- Kappa_eta_u %*% ginv(Kappa_u_u) %*% as.vector(m_u_r[n,])
    
    # -sum(
    #   sapply(1:T, function(t){
    #     sum(dtnorm(as.vector(X[[t]][n,]), mean = z_mode_t[[t]][n,] %*% t(exp_w[[t]]), sd = sigma_d, lower = 0, log = T))
    #   })) - dmvnorm(x = as.vector(t(eta_n)),
    #            mean = as.vector(t(m_mu[n,,])),
    #            sigma = Sigma_eta, log = T)
    
    -sum(
      sapply(1:T, function(t){
        sum(dtnorm(as.vector(X[[t]][n,]), mean = z_mode_t[[t]][n,] %*% t(exp_w[[t]]), sd = sigma_d, lower = 0, log = T))
      })) - dmvnorm(x = as.vector(eta_n),
                    mean = rep(0, K * T),
                    sigma = (Sigma_eta + Sigma_mu), log = T)
  }
  
  dis_list<- lapply(1:N, function(n){
    init<- m_eta[n,,]
    cost_fun_n<- cost_fun(n)
    
    list(init = init, cost_fun_n = cost_fun_n)
  })
  
  res<- pbmclapply(dis_list, function(l){
    init<- l$init
    cost_fun_n<- l$cost_fun_n
    
    optim(init, cost_fun_n, method = "SANN", control = list(maxit = 10))$par
  }, mc.cores = mc.cores)
  
  for(n in 1:N){
    m_eta[n,,]<- res[[n]]
  }
  
  var_params$m_eta<- m_eta
  var_params
}

update_eta_interp<- function(var_params, hyper_params, inter_t, mc.cores = 10){
  m_eta<- var_params$m_eta
  m_w<- var_params$m_w
  r_u<- var_params$r_u
  m_u_r<- var_params$m_u_r
  exp_w<- var_params$exp_w
  m_mu<- var_params$m_mu
  lambda_eta<- var_params$lambda_eta
  Kappa_eta<- var_params$Kappa_eta
  Kappa_eta_u<- var_params$Kappa_eta_u
  
  sigma_d<- hyper_params$sigma_d
  epsilon_eta<- hyper_params$epsilon_eta
  K<- hyper_params$K
  T<- hyper_params$T
  N<- hyper_params$N
  
  prod_const<- (as.vector(t(matrix(var_params$m_r, K, T))) %*% 
                  t(as.vector(t(matrix(var_params$m_r, K, T)))) + hyper_params$epsilon_eta * diag(T * K)) 
  
  time_decay<- kronecker(
    matrix(1, K, K), sapply(1:T, function(t_i){
      sapply(1:T, function(t_j){
        exp(-var_params$lambda_eta * (t_i - t_j)^2)
      })
    })
  ) 
  
  Sigma_eta<- prod_const * time_decay
  Sigma_mu<- kronecker(diag(1, K), hyper_params$Kappa_mu + diag(T) * hyper_params$epsilon_mu)
  
  cost_fun<- function(n) function(eta_n){
    m_eta_sub<- m_eta
    eta_n<- matrix(eta_n, K, T)
    
    z_mode_t<- lapply(1:T, function(t){
      m_eta_sub[n,,t]<- eta_n[,t]
      eta_t_old<- m_eta[n,,t]
      eta_t_n<- eta_n[,t]
      diff_exp<- exp(eta_t_n) - exp(eta_t_old)
      denom<- apply(m_eta[,,t], 2, function(x) sum(exp(x)))
      t(t(exp(m_eta_sub[,,t])) / (denom + diff_exp))
    })
    
    ## only the prior part influences the cost function
    -sum(
      sapply((1:T)[-inter_t], function(t){
        sum(dtnorm(as.vector(X[[t]][n,]), mean = z_mode_t[[t]][n,] %*% t(exp_w[[t]]), sd = sigma_d, lower = 0, log = T))
      })) - dmvnorm(x = as.vector(eta_n),
                    mean = rep(0, K * T),
                    sigma = (Sigma_eta + Sigma_mu), log = T)
  }
  
  dis_list<- lapply(1:N, function(n){
    init<- m_eta[n,,]
    cost_fun_n<- cost_fun(n)
    
    list(init = init, cost_fun_n = cost_fun_n)
  })
  
  res<- pbmclapply(dis_list, function(l){
    init<- l$init
    cost_fun_n<- l$cost_fun_n
    
    optim(init, cost_fun_n, method = "SANN", control = list(maxit = 10))$par
  }, mc.cores = mc.cores)
  
  for(n in 1:N){
    m_eta[n,,]<- res[[n]]
  }
  
  var_params$m_eta<- m_eta
  var_params
}

update_r<- function(var_params, hyper_params, mc.cores = 10, update_lambda = F){
  m_r<- var_params$m_r
  r_u<- var_params$r_u
  m_eta<- var_params$m_eta
  m_mu<- var_params$m_mu
  lambda_eta<- var_params$lambda_eta
  
  epsilon_eta<- hyper_params$epsilon_eta
  T<- hyper_params$T
  N<- hyper_params$N
  K<- hyper_params$K
  epsilon_r<- hyper_params$epsilon_r
  Kappa_r<- hyper_params$Kappa_r
  a<- hyper_params$a
  b<- hyper_params$b
  Kappa_r<- hyper_params$Kappa_r
  
  Sigma_r_rearranged<- matrix(0, K*T, K*T)
  for(t_i in 1:T){
    for(t_j in 1:T){
      Sigma_r_rearranged[((t_i-1) * K + 1) : (t_i * K), ((t_j-1) * K + 1) : (t_j * K)]<- matrix(Kappa_r[t_i, t_j])
    }
  }
  Sigma_r_rearranged<- Sigma_r_rearranged + diag(nrow(Sigma_r_rearranged)) * epsilon_r
  
  diag_S_r<- var_params$diag_S_r
  S_r<- ginv(ginv(Sigma_r_rearranged) + diag(diag_S_r))
  
  D_eta_eta<- matrix(0, length(m_r), length(m_r))
  for(t_i in 1:T){
    for(t_j in 1:T){
      decay<- exp(-lambda_eta * abs(t_i - t_j))
      D_eta_eta[(((t_i-1) * K) + 1) : (t_i * K) , (((t_j-1) * K) + 1) : (t_j * K)]<- matrix(decay, K, K)
    }
  }
  
  D_eta_u<- matrix(0, length(m_r), length(r_u))
  for(t in 1:T){
    decay<- exp(-lambda_eta * (t - 1))
    D_eta_u[(((t-1) * K) + 1) : (t * K) , t]<- matrix(decay, K, 1)
  }
  
  F_r<- function(m_r, S_r, r_u, lambda_eta){
    psi_0<- t(m_r) %*% m_r + sum(diag(S_r))
    
    Psi_1<- (m_r %*% t(r_u)) * D_eta_u
    
    Psi_2<-  t(D_eta_u * (m_r %*% t(r_u))) %*% (D_eta_u * (m_r %*% t(r_u))) + 
      t(D_eta_u * (diag(S_r) %*% t(r_u))) %*% (D_eta_u * (diag(S_r) %*% t(r_u)))
    
    K_uu<- r_u %*% t(r_u) * diag(sapply(1:T, function(t) exp(-lambda_eta * (t - 1))))
    
    W<- epsilon_eta^(-2) * diag(K * T) - 
      epsilon_eta^(-4) * Psi_1 %*% ginv(epsilon_eta^(-2) * Psi_2 + K_uu) %*% t(Psi_1)
    
    diff_eta_mu<- apply(
      sapply(1:N, function(n) as.vector(m_eta[n,,])), 1, sum
    )
    
    -((-t(diff_eta_mu) %*% W %*% diff_eta_mu) / 2 +
        N * (
          (-K*T) * log(epsilon_eta) + .5 * log(det(K_uu)) - 
            (((K*T)/2) * log(2*pi) + .5 * log(det(epsilon_eta^(-2) * Psi_2 + K_uu)))
        ) +
        (N/(2 * (epsilon_eta^(-2)))) * (sum(diag((ginv(K_uu) %*% Psi_2))) - psi_0))
  }
  
  KL<- function(m_r, S_r){
    .5 * (
      sum(diag(ginv(Sigma_r_rearranged) %*% (S_r + m_r %*% t(m_r)))) + log(det(Sigma_r_rearranged)) - log(det(S_r))
    )
  }
  
  F_r_S<- function(diag_S_root){
    diag_S<- diag_S_root^2
    S<- ginv(ginv(Sigma_r_rearranged) + diag(diag_S))
    F_r(m_r, S, r_u, lambda_eta) + KL(m_r, S)
  } 
  
  ## update S_r
  diag_S_r<- optim(sqrt(diag_S_r), F_r_S, method = "SANN", control = list(maxit = 3))$par^2
  
  S_r<- ginv(ginv(Sigma_r_rearranged) + diag(diag_S_r))
  r_u<- optim(r_u, function(r_u) F_r(m_r, S_r, r_u, lambda_eta) + KL(m_r, S_r), method = "SANN", control = list(maxit = 5))$par
  m_r<- optim(m_r, function(m_r) F_r(m_r, S_r, r_u, lambda_eta) + KL(m_r, S_r), method = "SANN", control = list(maxit = 5))$par
  
  ## update u_r
  Psi_1<- (m_r %*% t(r_u)) * D_eta_u
  
  Psi_2<-  t(D_eta_u * (m_r %*% t(r_u))) %*% (D_eta_u * (m_r %*% t(r_u))) + 
    t(D_eta_u * (diag(S_r) %*% t(r_u))) %*% (D_eta_u * (diag(S_r) %*% t(r_u)))
  
  K_uu<- r_u %*% t(r_u) * diag(sapply(1:T, function(t) exp(-lambda_eta * (t - 1))))
  
  m_u_r<- t(sapply(1:N, function(n){
    K_uu %*% ginv(epsilon_eta^2 * K_uu + Psi_2) %*% t(Psi_1) %*% as.vector(m_eta[n,,])
  }))
  
  if(update_lambda){
    print("update lambda_eta")
    lambda_eta<- optim(lambda_eta, function(lambda_eta){
      K_uu<- r_u %*% t(r_u) * diag(sapply(1:T, function(t) exp(-lambda_eta * (t - 1))))
      
      dis_list<- lapply(seq_len(ncol(m_u_r)), function(i) m_u_r[i,])
      
      log_lkhd_rmv<- sum(
        unlist(
          mclapply(
            dis_list, function(l) dtmvnorm(l, rep(0, T), K_uu, log = T), 
            mc.cores = mc.cores)
        )
      )
      
      F_r(m_r, S_r, r_u, lambda_eta) -
        log_lkhd_rmv - dgamma(lambda_eta, a, b, log = T)
    }, method = "Brent", lower = 0, upper = 10, control = list(maxit = 5))$par
    
    var_params$lambda_eta<- lambda_eta
  }
  
  var_params$r_u<- r_u
  var_params$diag_S_r<- diag_S_r
  var_params$S_r<- S_r
  var_params$m_r<- m_r
  var_params$m_u_r<- m_u_r
  
  Kappa_eta<- (S_r + m_r %*% t(m_r)) * D_eta_eta
  var_params$Kappa_eta<- Kappa_eta
  
  Kappa_eta_u<- D_eta_u * m_r %*% t(r_u)
  var_params$Kappa_eta_u<- Kappa_eta_u
  
  var_params$Sigma_eta<- m_r %*% t(m_r) + S_r
  
  var_params
}

update_s<- function(var_params, hyper_params, mc.cores = 10, update_lambda = F){
  m_s<- var_params$m_s
  S_s<- var_params$S_s
  s_u<- var_params$s_u
  T<- hyper_params$T
  N<- hyper_params$N
  M<- hyper_params$M
  K<- hyper_params$K
  m_alpha<- var_params$m_alpha
  m_l<- var_params$m_l
  diag_S_s<- var_params$diag_S_s
  lambda_alpha<- var_params$lambda_alpha
  
  epsilon_s<- hyper_params$epsilon_s
  epsilon_alpha<- hyper_params$epsilon_alpha
  c<- hyper_params$c
  d<- hyper_params$d
  Kappa_s<- hyper_params$Kappa_s
  
  D_alpha_u<- matrix(0, length(m_s), length(s_u))
  for(t in 1:T){
    decay<- exp(-lambda_alpha * (t - 1))
    D_alpha_u[(((t-1) * M) + 1) : (t * M) , t]<- matrix(decay, M, 1)
  }
  
  D_alpha_alpha<- matrix(0, length(m_s), length(m_s))
  for(t_i in 1:T){
    for(t_j in 1:T){
      decay<- exp(-lambda_alpha * abs(t_i - t_j))
      D_alpha_alpha[(((t_i-1) * M) + 1) : (t_i * M) , (((t_j-1) * M) + 1) : (t_j * M)]<- matrix(decay, M, M)
    }
  }
  
  F_s<- function(m_s, S_s, s_u, lambda_alpha){
    psi_0<- t(m_s) %*% m_s + sum(diag(S_s))
    
    Psi_1<- (m_s %*% t(s_u)) * D_alpha_u
    
    Psi_2<-  t(D_alpha_u * (m_s %*% t(s_u))) %*% (D_alpha_u * (m_s %*% t(s_u))) + 
      t(D_alpha_u * (diag(S_s) %*% t(s_u))) %*% (D_alpha_u * (diag(S_s) %*% t(s_u)))
    
    K_uu<- s_u %*% t(s_u) * diag(sapply(1:T, function(t) exp(-lambda_alpha * (t - 1))))
    
    W<- epsilon_alpha^(-2) * diag(M * T) - 
      epsilon_alpha^(-4) * Psi_1 %*% ginv(epsilon_alpha^(-2) * Psi_2 + K_uu) %*% t(Psi_1)
    
    diff_alpha_zeta<- apply(
      sapply(1:K, function(k) as.vector(m_alpha[k,,])), 1, sum
    )
    
    -((-t(diff_alpha_zeta) %*% W %*% diff_alpha_zeta) / 2 +
        N * (
          (-K*T) * log(epsilon_alpha) + .5 * log(det(K_uu)) - 
            (((K*T)/2) * log(2*pi) + .5 * log(det(epsilon_alpha^(-2) * Psi_2 + K_uu)))
        ) +
        (N/(2 * (epsilon_alpha^(-2)))) * (sum(diag((ginv(K_uu) %*% Psi_2))) - psi_0))
    
  }
  
  Sigma_s_rearranged<- matrix(0, M*T, M*T)
  for(t_i in 1:T){
    for(t_j in 1:T){
      Sigma_s_rearranged[((t_i-1) * M + 1) : (t_i * M), ((t_j-1) * M + 1) : (t_j * M)]<- matrix(Kappa_s[t_i, t_j])
    }
  }
  Sigma_s_rearranged<- Sigma_s_rearranged + diag(nrow(Sigma_s_rearranged)) * epsilon_s
  
  S_s<- ginv(ginv(Sigma_s_rearranged) + diag(diag_S_s))
  
  KL<- function(m_s, S_s){
    .5 * (
      sum(diag(ginv(Sigma_s_rearranged) %*% (S_s + m_s %*% t(m_s)))) + log(det(Sigma_s_rearranged)) - log(det(S_s))
    )
  }
  
  F_s_S<- function(diag_S_root){
    diag_S<- diag_S_root^2
    S<- ginv(ginv(Sigma_s_rearranged) + diag(diag_S))
    F_s(m_s, S, s_u, lambda_alpha) + KL(m_s, S)
  }
  
  diag_S_s<- optim(sqrt(diag_S_s), F_s_S, method = "SANN", control = list(maxit = 5))$par^2
  S_s<- ginv(ginv(Sigma_s_rearranged) + diag(diag_S_s))
  s_u<- optim(s_u, function(s_u) F_s(m_s, S_s, s_u, lambda_alpha), method = "SANN", control = list(maxit = 5))$par
  m_s<- optim(m_s, function(m_s) F_s(m_s, S_s, s_u, lambda_alpha) + KL(m_s, S_s), method = "SANN", control = list(maxit = 5))$par
  
  ## update u_s
  Psi_1<- (m_s %*% t(s_u)) * D_alpha_u
  
  Psi_2<-  t(D_alpha_u * (m_s %*% t(s_u))) %*% (D_alpha_u * (m_s %*% t(s_u))) + 
    t(D_alpha_u * (diag(S_s) %*% t(s_u))) %*% (D_alpha_u * (diag(S_s) %*% t(s_u)))
  
  K_uu<- s_u %*% t(s_u) * diag(sapply(1:T, function(t) exp(-lambda_alpha * (t - 1))))
  m_u_s<- t(
    sapply(1:K, function(k){
      K_uu %*% ginv(epsilon_alpha^2 * K_uu + Psi_2) %*% t(Psi_1) %*% as.vector(m_alpha[k,,])
    })
  )
  
  ## update lambda_alpha
  if(update_lambda){
    print("update lambda_alpha")
    lambda_alpha<- optim(lambda_alpha, function(lambda_alpha){
      K_uu<- s_u %*% t(s_u) * diag(sapply(1:T, function(t) exp(-lambda_alpha * (t - 1))))
      
      dis_list<- lapply(seq_len(ncol(m_u_s)), function(i) m_u_s[i,])
      
      log_lkhd_rmv<- sum(
        unlist(
          mclapply(
            dis_list, function(l) dtmvnorm(l, rep(0, T), K_uu, log = T), 
            mc.cores = mc.cores)
        )
      )
      
      F_s(m_s, S_s, s_u, lambda_alpha) -
        log_lkhd_rmv - dgamma(lambda_alpha, c, d, log = T)
    }, method = "Brent", lower = 0, upper = 5, control = list(maxit = 5))$par
    
    var_params$lambda_alpha<- lambda_alpha
  }
  
  var_params$diag_S_s<- diag_S_s
  var_params$S_s<- S_s
  var_params$s_u<- s_u
  var_params$m_s<- m_s
  var_params$m_u_s<- m_u_s
  
  Kappa_alpha<- (S_s + m_s %*% t(m_s)) * D_alpha_alpha
  var_params$Kappa_alpha<- Kappa_alpha
  
  Kappa_alpha_u<- D_alpha_u * m_s %*% t(s_u)
  var_params$Kappa_alpha_u<- Kappa_alpha_u
  
  var_params$Sigma_alpha<- m_s %*% t(m_s) + S_s
  
  var_params
}

model_setup<- function(X, m_len, K, M = NULL, T = NULL, hyper_params = NULL){
  # X: a list of term-document matrices (tm package) at each time stamp, note that all the matrices should share the same set of vocabularies.
  # m_len: a list of number of documents from each source at each time stamp
  
  if(is.null(hyper_params)){
    print("initializing hyper parameters..")
    
    if(is.null(M)){
      M<- length(m_len[[1]])
    }
    
    if(is.null(T)){
      T<- length(X)
    }
    
    N<- nrow(X[[1]])
    
    hyper_params<- list(
      sigma_d = .001,
      epsilon_eta = .1,
      epsilon_alpha = .1,
      N = N,
      K = K,
      T = T,
      M = M,
      m_len = m_len,
      Kappa_s = rbf(1, 100)(1:T),
      epsilon_s = 1,
      Kappa_r = rbf(1, 100)(1:T),
      epsilon_r = 1,
      Kappa_l = rbf(1, 100)(1:T),
      epsilon_l = 0,
      Kappa_mu = rbf(1, 100)(1:T),
      epsilon_mu = 0,
      a = 1,
      b = 10,
      c = 1,
      d = 10
    )
  }
  
  print("initializing variational parameters..")
  var_params<- init_var_params(N, K, T, M, m_len)
  
  list(X = X, m_len = m_len, var_params = var_params, hyper_params = hyper_params)
}
  
print_topics<- function(model, top_n = 10){
  X<- model$X
  hyper_params<- model$hyper_params
  var_params<- model$var_params
  K<- hyper_params$K
  T<- hyper_params$T
  
  Terms<- attributes(X[[1]])$dimnames[[1]]
  
  for(k in 1:K){
    print(paste("topic", k, sep = ": "))
    for(t in 1:T){
      print(paste("time_stemp", t, sep = ": "))
      print(Terms[order(var_params$m_eta[,k,t], decreasing = T)[1:top_n]])
    }
  }
}
  
update_model<- function(model, iter = 1000, mc.cores = detectCores() - 1, update_lambda = F, print_topics = F){
  X<- model$X
  m_len<- model$m_len
  var_params<- model$var_params
  hyper_params<- model$hyper_params
  
  for(iter_t in 1:iter){
    print(paste("iteration", iter_t, sep = ": "))
    
    print("E-step..")
    print("update eta..")
    var_params<- update_eta(var_params, hyper_params, mc.cores = mc.cores)
    
    print("update alpha..")
    var_params<- update_alpha(var_params, hyper_params, mc.cores = mc.cores)
    
    print("update W..")
    var_params<- update_w(var_params, hyper_params, mc.cores = mc.cores)
    
    print("M-step..")
    print("update r..")
    var_params<- update_r(var_params, hyper_params, mc.cores = mc.cores, update_lambda = update_lambda)
    
    print("update s..")
    var_params<- update_s(var_params, hyper_params, mc.cores = mc.cores, update_lambda = update_lambda)
    
    model$var_params<- var_params
    
    if(print_topics){
      print_topics(model)
    }
  }
  
  return(model)
}

update_model_interp<- function(model, iter = 1000, inter_t, mc.cores = detectCores() - 1, update_lambda = F){
  X<- model$X
  var_params<- model$var_params
  hyper_params<- model$hyper_params
  
  for(iter_t in 1:iter){
    print(paste("iteration", iter_t, sep = ": "))
    
    print("E-step..")
    print("update eta..")
    var_params<- update_eta_interp(var_params, hyper_params, inter_t, mc.cores = mc.cores)
    
    print("update alpha..")
    var_params<- update_alpha(var_params, hyper_params, mc.cores = mc.cores)
    
    print("update W..")
    var_params<- update_w_interp(var_params, hyper_params, inter_t, mc.cores = mc.cores)
    
    print("M-step..")
    print("update r..")
    var_params<- update_r(var_params, hyper_params, mc.cores = mc.cores, update_lambda = update_lambda)
    
    print("update s..")
    var_params<- update_s(var_params, hyper_params, mc.cores = mc.cores, update_lambda = update_lambda)
    
    model$var_params<- var_params
  } 
  
  return(model)
}
