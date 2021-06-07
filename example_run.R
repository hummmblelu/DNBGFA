rm(list = ls())

## loading required functions
source("DNBGFA_var_functions.R")

## loading example data set (Contains objects X, and m_len)
load("example_data.RData")

## model inference 
print("inference model")
model<- model_setup(X, m_len, 20)
update_model(model, iter = 1000, mc.cores = 10, print_topics = T)

## model inference with interpolation (the second time stamp)
print("inference model")
X_remvd<- X[[2]]
X[[2]]<- matrix(NA, nrow(X[[2]]), ncol(X[[2]])) ## make the second matrix NA valued
model<- model_setup(X, m_len, 20)
model<- update_model_interp(model, iter = 1000, inter_t = 2, mc.cores = 10)

## compute expected Z at the time stamp 2
exp_z<- apply(model$var_params$m_eta[,,2], 2, function(x) exp(x) / sum(exp(x)))

## compute mean square error
print(sqrt(mean((X_remvd - exp_z %*% t(model$var_params$exp_w[[2]]))^2)))
