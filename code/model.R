if(require(MCMCpack) == FALSE){
    install.packages("MCMCpack", dependencies = TRUE)
}
library(MCMCpack)

gibbs_sampling_beta = function(y, X, beta0, sigma, Sigma0){
    d = ncol(X)
    Sigma = solve(sigma^2*t(X)%*%X + Sigma0)
    beta = Sigma%*%(sigma^2*t(X)%*%y + solve(Sigma0)%*%beta0)
    beta = as.matrix(beta, nrow = d, ncol = 1)
    return(beta)
}

gibbs_sampling_sigma = function(y, X, beta, s0, nu0){
    n = nrow(X)
    s = s0 + t(y - X%*%beta)%*%(y - X%*%beta)
    nu = nu0 + n
    sigma = sqrt(rinvgamma(n = 1, shape = nu/2, scale = s/2))
    return(sigma)
}

regression = function(
    y, x, iter_max = 1000, burn_in = 100,
    d = ncol(x) + 1, beta0 = rep(0, length = d),
    sima0 = 1, Sigma0 = diag(d)
){
    y = as.matrix(y)
    X = as.matrix(cbind(1, x))

    d = ncol(X)
    n = nrow(X)

    beta_array = array(NA, dim = c(d, iter_max + 1))
    beta_array[, 1] = beta0

    sigma_array = array(NA, dim = c(iter_max + 1))
    sigma_array[1] = sigma0

    for(s in 1:iter_max){
        beta_array[, s+1] = gibbs_sampling_beta(y = y, X = X, beta = beta_array[, s], sigma = sigma_array[s], Sigma0)
        sigma_array[s+1] = gibbs_sampling_sigma(y = y, X = X, beta = beta, s0 = s0, nu0 = nu0)
        if(s %% 100 == 0){
            cat("number of iteration is:", s, "\n")
        }
    }

    args_list = list(
        beta = beta_array[, burn_in:iter_max + 1],
        sigma = sigma_array[burn_in:iter_max + 1]
    )
    return(args_list)
    message("the process has finished")
}
