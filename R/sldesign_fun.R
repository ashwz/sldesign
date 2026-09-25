
#' @title Get the two-parameter logistic regression given the ranges 
#' of the response rates.
#' 
#' @description Get the response rates and the values of the parameters 
#' under logistic regression given the ranges of the response rates.
#' 
#' @param dose     Vector of values of doses.
#' @param q_range  Matrix of ranges of the response rates q (a matrix with 
#'                 2 columns, each row represents the range of q).
#' 
#' @return The list of the response rates and the values of the parameters.
#' 
#' @examples
#' dose    = (1:4)*0.1
#' q_range = rbind(c(0.2, 0.6), c(0.3, 0.7), c(0.4, 0.8))
#' sl_scenario_range(dose, q_range)
#' 
#' @export
sl_scenario_range = function(dose    = (1:4)*0.1, 
                             q_range = rbind(c(0.2, 0.6), c(0.3, 0.7))){
    
    if (sum(diff(dose) <= 0) > 0){
        stop("dose should be in increasing order and no tie.")
    }
    
    q_all   = NULL
    par_all = NULL
    for (i in 1:nrow(q_range)){
        logit_range = log(q_range[i, ]/(1-q_range[i, ]))
        x_range     = range(dose)
        
        beta  = diff(logit_range) / diff(x_range)
        alpha = logit_range[1] - beta*x_range[1]
        
        q_all   = rbind(q_all, 1/(1+exp(-alpha-beta*dose)))
        par_all = rbind(par_all, data.frame(alpha = alpha, beta = beta))
        
    }
    
    return(list(q_all = q_all, par_all = par_all))
}


#' @title Generate dose-response data.
#' 
#' @description Generate dose-response data given sample size and 
#' response rate per dose.
#' 
#' @param x    Vector of values (or transformed values; e.g. log) of doses.
#' @param q    Vector of response rates.
#' @param N_n  Vector of cumulative sample sizes at both escalation and 
#'             expansion phases.
#' 
#' @return A data frame with three columns: x (values or transformed 
#' values of doses), 
#' n (sample sizes), y (number of responders).
#' 
#' @examples
#' set.seed(1000)
#' 
#' x   = (1:4)*0.1
#' q   = c(0.20, 0.31 ,0.45, 0.60)
#' N_n = c(16, 10, 10, 16)
#' 
#' sl_data_sim(x, q, N_n)
#' 
#' @export
sl_data_sim = function(x, q, N_n){
    
    y = rbinom(length(N_n), N_n, q)
    data.frame(x = x, N_n = N_n, y = y)
    
}


#' @title Sufficient range of response rates that the skip-levels (SL) design is 
#' more efficient than the uniformly-assign (UA) design under two-parameter 
#' logistic regression.
#' 
#' @description Get the sufficient range of response rates of all doses that 
#' the SL design is more efficient than the UA design under two-parameter 
#' logistic regression in terms of the determinant of the fisher information.
#' 
#' @param x    Vector of values (or transformed values; e.g. log) of doses.
#' @param n    Total sample size at the expansion phase.
#' @param N    Vector of sample size at each dose level.
#' 
#' @return Sufficient range of response rates.
#' 
#' @examples
#' x = (1:4)*0.1
#' n = 28
#' N = c(3, 3, 3, 6)
#' 
#' sl_q_range(x, n, N)
#' 
#' @export
sl_q_range = function(x = (1:4)*0.1, n = 28, N = rep(0, length(x))){
    
    m  = length(N)
    
    temp   = 0
    temp_2 = 0
    temp_4 = 0
    for (i in 2:(m-1)){
        temp_2 = temp_2 + N[i]*(1/2-1/m) * ((x[1]-x[i])^2 + (x[i]-x[m])^2)
        temp_4 = temp_4 + ((n/m+N[1]) * (x[1]-x[i])^2 + (n/m+N[m]) * (x[i]-x[m])^2) / 16/m
    }
    
    if (m > 3){
        for (i in 2:(m-2)){
            for (j in (i+1):(m-1)){
                temp  = temp  + (N[i]+N[j]+n/m) * (x[i]-x[j])^2 /16/m
            }
        }
    }
    
    temp_3 = (1/2-1/m) * (N[1]+N[m]+n/2+n/m) * (x[1]-x[m])^2
    
    if (temp + temp_4 == 0){
        ss = 0
    } else {
        ss = (temp + temp_4)/(temp_3 + temp_2)
    }
    c  = sqrt(ss)
    
    q_range = c(1/2-sqrt(1/4-c), 1/2+sqrt(1/4-c))
    
    q_range
    
}

#' @title Sufficient range of response rates that the skip-levels (SL) design is 
#' more efficient than the uniformly-assign (UA) design under two-parameter 
#' logistic regression when the values of doses (or transformed doses) are 
#' equally spaced.
#' 
#' @description Get the sufficient range of response rates of all doses that 
#' the SL design is more efficient than the UA design under two-parameter 
#' logistic regression in terms of the determinant of the fisher information 
#' when the values of doses (or transformed doses) are equally spaced.
#'
#' @param n    Total sample size at the expansion phase.
#' @param N    Vector of sample size at each dose level.
#' 
#' @return Sufficient range of response rates.
#' 
#' @examples
#' n    = 28
#' N    = c(3, 3, 3, 6)
#' sl_q_range_es(n, N)
#' 
#' @export
sl_q_range_es = function(n = 28, N = rep(0, 4)){
    
    m  = length(N)
    
    temp   = 0
    temp_2 = 0
    temp_4 = 0
    for (i in 2:(m-1)){
        temp_2 = temp_2 + N[i]*(1/2-1/m) * ((i-1)^2 + (i-m)^2)
        temp_4 = temp_4 + ((n/m+N[1]) * (i-1)^2 + (n/m+N[m]) * (i-m)^2) / 16/m
    }
    
    if (m > 3){
        for (i in 2:(m-2)){
            for (j in (i+1):(m-1)){
                temp  = temp  + (N[i]+N[j]+n/m) * (i-j)^2 /16/m
            }
        }
    }
    
    temp_3 = (1/2-1/m) * (N[1]+N[m]+n/2+n/m) * (m-1)^2
    
    if (temp + temp_4 == 0){
        ss = 0
    } else {
        ss = (temp + temp_4)/(temp_3 + temp_2)
    }
    c  = sqrt(ss)
    
    q_range = c(1/2-sqrt(1/4-c), 1/2+sqrt(1/4-c))
    
    q_range
    
}

#' @title Posterior mode of parameters of Bayesian two-parameter 
#' logistic regression.
#' 
#' @description Calculate posterior mode of parameters of Bayesian 
#' two-parameter logistic regression.
#'
#' @param data        Data from function sl_data_sim().
#' @param mu_prior    Mean of the prior normal distribution.
#' @param sigma_prior Standard error of the prior normal distribution.
#' 
#' @return Posterior mode of parameters.
#' 
#' @examples
#' x    = (1:4)*0.1
#' q    = c(0.20, 0.31 ,0.45, 0.60)
#' N_n  = c(16, 10, 10, 16)
#' data = sl_data_sim(x, q, N_n)
#' 
#' mu_prior    = c(0, 0)
#' sigma_prior = c(10, 10)
#' sl_post_mode_logistic(data, mu_prior, sigma_prior)
#' 
#' @export
sl_post_mode_logistic = function(data, 
                                 mu_prior    = c(0, 0), 
                                 sigma_prior = c(100, 100)){
    
    post_dist_logistic = function(data, par, mu_prior, sigma_prior){
        
        q = 1/(1+exp(- par[1] - par[2]*data$x))
        
        like_binom = dbinom(data$y, data$N_n, q)
        like_norm  = dnorm(par, mu_prior, sigma_prior)
        
        - ( sum(log(like_binom)) + sum(log(like_norm)) )
        
    }
    
    rst = optim(c(0, 0), post_dist_logistic, 
                data        = data, 
                mu_prior    = mu_prior, 
                sigma_prior = sigma_prior)
    
    rst$par
    
}
