
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
#' logistic regression in terms of the determinant of the Fisher information.
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
    
    S_B  = 0
    S_C1 = 0
    S_C2 = 0
    for (i in 2:(m-1)){
        S_C1 = S_C1 + N[i] * ((x[i]-x[1])^2 + (x[m]-x[i])^2)
        S_C2 = S_C2 + ((n/m+N[1]) * (x[i]-x[1])^2 + (n/m+N[m]) * (x[m]-x[i])^2)
    }
    
    if (m > 3){
        for (i in 2:(m-2)){
            for (j in (i+1):(m-1)){
                S_B = S_B + (N[i]+N[j]+n/m) * (x[j]-x[i])^2
            }
        }
    }
    
    S_A = (N[1]+N[m]+n/2+n/m) * (x[m]-x[1])^2
    
    if (S_B + S_C2 == 0){
        ss = 0
    } else {
        ss = (S_B + S_C2)/(S_A + S_C1)/(m-2)/8
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
#' logistic regression in terms of the determinant of the Fisher information 
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
    
    S_B  = 0
    S_C1 = 0
    S_C2 = 0
    for (i in 2:(m-1)){
        S_C1 = S_C1 + N[i] * ((i-1)^2 + (m-i)^2)
        S_C2 = S_C2 + ((n/m+N[1]) * (i-1)^2 + (n/m+N[m]) * (m-i)^2)
    }
    
    if (m > 3){
        for (i in 2:(m-2)){
            for (j in (i+1):(m-1)){
                S_B = S_B + (N[i]+N[j]+n/m) * (j-i)^2
            }
        }
    }
    
    S_A = (N[1]+N[m]+n/2+n/m) * (m-1)^2
    
    if (S_B + S_C2 == 0){
        ss = 0
    } else {
        ss = (S_B + S_C2)/(S_A + S_C1)/(m-2)/8
    }
    c  = sqrt(ss)
    
    q_range = c(1/2-sqrt(1/4-c), 1/2+sqrt(1/4-c))
    
    q_range
    
}

#' @title Empirical region of response rates that the skip-levels (SL) design is 
#' more efficient than the uniformly-assign (UA) design under two-parameter 
#' logistic regression.
#' 
#' @description Get the grid searched empirical region of response rates of 
#' all doses that the SL design is more efficient than the UA design 
#' under two-parameter logistic regression in terms of the determinant of 
#' the Fisher information.
#'
#' @param x     Vector of values (or transformed values; e.g. log) of doses.
#' @param n     Total sample size at the expansion phase.
#' @param N     Vector of sample size at each dose level.
#' @param by    Space between the grid search points.
#' 
#' @return Sufficient range of response rates.
#' 
#' @examples
#' x  = (1:4)*0.1
#' n  = 28
#' N  = c(3, 3, 3, 6)
#' by = 0.01
#' 
#' sl_q_range_emp(x, n, N, by)
#' 
#' @export
sl_q_range_emp = function(x, n, N, by = 0.01){
    
    m = length(x)
    
    q_1 = seq(by, 1-by, by = by)
    q_m = seq(by, 1-by, by = by)
    
    q_range_emp_rst = expand_grid(q_1, q_m) %>%
        mutate(beta_q1m  = (log(q_m/(1-q_m)) - log(q_1/(1-q_1))) / (x[m] - x[1]),
               alpha_q1m = log(q_m/(1-q_m)) - beta_q1m*x[m],
               rst       = 
                   mapply(det_info_logistic, MoreArgs = list(x = x, N_n = N + c(n/2, rep(0, m-2), n/2)), alpha = alpha_q1m, beta = beta_q1m) -
                   mapply(det_info_logistic, MoreArgs = list(x = x, N_n = N + n/m), alpha = alpha_q1m, beta = beta_q1m))
    
    q_range_emp_rst = q_range_emp_rst %>% 
        filter(rst >= 0)
    
    return(q_range_emp_rst)
    
}

#' @title Independent empirical range of response rates that the skip-levels 
#' (SL) design is more efficient than the uniformly-assign (UA) design 
#' under two-parameter logistic regression.
#' 
#' @description Get the grid searched independent empirical range of 
#' response rates of all doses (symmetric at q = 0.5) that 
#' the SL design is more efficient than 
#' the UA design under two-parameter logistic regression in terms of 
#' the determinant of the Fisher information.
#'
#' @param x     Vector of values (or transformed values; e.g. log) of doses.
#' @param n     Total sample size at the expansion phase.
#' @param N     Vector of sample size at each dose level.
#' @param by    Space between the grid search points.
#' 
#' @return Sufficient range of response rates.
#' 
#' @examples
#' x  = (1:4)*0.1
#' n  = 28
#' N  = c(3, 3, 3, 6)
#' by = 0.01
#' 
#' sl_q_range_emp_ind(x, n, N, by)
#' 
#' @export
sl_q_range_emp_ind = function(x, n, N, by = 0.01){
    
    m = length(x)
    
    q = seq(by, 1-by, by = by)
    
    q_range_emp_rst = expand_grid(q) %>%
        mutate(beta_q1m  = (log((1-q)/q) - log(q/(1-q))) / (x[m] - x[1]),
               alpha_q1m = log((1-q)/q) - beta_q1m*x[m],
               rst       = 
                   mapply(det_info_logistic, MoreArgs = list(x = x, N_n = N + c(n/2, rep(0, m-2), n/2)), alpha = alpha_q1m, beta = beta_q1m) -
                   mapply(det_info_logistic, MoreArgs = list(x = x, N_n = N + n/m), alpha = alpha_q1m, beta = beta_q1m))
    
    q_range_emp_rst = q_range_emp_rst %>% 
        filter(rst >= 0)
    
    q_range_emp = range(q_range_emp_rst$q)
    
    return(q_range_emp)
    
}

#' @title Region plot for sufficient region of the theorem, empirical region
#' and independent empirical region that the SL design is more efficient than 
#' the UA design under two-parameter logistic regression in terms of 
#' the determinant of the Fisher information. 
#' 
#' @description Compare the regions of the theorem, empirical grid search and 
#' the empirical grid search with similar structure as that of the theorem that 
#' the SL design is more efficient than 
#' the UA design under two-parameter logistic regression in terms of 
#' the determinant of the Fisher information. 
#'
#' @param x     Vector of values (or transformed values; e.g. log) of doses.
#' @param n     Total sample size at the expansion phase.
#' @param N     Vector of sample size at each dose level.
#' @param by    Space between the grid search points.
#' 
#' @return Sufficient range of response rates.
#' 
#' @examples
#' x  = (1:4)*0.1
#' n  = 28
#' N  = c(3, 3, 3, 6)
#' by = 0.01
#' 
#' sl_q_range_plot(x, n, N, by)
#' 
#' @export
sl_q_range_plot = function(x, n, N, by = 0.01){
    
    q_range             = sl_q_range(x, n, N)
    q_range_emp_rst     = sl_q_range_emp(x, n, N, by)
    q_range_emp_ind_rst = sl_q_range_emp_ind(x, n, by)
    
    q_range_emp = rbind(q_range_emp_rst %>% 
                            group_by(q_1) %>% 
                            filter(q_m == max(q_m)),
                        q_range_emp_rst %>% 
                            group_by(q_1) %>% 
                            filter(q_m == min(q_m)) %>% 
                            arrange(desc(q_1)))
    
    q_range_emp_ind = data.frame(
        x = c(q_range_emp_ind_rst[1], q_range_emp_ind_rst[1], q_range_emp_ind_rst[2], q_range_emp_ind_rst[2]),
        y = c(q_range_emp_ind_rst[1], q_range_emp_ind_rst[2], q_range_emp_ind_rst[2], q_range_emp_ind_rst[1])
    )
    
    q_range_the = data.frame(
        x = c(q_range[1], q_range[1], q_range[2], q_range[2]),
        y = c(q_range[1], q_range[2], q_range[2], q_range[1])
    )
    
    plot_q_range = ggplot() +
        geom_polygon(
            aes(x = q_1, y = q_m, fill = "Empirical"),
            data = q_range_emp,
            alpha = 0.3,
            color = "green3"
        ) + geom_polygon(
            aes(x = x, y = y, fill = "Empirical (Independent)"),
            data = q_range_emp_ind,
            alpha = 0.3,
            color = "orange3"
        ) +
        geom_polygon(
            aes(x = x, y = y, fill = "Theorem"),
            data = q_range_the,
            alpha = 0.7,
            color = "blue"
        ) +
        scale_fill_manual(
            name = NULL,
            values = c("Empirical" = "#00BA38", 
                       "Empirical (Independent)" = "orange",
                       "Theorem" = "skyblue")
        ) +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
        scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
        theme_bw() +
        labs(
            x = expression(tilde(q)[1]),
            y = expression(tilde(q)[m])
        ) +
        theme(legend.position = "bottom",
              text = element_text(size = 15))
    
    return(plot_q_range)
    
}

#' @title Determinant of the Fisher information matrix of the 
#' two-parameter logistic regression.
#' 
#' @description Calculate the determinant of the Fisher information matrix of 
#' the two-parameter logistic regression
#'
#' @param x     Vector of values (or transformed values; e.g. log) of doses.
#' @param N_n   Total sample size at both the escalation and the expansion phases.
#' @param alpha Intercept of the logistic regression
#' @param beta  Slope of the logistic regression
#' 
#' @return Sufficient range of response rates.
#' 
#' @examples
#' x     = (1:4)*0.1
#' N_n   = c(6, 3, 3, 6)
#' alpha = -1.4
#' beta  = 5.6
#' 
#' det_info_logistic(x, N_n, alpha, beta)
#' 
#' @export
det_info_logistic = function(x, N_n, alpha, beta){
    
    q = 1/(1+exp(-alpha-beta*x))
    
    v = N_n*q*(1-q)
    
    Info      = matrix(NA, 2, 2)
    Info[1,1] = sum(v)
    Info[1,2] = sum(v*x)
    Info[2,1] = Info[1,2]
    Info[2,2] = sum(v*x^2)
    
    det(Info)
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
