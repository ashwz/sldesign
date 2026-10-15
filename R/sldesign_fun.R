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
    q_range_emp_ind_rst = sl_q_range_emp_ind(x, n, N, by)
    
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





#' @title Fit the Bayesian model for all combinations of dose-expansion responses
#' (SL or UA design).
#'
#' @description In the SL design, the dose-expansion participants are split equally
#' between the lowest (DL1) and the highest (DLm) dose levels; in the UA design
#' they are split equally across all m dose levels. For every possible
#' combination of expansion responses, the Bayesian two-parameter logistic
#' regression model is fitted to the combined dose-escalation and
#' dose-expansion data, and the posterior draws of the ratios q_i / q_m,
#' i = 1, ..., m - 1, are kept. Use \code{sl_noninf_table()} to create the
#' non-inferiority table from the fit for any delta and tau.
#'
#' @param N            Vector of sample sizes per dose level at the
#'                     dose-escalation phase.
#' @param Y            Vector of numbers of responses per dose level at the
#'                     dose-escalation phase.
#' @param n            Total sample size of the dose-expansion phase, split
#'                     equally across the dose levels in the expansion
#'                     (SL: n / 2 at DL1 and at DLm, n should be even;
#'                     UA: n / m at each dose level, n should be a multiple of m).
#' @param method       "SL" or "UA".
#' @param x            Vector of values (or transformed values; e.g. log) of
#'                     doses, in increasing order.
#' @param mu_prior     Vector of the normal mean for two-parameter logistic
#'                     regression.
#' @param sigma_prior  Vector of the normal standard error for two-parameter
#'                     logistic regression.
#' @param DL           Labels of the dose levels used in the results, e.g.
#'                     c("48mg", "96mg", "200mg"). Default "DL1", ..., "DLm".
#'
#' @details The progress ("combination k / total") is shown in the console, and
#' also in a shiny progress bar when called inside \code{shiny::withProgress()}.
#'
#' @return A list: grid (the combinations of expansion responses, one row each),
#' ratio_draws (array: combination x MCMC draw x i), rates (observed response
#' rate per dose level, escalation + expansion, for each combination), dl_exp
#' (dose levels in the expansion), n_per (expansion sample size per dose level in
#' the expansion), m (number of dose levels), n_combo (number of combination), 
#' method and DL (labels of the dose levels).
#'
#' @examples
#' set.seed(1000)
#' 
#' N = rep(10, 4)       # dose escalation: sample size per dose level
#' Y = c(3, 4, 5, 6)    # dose escalation: responses per dose level
#' n = 6                # dose expansion: total sample size
#' 
#' fit = sl_logistic_expan_comb(N, Y, n, method = "SL")
#' 
#' @export
sl_logistic_expan_comb = function(N, Y, n,
                                  method      = c("SL", "UA"),
                                  x           = (1:length(N))*0.1,
                                  mu_prior    = c(0, 0),
                                  sigma_prior = c(100, 100),
                                  DL          = paste0("DL", seq_along(N))){
    
    method = match.arg(method)
    
    m = length(N)
    if (m < 2)                stop("At least 2 dose levels are required.")
    if (length(Y) != m)       stop("Y should have the same length as N.")
    if (length(x) != m)       stop("x should have the same length as N.")
    if (any(Y > N))           stop("Y cannot exceed N.")
    if (any(diff(x) <= 0))    stop("x should be in increasing order and no tie.")
    DL = as.character(DL)
    if (length(DL) != m || anyDuplicated(DL) > 0) stop("DL should have the same length as N and no tie.")
    
    ## dose levels in the expansion: SL = lowest and highest; UA = all
    dl_exp = if (method == "SL") c(1, m) else seq_len(m)
    k_exp  = length(dl_exp)
    if (n < k_exp || n %% k_exp != 0){
        stop(sprintf("n should be a multiple of %d (the number of dose levels in the %s expansion).",
                     k_exp, method))
    }
    
    n_per = n / k_exp                              # expansion sample size per dose level
    grid  = expand.grid(rep(list(0:n_per), k_exp)) # all combinations of expansion responses
    n_add = rep(0, m)
    n_add[dl_exp] = n_per
    n_tot = N + n_add
    
    ## ---- fit the model for every combination of expansion responses ----
    n_combo     = nrow(grid)
    ratio_draws = NULL     # array: combination x MCMC draw x i (i = 1, ..., m - 1)
    rates       = matrix(NA_real_, nrow = n_combo, ncol = m)
    
    ## progress: "k / n_combo" in the console; inside a shiny withProgress() the
    ## bar and its detail text are updated as well
    in_shiny = requireNamespace("shiny", quietly = TRUE) &&
        !is.null(shiny::getDefaultReactiveDomain())
    
    for (k in seq_len(n_combo)){
        y_exp         = rep(0, m)
        y_exp[dl_exp] = as.numeric(grid[k, ])
        y             = Y + y_exp
        
        fit = sl_bayes_logistic_model(data.frame(x = x, N_n = n_tot, y = y),
                                      mu_prior    = mu_prior,
                                      sigma_prior = sigma_prior)
        post_q_all = fit$post_q_all      # posterior draws x dose levels
        
        if (is.null(ratio_draws)){
            ratio_draws = array(NA_real_, dim = c(n_combo, nrow(post_q_all), m - 1))
        }
        ratio_draws[k, , ] = post_q_all[, 1:(m - 1), drop = FALSE] / post_q_all[, m]
        rates[k, ]         = y / n_tot
        
        cat(sprintf("\rFitting the model: combination %d / %d", k, n_combo))
        utils::flush.console()
        if (in_shiny){
            try(shiny::incProgress(1 / n_combo,
                                   detail = sprintf("%d / %d", k, n_combo)),
                silent = TRUE)
        }
    }
    cat("\n")
    
    list(grid = grid, ratio_draws = ratio_draws, rates = rates,
         dl_exp = dl_exp, n_per = n_per, m = m, n_combo = n_combo,
         method = method, DL = DL)
}


#' @title Non-inferiority decision table for all combinations of dose-expansion responses.
#'
#' @description Create the non-inferiority decision table from the model fitted by
#' \code{sl_logistic_expan_comb()}. The non-inferiority decision table has 
#' one row per combination of
#' expansion responses, with the combination (e.g. "(5, 7)"), the observed
#' response rate per dose level (escalation + expansion), the posterior
#' probabilities P(q_i / q_m > delta), i = 1, ..., m - 1, and the lowest DL in
#' the admissible region (AR): the lowest dose level i (i < m) with
#' P(q_i / q_m > delta) > tau; if no lower dose level meets the threshold, the
#' highest dose level m. No model is fitted here, so the table can be created
#' quickly for different delta and tau.
#'
#' @param fit    Output of \code{sl_logistic_expan_comb()}.
#' @param delta  Non-inferiority threshold (ratio of response rates).
#' @param tau    Posterior probability threshold.
#'
#' @return A data frame (numbers are not rounded). The last column,
#' "Lowest DL in AR", is a factor with levels DL1, ..., DLm. The attribute
#' "fit" holds the model fit, used by \code{sl_noninf_plot()} and
#' \code{sl_noninf_prob()}.
#'
#' @examples
#' set.seed(1000)
#' 
#' N = rep(10, 4)       # dose escalation: sample size per dose level
#' Y = c(3, 4, 5, 6)    # dose escalation: responses per dose level
#' n = 6                # dose expansion: total sample size
#' 
#' fit = sl_logistic_expan_comb(N, Y, n, method = "SL")
#' 
#' sl_noninf_table(fit, delta = 0.85, tau = 0.5)
#' 
#' @export
sl_noninf_table = function(fit, delta = 0.9, tau = 0.8){
    
    m = fit$m
    
    ## ---- posterior probabilities P(q_i / q_m > delta): combinations x (m - 1) ----
    probs = apply(fit$ratio_draws > delta, c(1, 3), mean)
    
    ## ---- lower bound in AR ----
    sel = apply(probs, 1, function(q){
        ok = which(q > tau)
        if (length(ok) > 0) min(ok) else m
    })
    
    ## ---- non-inferiority decision table ----
    ## expansion responses in ONE column as a tuple, e.g. "(5, 7)"
    resp = data.frame(
        tmp = apply(as.matrix(fit$grid), 1, function(v)
            paste0("(", paste(sprintf("%d", as.integer(round(v))), collapse = ", "), ")")),
        check.names = FALSE)
    colnames(resp) = paste0("Expansion responses (n=", fit$n_per, " each)")
    
    ## dose level labels: "DL1", ..., "DLm", or the labels given in DL
    lab = if (is.null(fit$DL)) paste0("DL", seq_len(m)) else fit$DL
    
    ## observed response rate per DL (escalation + expansion)
    rt = as.data.frame(fit$rates)
    colnames(rt) = paste0("Observed rate ", lab, " (esc + exp)")
    
    pm = as.data.frame(probs)
    colnames(pm) = paste0("P(", lab[1:(m - 1)], "/", lab[m], " > \u03b4)")
    
    sel_df = data.frame(`Lowest DL in AR` = factor(lab[sel], levels = lab),
                        check.names = FALSE)
    
    out = cbind(resp, rt, pm, sel_df)
    rownames(out) = NULL
    
    attr(out, "fit") = fit
    
    out
}


#' @title Non-inferiority decision plot (SL design).
#'
#' @description Draw the non-inferiority decision plot: each block
#' is one combination of dose-expansion responses (x: responses at the lowest
#' dose level, y: responses at the highest dose level), labeled and colored by
#' the lowest DL in AR (no legend; the plot title is "Lowest dose that is non-inferior to the high dose").
#' Grey block: no lower dose level is non-inferior to the highest dose level
#' (the highest dose level is selected). Other colors: the lowest dose level
#' that is non-inferior to the highest dose level.
#'
#' @param tab  Non-inferiority table from \code{sl_noninf_table()} with method = "SL".
#'
#' @return A ggplot object.
#'
#' @examples
#' set.seed(1000)
#' 
#' fit = sl_logistic_expan_comb(N = rep(10, 4), Y = c(3, 4, 5, 6), n = 6)
#' tab = sl_noninf_table(fit, delta = 0.85, tau = 0.5)
#' 
#' sl_noninf_plot(tab)
#'
#' @export
sl_noninf_plot = function(tab){
    
    fit = attr(tab, "fit")
    if (is.null(fit) || length(fit$dl_exp) != 2 ||
        (!is.null(fit$method) && fit$method != "SL")){
        stop("The non-inferiority table plot is only available for the SL design (use the output of sl_noninf_table()).")
    }
    
    grid = as.data.frame(fit$grid)
    lev  = if (is.null(fit$DL)) paste0("DL", seq_len(fit$m)) else fit$DL
    
    df = data.frame(x   = grid[[1]],
                    y   = grid[[2]],
                    sel = tab[["Lowest DL in AR"]])
    
    cols = setNames(c(scales::hue_pal()(fit$m)), lev)
    
    ## size of the labels shrinks when there are many blocks
    txt_size = max(2.5, min(5, 45 / (fit$n_per + 1)))
    
    ggplot(df, aes(x = x, y = y, fill = sel)) +
        geom_tile(color = "white", linewidth = 0.3) +
        # geom_text(aes(label = sel), size = txt_size) +
        scale_fill_manual(values = cols, drop = FALSE, name = "Lowest DL in AR") +
        scale_x_continuous(breaks = 0:fit$n_per, expand = c(0.01, 0.01)) +
        scale_y_continuous(breaks = 0:fit$n_per, expand = c(0.01, 0.01)) +
        labs(title = paste0("Lowest dose that is non-inferior to the high dose (", lev[fit$dl_exp[2]], ")"),
             x = paste0("Expansion responses at ", lev[fit$dl_exp[1]], " (n=", fit$n_per, ")"),
             y = paste0("Expansion responses at ", lev[fit$dl_exp[2]], " (n=", fit$n_per, ")")) +
        theme_bw(base_size = 16) +
        theme(panel.grid      = element_blank(),
              axis.title      = element_text(),
              plot.title      = element_text(hjust = 0.5),
              legend.position = "bottom",
              legend.title = element_blank())
}


#' @title Probability of being the lowest non-interior dose.
#'
#' @description Given the non-inferiority table and the assumed true response rates
#' at the dose levels in the expansion, the probability of observing each
#' combination of expansion responses is the product of the binomial
#' probabilities over those dose levels. The probability of each dose level to
#' be the lower bound of the admissible region (AR) is the sum of these
#' probabilities over the combinations whose lower bound in the AR is that dose level.
#'
#' @param tab  Non-inferiority table from \code{sl_noninf_table()}.
#' @param q    Assumed true response rates at the dose levels in the expansion
#'             (SL: 2 numbers, the lowest then the highest dose level;
#'             UA: m numbers, one per dose level).
#'
#' @return A one-row data frame with columns DL1, ..., DLm: the probability of
#' each dose level to be the lowest DL in the AR (they sum to 1).
#'
#' @examples
#' set.seed(1000)
#' 
#' fit = sl_logistic_expan_comb(N = rep(10, 4), Y = c(3, 4, 5, 6), n = 6)
#' tab = sl_noninf_table(fit, delta = 0.85, tau = 0.5)
#' 
#' sl_noninf_prob(tab, q = c(0.4, 0.6))
#'
#' @export
sl_noninf_prob = function(tab, q){
    
    fit = attr(tab, "fit")
    if (is.null(fit)){
        stop("tab should be the output of sl_noninf_table().")
    }
    
    m = fit$m
    k = length(fit$dl_exp)
    
    ## SL design: allow the response rates of all m dose levels
    if (length(q) == m && m != k) q = q[fit$dl_exp]
    
    if (length(q) != k || any(is.na(q)) || any(q < 0 | q > 1)){
        stop(sprintf("q should be %d response rates between 0 and 1 (one per dose level in the expansion).", k))
    }
    
    ## probability of each combination of expansion responses
    prob = apply(as.matrix(fit$grid), 1, function(y) prod(dbinom(y, fit$n_per, q)))
    
    ## sum over the combinations with the same lowest DL in AR
    sel = tab[["Lowest DL in AR"]]
    out = sapply(levels(sel), function(l) sum(prob[sel == l]))
    
    out = as.data.frame(t(out))
    rownames(out) = NULL
    
    out
}

