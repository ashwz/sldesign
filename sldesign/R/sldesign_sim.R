
#' @title Parameter values of response-toxicity multinomial distribution.
#' 
#' @description Get the parameter values of response-toxicity 
#' multinomial distribution given the toxicity rate, the response rate and 
#' the odds ratio.
#' 
#' @param gamma Odds ratio of the response rates in toxic group versus 
#'              that in non-toxic group.
#' @param p     Vector of the toxicity rates.
#' @param q     Vector of the response rates.
#' 
#' @return P: the parameters p_11, p_10, p_01, p_00 (11: a toxicity 
#' and a response; 10: a toxicity and no response; 01: no toxicity 
#' and a response; 00: no toxicity and no response); rho: the correlation 
#' of the binary toxicity indicator and the binary response indicator.
#' 
#' @examples
#' gamma = 0.5
#' p     = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q     = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' 
#' sl_multinom_par(gamma, p, q)
#' 
#' @export
sl_multinom_par = function(gamma, p, q){
    
    num_DL = length(p)
    
    a = 1 - gamma; 
    b = 1 - (p+q)*(1-gamma); 
    c = -gamma*p*q
    
    if (gamma==1) {
        p_11 = p*q
    } else {
        p_11 = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
    }
    p_10 = p - p_11
    p_01 = q - p_11
    p_00 = 1 - p - q + p_11
    rho  = (p_11-p*q)/sqrt(p*q*(1-p)*(1-q))
    
    P           = cbind(p_11, p_10, p_01, p_00)
    rownames(P) = lapply(1:num_DL, function(x){paste("DL",x,sep="")})
    colnames(P) = c("p_11","p_10","p_01","p_00")
    
    return(list(P = P, rho = rho))
}

#' @title Scenarios of dose-response curves under logistic regression.
#' 
#' @description Get the dose-response curves under logistic regression, where 
#' the curve is determined by the dose-response data point of the maximum 
#' tolerated dose (MTD) and the MTD - (k+1) with response rate equals 
#' delta*q_MTD (k: number of dose levels (DLs) to be skipped; q_MTD: 
#' respnose rate at MTD; delta: a large constant in (0, 1)).
#' 
#' @param max_skip  Maximum DLs to be skipped. The results will show the case 
#'                  when skipping 0 to max_skip DLs.
#' @param dose      Vector of values of doses.
#' @param q_MTD     Response rate at the MTD.
#' @param MTD       The DL of the MTD.
#' @param delta     Relative efficiency of the response rates of the 
#'                  non-inferiority dose and the MTD. It is a 
#'                  constant in (0, 1).
#' 
#' @return q_all: response rates matrix, each row represents the response 
#' rates scenario for each skip value; par_all: parameters matrix, each row 
#' represents the values of parameters per each skip value.
#' 
#' @examples
#' sl_scenario_dede(max_skip = 2, 
#'                  dose     = (1:6)*0.1, 
#'                  q_MTD    = 0.6, 
#'                  MTD      = 5, 
#'                  delta    = 0.9)
#' 
#' @export
sl_scenario_dede = function(max_skip = 2, 
                            dose     = 1:6, 
                            q_MTD    = 0.6, 
                            MTD      = 5, 
                            delta    = 0.9){
    
    if (sum(diff(dose) <= 0) > 0){
        stop("dose should be in increasing order and no tie.")
    }
    
    q_range     = c(q_MTD*delta, q_MTD)
    logit_range = log(q_range/(1-q_range))
    
    q_all   = NULL
    par_all = NULL
    for (sc_skip in 0:max_skip){
        
        x_range = dose[c(MTD - sc_skip - 1, MTD)]
        
        beta  = diff(logit_range) / diff(x_range)
        alpha = logit_range[1] - beta*x_range[1]
        
        q_all   = rbind(q_all, round(1/(1+exp(-alpha-beta*dose)),3))
        par_all = rbind(par_all, data.frame(alpha = alpha, beta = beta))
    }
    
    return(list(q_all = q_all, par_all = par_all))
}

#' @title Simulate dose-escalation phase data using BOIN design.
#' 
#' @description Simulate dose-escalation phase data using BOIN design under 
#' multinomial distribution assumption for response and toxicity.
#' 
#' @param p      Vector of the toxicity rates.
#' @param q      Vector of the response rates.
#' @param gamma  Odds ratio of the response rates in toxic group versus that 
#'               in non-toxic group.
#' @param dose   Vector of values of doses.
#' @param x      Vector of transformed values of doses (e.g. log).
#' @param p_DLT  Target toxicity rate.
#' @param n_size Cohort size.
#' @param n_max  Maximum sample size per dose.
#' @param ...    Other inputs for BOIN design. See the inputs 
#'               from get.boundary().
#' 
#' @return Dose-escalation phase data.
#' 
#' @examples
#' set.seed(1000)
#' 
#' gamma  = 0.5
#' p      = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q      = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose   = (1:6)*0.1
#' x      = dose
#' p_DLT  = 0.3
#' n_size = 3
#' n_max  = 15
#' 
#' sl_data_escal_sim(p, q, gamma, dose, x, p_DLT, n_size, n_max)
#' 
#' @export
sl_data_escal_sim = function(p, q, gamma, dose, x, 
                             p_DLT  = 0.3, 
                             n_size = 3, 
                             n_max  = 15, ...){
    
    num_DL   = length(x) # number of DL's
    n_cohort = n_max/n_size*num_DL
    
    # BOIN design decision table
    boundary = get.boundary(target      = p_DLT, 
                            ncohort     = n_cohort, 
                            cohortsize  = n_size, 
                            n.earlystop = n_max, ...)
    
    DT = boundary$boundary_tab
    
    P_rho = sl_multinom_par(gamma, p, q)
    P     = P_rho$P
    
    num_DL     = length(x)
    DL_cur     = 1
    DL_U_cur   = num_DL + 1
    data_escal = NULL
    for (j in 1:n_cohort){
        
        # data for each run
        Y            = t(rmultinom(n = 1, size = n_size, prob = P[DL_cur,]))
        colnames(Y)  = c("y_11","y_10","y_01","y_00")
        row.names(Y) = NULL
        y_T          = Y[1] + Y[2]
        y_E          = Y[1] + Y[3]
        
        data_each = data.frame(DL       = DL_cur, 
                               dose     = dose[DL_cur],
                               x        = x[DL_cur],
                               p        = p[DL_cur],
                               q        = q[DL_cur],
                               gamma    = gamma,
                               Y, y_T, y_E, n_size, 
                               y_T_cum  = NA,
                               y_E_cum  = NA,
                               n_cum    = NA,
                               decision = "", 
                               DL_U     = DL_U_cur,
                               DL_next  = NA,
                               stage    = 1)
        
        # include cumsum
        data_escal = rbind(data_escal, data_each) %>% 
            group_by(DL) %>% 
            mutate(y_T_cum = cumsum(y_T),
                   y_E_cum = cumsum(y_E),
                   n_cum   = cumsum(n_size)) %>%
            group_by(DL, stage) %>% 
            mutate(y_T_cum_s = cumsum(y_T),
                   y_E_cum_s = cumsum(y_E),
                   n_cum_s   = cumsum(n_size)) %>%
            ungroup(DL)
        
        # stop and remove the last row since > n_max at current DL
        if ((data_escal %>% slice(n()))$n_cum > n_max){
            data_escal = data_escal %>% slice(-n())
            break
        }
        
        data_escal = data_escal %>% 
            slice(-n()) %>% 
            rbind(data_escal %>% 
                      slice(n()) %>%
                      mutate(decision = case_when(y_T_cum <= DT[2, DT[1,] == n_cum] ~ "E",
                                                  y_T_cum >  DT[2, DT[1,] == n_cum] & y_T_cum < DT[3, DT[1,] == n_cum] ~ "S",
                                                  y_T_cum >= DT[3, DT[1,] == n_cum] & y_T_cum < DT[4, DT[1,] == n_cum] ~ "D",
                                                  y_T_cum >= DT[4, DT[1,] == n_cum] ~ "U"),
                             DL_U     = case_when(decision == "U" ~ DL,
                                                  .default = DL_U),
                             DL_next  = case_when(decision == "E" ~ DL + 1,
                                                  decision == "S" ~ DL,
                                                  decision == "D" | decision == "U" ~ DL - 1),
                             DL_next  = case_when(DL_next == DL_U |  DL_next > num_DL ~ DL_U - 1,
                                                  DL_next == 0 ~ 1,
                                                  .default = DL_next)))
        
        DL_cur   = (data_escal %>% slice(n()))$DL_next
        DL_U_cur = (data_escal %>% slice(n()))$DL_U
        
        if (DL_U_cur == 1){
            break
        } 
        
    }
    
    return(data_escal)
}

#' @title Select maximum tolerated dose (MTD).
#' 
#' @description Select MTD give the dose-escalation data.
#' 
#' @param data_escal Dose-escalation data from function sl_data_escal_sim().
#' @param p_DLT      Target toxicity rate, which should be the same as the 
#'                   input in function sl_data_escal_sim().
#' 
#' @return Selected dose level (DL) of the MTD.
#' 
#' @examples
#' set.seed(1000)
#' 
#' gamma  = 0.5
#' p      = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q      = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose   = (1:6)*0.1
#' x      = dose
#' p_DLT  = 0.3
#' n_size = 3
#' n_max  = 15
#' 
#' data_escal = sl_data_escal_sim(p, q, gamma, dose, x, p_DLT, n_size, n_max)
#' 
#' sl_select_MTD(data_escal, p_DLT)
#' 
#' @export
sl_select_MTD = function(data_escal, p_DLT = 0.3){
    
    # dose upper bound: all doses >= DL_U should be eliminated
    DL_U = (data_escal %>% slice(n()))$DL_U
    
    data_last = data_escal %>% 
        group_by(DL) %>% 
        slice(n()) %>% 
        dplyr::select(DL, dose, x, y_T_cum, y_E_cum, n_cum)
    
    select_mtd = select.mtd(target = p_DLT, 
                            npts   = data_last$n_cum, 
                            ntox   = data_last$y_T_cum)
    
    MTD = min(select_mtd$MTD, DL_U - 1)
    
    return(MTD)
}


#' @title Uniformly assign subjects to each dose level (DL).
#' 
#' @description Get the number of subjects assigned to each DL uniformly.
#' 
#' @param n  Total number of subjects.
#' @param m  Number of DLs.
#' 
#' @return Number of subjects assigned to each DL.
#' 
#' @examples
#' sl_ss_assign_fun(30, 4)
#' 
#' @export
sl_ss_assign_fun = function(n, m){
    
    n_each     = c()
    n_2        = n
    n_each_cur = ceiling(n_2 / m)
    for (i in 1:(m - 1)){
        n_each = c(n_each, n_each_cur)
        
        n_2    = n_2 - n_each_cur
        div    = m - i
        if (n_2 %% div == 0){
            n_each = c(n_each, rep(n_2 / div, div))
            break
        }
    }
    
    return(n_each)
    
}


#' @title Simulate dose-expansion phase data using skip-levels (SL) design 
#' or uniformly-assign (UA) design.
#' 
#' @description Simulate dose-expansion phase data using skip-levels (SL) 
#' design or uniformly-assign (UA) design given the dose level (DL) of the 
#' maximum tolerated dose (MTD) and the number of skip.
#' 
#' @param p      Vector of the toxicity rates.
#' @param q      Vector of the response rates.
#' @param gamma  Odds ratio of the response rates in toxic group versus that 
#'               in non-toxic group.
#' @param dose   Vector of values of doses.
#' @param x      Vector of transformed values of doses (e.g. log).
#' @param n      Total number of subjects assigned to the expansion phase.
#' @param MTD    The DL of the MTD. Get it from function sl_select_MTD().
#' @param skip   The number of DLs to be skipped.
#' @param method Two options: "SL" (the SL design) or "UA" (the UA design).
#' 
#' @return data_expan: Dose-expansion phase data; CDR: RP2D candidates region.
#' 
#' @examples
#' set.seed(1000)
#' 
#' gamma  = 0.5
#' p      = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q      = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose   = (1:6)*0.1
#' x      = dose
#' n      = 36
#' MTD    = 5
#' skip   = 2
#' method = "SL"
#' 
#' sl_data_expan_sim(p, q, gamma, dose, x, MTD, n, skip, method)
#' 
#' @export
sl_data_expan_sim = function(p, q, gamma, dose, x, MTD, 
                             n      = 36, 
                             skip   = 2, 
                             method = c("SL", "UA")){
    
    data_expan = NULL
    CDR        = MTD
    
    if (MTD > 1){
        # update number of DLs to be skipped if MTD is low.
        skip_update = min(skip, MTD - 2) 
        MTD_skip    = MTD - (skip_update + 1)
        
        # get expansion data
        if (method == "SL"){
            loc = c(MTD_skip, MTD)
            
        } else if (method == "UA"){
            loc = MTD_skip : MTD
            
        } else {
            stop("The method must be either 'SL' or 'UA'.")
            
        }
        
        # candidate region
        CDR = loc[1]:loc[length(loc)]
        
        # equally assign subjects
        n_each = sl_ss_assign_fun(n, length(loc))
        
        P_rho = sl_multinom_par(gamma, p, q)
        P     = P_rho$P
        
        # skip levels
        Y = t(apply(cbind(P[loc,], n_each),1,function(x){
            rmultinom(n = 1, size = x[5], prob = x[1:4])
        }))
        Y = data.frame(Y)
        colnames(Y) = c("y_11","y_10","y_01","y_00")
        row.names(Y) = NULL
        
        y_T = Y[,1] + Y[,2]
        y_E = Y[,1] + Y[,3]
        
        data_expan = tibble(DL        = loc,
                            dose      = dose[loc],
                            x         = x[loc],
                            p         = p[loc],
                            q         = q[loc],
                            gamma     = gamma,
                            Y, y_T, y_E, 
                            n_size    = n_each,
                            y_T_cum   = NA,
                            y_E_cum   = NA,
                            n_cum     = NA,
                            decision  = "",
                            DL_U      = NA,
                            DL_next   = NA,
                            stage     = 2,
                            y_T_cum_s = NA,
                            y_E_cum_s = NA, 
                            n_cum_s   = NA)
        
    }
    
    list(data_expan = data_expan,
         CDR        = CDR)
    
}


#' @title Combined data of the dose-escalation phase data and 
#' the dose-expansion phase data using the skip-levels (SL) design or 
#' the uniformly-assign (UA) design.
#' 
#' @description Combine the dose-escalation phase data and 
#' the dose-expansion phase data using the skip-levels (SL) design or 
#' the uniformly-assign (UA) design.
#' 
#' @param data_escal  Dose-escalation phase data. Simulated from function 
#'                    sl_data_escal_sim(). Users can use other 
#'                    dose-escalation design besides BOIN with the same 
#'                    structure as the output of function sl_data_escal_sim().
#' @param data_expan  Dose-expansion phase data. Simulated from function 
#'                    sl_data_expan_sim().
#' 
#' @return data_dede: dose-escalation and dose-expansion phases combined data;
#' data_dede_DL: the summary of the combined data per DL.
#' 
#' @examples
#' set.seed(1000)
#' 
#' gamma  = 0.5
#' p      = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q      = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose   = (1:6)*0.1
#' x      = dose
#' p_DLT  = 0.3
#' n_size = 3
#' n_max  = 15
#' n      = 36
#' skip   = 2
#' method = "SL"
#' 
#' # escalation data
#' data_escal = sl_data_escal_sim(p, q, gamma, dose, x, p_DLT, n_size, n_max)
#' MTD        = sl_select_MTD(data_escal, p_DLT)
#' 
#' # expansion data
#' rst_expan  = sl_data_expan_sim(p, q, gamma, dose, x, MTD, n, skip, method)
#' data_expan = rst_expan$data_expan
#' 
#' # combined data
#' sl_data_dede(data_escal, data_expan, dose, x, p, q, gamma)
#' 
#' @export
sl_data_dede = function(data_escal, data_expan,
                        dose, x, p, q, gamma){
    
    # full data: combine both phases data
    data_dede = rbind(data_escal, data_expan) %>%
        group_by(DL) %>% 
        mutate(y_T_cum = cumsum(y_T),
               y_E_cum = cumsum(y_E),
               n_cum   = cumsum(n_size)) %>%
        group_by(DL, stage) %>% 
        mutate(y_T_cum_s = cumsum(y_T),
               y_E_cum_s = cumsum(y_E),
               n_cum_s   = cumsum(n_size)) %>%
        ungroup(DL)
    
    # output final data per DL
    data_dede_DL = data_dede %>%
        ungroup(stage) %>%
        dplyr::select(DL, dose, x, p, q, gamma, y_T_cum, y_E_cum, n_cum) %>%
        group_by(DL) %>%
        slice(n()) %>%
        ungroup() %>%
        right_join(tibble(DL = 1:length(dose), dose, x, p, q, gamma))
    
    data_dede_DL[is.na(data_dede_DL)] = 0
    
    list(data_dede    = data_dede,
         data_dede_DL = data_dede_DL)
    
}


#' @title Bayesian two-parameter logistic regression model with 
#' normal prior distribution.
#' 
#' @description Fit the Bayesian two-parameter logistic regression 
#' model with normal prior distribution.
#' 
#' @param stan_data    Data frame including x: vector of values (or 
#'                     transformed values; e.g. log) of doses; 
#'                     y: vector of number of responses; 
#'                     N_n: vector of total number of subjects.
#' @param mu_prior     Vector of the normal mean for two-parameter 
#'                     logistic regression.
#' @param sigma_prior  Vector of the normal standard error for two-parameter 
#'                     logistic regression.
#' 
#' @return  post_alpha: posterior sample of intercept; post_beta: 
#' posterior sample of slope; post_q_all: posterior sample of 
#' response rates for all doses.
#' 
#' @examples
#' stan_data = data.frame(x   = (1:4)*0.1,
#'                        y   = c(1, 3, 4, 5),
#'                        N_n = c(16, 10, 10, 16))
#' mu_prior    = c(0, 0)
#' sigma_prior = c(100, 100)
#' 
#' sl_bayes_logistic_model(stan_data, mu_prior, sigma_prior)
#' 
#' @export
sl_bayes_logistic_model = function(stan_data, 
                                   mu_prior    = c(0, 0), 
                                   sigma_prior = c(100, 100)){
    
    logistic_fun = function(x, alpha, beta){
        1/(1 + exp(- alpha - beta*x))
    }
    
    lst_data = list(num_DL  = nrow(stan_data),
                    y       = stan_data$y,
                    N_n     = stan_data$N_n,
                    x       = stan_data$x,
                    mu_a    = mu_prior[1], 
                    mu_b    = mu_prior[2], 
                    sigma_a = sigma_prior[1],
                    sigma_b = sigma_prior[2])
    
    fit = rstan::sampling(object = stanmodels$Logistic, data = lst_data, refresh = 0)
    
    # posterior sample
    post_par   = rstan::extract(fit)
    
    post_alpha = post_par$alpha
    post_beta  = post_par$beta
    post_q_all = sapply(stan_data$x, function(x){ 
        logistic_fun(x, post_alpha, post_beta) 
    })
    
    list(post_alpha = post_alpha,
         post_beta  = post_beta,
         post_q_all = post_q_all)
    
}


#' @title Bayesian two-parameter logistic regression model with 
#' normal prior distribution given the combined dose-escalation
#' and dose-expansion phases data and the RP2D candidates region (CDR).
#' 
#' @description Fit the Bayesian two-parameter logistic regression model with 
#' normal prior distribution given the combined dose-escalation
#' and dose-expansion phases data and the CDR. The data with DLs within CDR
#' will be used to build the model.
#' 
#' @param data_dede_DL Data from function sl_data_dede().
#' @param CDR          RP2D candidates region; i.e. 
#'                     DLs within [MTD-(skip+1), MTD].
#' @param mu_prior     Vector of the normal mean for two-parameter 
#'                     logistic regression.
#' @param sigma_prior  Vector of the normal standard error for two-parameter 
#'                     logistic regression.
#' 
#' @return  post_alpha: posterior sample of intercept; post_beta: 
#' posterior sample of slope; post_q_all: posterior sample of 
#' response rates for all doses.
#' 
#' @examples
#' set.seed(1000)
#' 
#' gamma       = 0.5
#' p           = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q           = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose        = (1:6)*0.1
#' x           = dose
#' p_DLT       = 0.3
#' n_size      = 3
#' n_max       = 15
#' n           = 36
#' skip        = 2
#' method      = "SL"
#' mu_prior    = c(0, 0)
#' sigma_prior = c(100, 100)
#' 
#' # escalation data
#' data_escal = sl_data_escal_sim(p, q, gamma, dose, x, p_DLT, n_size, n_max)
#' MTD        = sl_select_MTD(data_escal, p_DLT)
#' 
#' # expansion data
#' rst_expan  = sl_data_expan_sim(p, q, gamma, dose, x, MTD, n, skip, method)
#' data_expan = rst_expan$data_expan
#' CDR        = rst_expan$CDR
#' 
#' # combined data
#' rst_data_dede = sl_data_dede(data_escal, data_expan, dose, x, p, q, gamma)
#' data_dede_DL  = rst_data_dede$data_dede_DL
#' 
#' # model
#' sl_bayes_logistic_model_dede(data_dede_DL, CDR, mu_prior, sigma_prior)
#' 
#' @export
sl_bayes_logistic_model_dede = function(data_dede_DL, CDR, 
                                        mu_prior    = c(0, 0), 
                                        sigma_prior = c(100, 100)){
    
    post_alpha = NULL
    post_beta  = NULL
    post_q_all = NULL
    
    if (length(CDR) > 1){
        data_model_filter = data_dede_DL %>% 
            filter(DL %in% CDR) %>%
            group_by(DL) %>% 
            slice(n())
        
        stan_data = data.frame(x    = data_model_filter$x, 
                               N_n  = data_model_filter$n_cum, 
                               y    = data_model_filter$y_E_cum)
        
        fit = sl_bayes_logistic_model(stan_data, 
                                      mu_prior, sigma_prior)
        
        post_alpha = fit$post_alpha
        post_beta  = fit$post_beta
        post_q_all = fit$post_q_all
        
    }
    
    return(list(post_alpha = post_alpha,
                post_beta  = post_beta,
                post_q_all = post_q_all))
    
}


#' @title Get the lowest dose level (DL) of the RP2D admissible region (AR).
#' 
#' @description Get the lowest dose level (DL) of the RP2D 
#' admissible region (AR) using the Bayesian two-parameter logistic
#' regression.
#' 
#' @param post_q_all  Posterior sample of response rates for all DLs from 
#'                    function sl_bayes_logistic_model_dede().
#' @param CDR         RP2D candidates region; i.e. 
#'                    DLs within [MTD-(skip+1), MTD].
#' @param delta       Relative efficiency of the response rates of the 
#'                    non-inferiority dose and the MTD. It is a 
#'                    constant in (0, 1).   
#' @param cc          Posterior probability cut to decide if the low dose is
#'                    non-inferior to the MTD. It is a constant in (0, 1).  
#' 
#' @return The lowest DL of the AR.
#' 
#' @examples 
#' set.seed(1000)
#' 
#' gamma       = 0.5
#' p           = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q           = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose        = (1:6)*0.1
#' x           = dose
#' p_DLT       = 0.3
#' n_size      = 3
#' n_max       = 15
#' n           = 36
#' skip        = 2
#' method      = "SL"
#' mu_prior    = c(0, 0)
#' sigma_prior = c(100, 100)
#' delta       = 0.9
#' cc          = 0.7
#' 
#' # escalation data
#' data_escal = sl_data_escal_sim(p, q, gamma, dose, x, p_DLT, n_size, n_max)
#' MTD        = sl_select_MTD(data_escal, p_DLT)
#' 
#' # expansion data
#' rst_expan  = sl_data_expan_sim(p, q, gamma, dose, x, MTD, n, skip, method)
#' data_expan = rst_expan$data_expan
#' CDR        = rst_expan$CDR
#' 
#' # combined data
#' rst_data_dede = sl_data_dede(data_escal, data_expan, dose, x, p, q, gamma)
#' data_dede_DL  = rst_data_dede$data_dede_DL
#' 
#' # model
#' post_par = sl_bayes_logistic_model_dede(data_dede_DL, CDR, 
#'                                         mu_prior, sigma_prior)
#' 
#' # The lowest DL of the AR
#' sl_get_ARLD(post_par$post_q_all, CDR, delta, cc)
#' 
#' @export
sl_get_ARLD = function(post_q_all, CDR, 
                       delta = 0.9, 
                       cc    = 0.7){
    
    ARLD = max(CDR)
    
    if (length(CDR) > 1){
        
        q_post_H   = post_q_all[, length(CDR)]
        for (j in 1:(length(CDR) - 1)){
            
            q_post_j = post_q_all[, j]
            
            if (mean(q_post_j/q_post_H >= delta) > cc){
                ARLD = CDR[j]
                break
                
            } else {
                next 
                
            }
        }
        
    }
    
    
    return(ARLD)
}


#' @title Simulation for dose-escalation and dose-expansion phases.
#' 
#' @description Simulate dose-escalation and dose-expansion phases data and 
#' calculate the maximum tolerated dose (MTD) and the admissible region 
#' (AR) for RP2Ds. 
#' 
#' @param p            Vector of the toxicity rates.
#' @param q            Vector of the response rates.
#' @param gamma        Odds ratio of the response rates in toxic group 
#'                     versus that in non-toxic group.
#' @param dose         Vector of values of doses.
#' @param x            Vector of transformed values of doses (e.g. log).
#' @param p_DLT        Target toxicity rate.
#' @param n_size       Cohort size.
#' @param n_max        Maximum sample size per dose.
#' @param n            Total number of subjects assigned to the expansion phase.
#' @param skip         The number of DLs to be skipped.
#' @param method       Two options: "SL" (the SL design) or "UA" 
#'                     (the UA design).
#' @param mu_prior     Vector of the normal mean prior for the slope and intercept 
#'                     of the Bayesian logistic regression.
#' @param sigma_prior  Vector of the normal standard error prior for the slope and 
#'                     intercept of the Bayesian logistic regression.
#' @param delta        Relative efficiency of the response rates of the 
#'                     non-inferiority dose and the MTD. It is a constant 
#'                     in (0, 1).   
#' @param cc           Posterior probability cut to decide if the low dose is 
#'                     non-inferior to the MTD. It is a constant in (0, 1). 
#' @param ...          Other inputs for BOIN design. See the inputs 
#'                     from get.boundary().
#' 
#' @return data_dede: combined dose-escalation and dose-expansion data; 
#' data_dede_DL: summary of the combined data; post_par: posterior sample of
#' intercept, slope and the response rates; rst_AR: the AR for each cc.
#' 
#' @examples
#' set.seed(1000)
#' 
#' gamma       = 0.5
#' p           = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55)
#' q           = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61)
#' dose        = (1:6)*0.1
#' x           = dose
#' p_DLT       = 0.3
#' n_size      = 3
#' n_max       = 15
#' n           = 36
#' skip        = 2
#' method      = "SL"
#' mu_prior    = c(0, 0)
#' sigma_prior = c(100, 100)
#' delta       = 0.9
#' cc          = c(0.7, 0.8)
#' 
#' sl_dede_full_sim(p, q, gamma, dose, x, 
#'                  p_DLT, n_size, n_max, 
#'                  n, skip, method, 
#'                  mu_prior, sigma_prior, 
#'                  delta, cc)
#' 
#' @export
sl_dede_full_sim = function(p, q, gamma, dose, x, 
                            p_DLT       = 0.3, 
                            n_size      = 3, 
                            n_max       = 15,
                            n           = 36, 
                            skip        = 2, 
                            method      = c("SL", "UA"), 
                            mu_prior    = c(0, 0), 
                            sigma_prior = c(100, 100), 
                            delta       = 0.9, 
                            cc          = c(0.7, 0.8),
                            ...){
    
    # simulate escalation phase data based on BOIN
    data_escal = sl_data_escal_sim(p, q, gamma, dose, x, 
                                   p_DLT, n_size, n_max, ...)
    
    # output MTD: isotonic regression PAVA
    MTD = sl_select_MTD(data_escal = data_escal,
                        p_DLT      = p_DLT)
    
    # simulate expansion phase data
    data_expan_CDR = sl_data_expan_sim(p, q, gamma, dose, x, MTD, 
                                       n, skip, method)
    
    data_expan = data_expan_CDR$data_expan
    
    # candidates region
    CDR = data_expan_CDR$CDR
    
    # combined data
    rst_data_dede = sl_data_dede(data_escal, data_expan,
                                 dose, x, p, q, gamma)
    
    data_dede    = rst_data_dede$data_dede
    data_dede_DL = rst_data_dede$data_dede_DL
    
    # fit Bayesian logistic regression
    post_par = sl_bayes_logistic_model_dede(data_dede_DL, CDR, 
                                            mu_prior, sigma_prior)
    
    # output the ARLD
    rst_AR = NULL
    for (j in 1:length(cc)){
        
        ARLD = sl_get_ARLD(post_par$post_q_all, CDR, delta, cc[j])
        
        AR   = ARLD:max(CDR)
        
        rst_AR = rbind(rst_AR,
                       tibble(cc     = cc[j],
                              MTD    = MTD,
                              CDR    = list(CDR),
                              method = method,
                              ARLD   = ARLD,
                              AR     = list(AR)))
        
    }
    
    list(data_dede    = data_dede,
         data_dede_DL = data_dede_DL,
         post_par     = post_par,
         rst_AR       = rst_AR)
    
}




