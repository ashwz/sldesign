## ---- packages -------------------------------------------------------
## The SL / UA design functions come from the 'sldesign' package:
##   sl_logistic_expan_comb(), sl_noninf_table(), sl_noninf_plot(), sl_noninf_prob(),
##   sl_q_range(), sl_q_range_emp_ind(), sl_q_range_plot(), sl_bayes_logistic_model_dede()
## The Stan model is compiled when the package is installed and is created once when the
## package is loaded (stanmodels), so the app compiles or loads no Stan code itself and
## every fit reuses the same model.
suppressPackageStartupMessages({
  library(sldesign)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

parse_vec <- function(txt) {
  # allow either commas or semicolons as separators, and drop empty pieces
  pieces <- strsplit(txt, "[,;]")[[1]]
  pieces <- trimws(pieces)
  pieces <- pieces[pieces != ""]         # drop blanks from trailing commas
  suppressWarnings(as.numeric(pieces))
}

## ---- Dose level labels used in all tables and plots (e.g. 48mg, 96mg, 200mg) ----
## txt: comma-separated labels; blank = DL1, DL2, ..., DLm
parse_dl <- function(txt, m) {
  if (is.null(txt)) return(paste0("DL", seq_len(m)))
  lab <- trimws(strsplit(txt, "[,;]")[[1]])
  lab <- lab[lab != ""]
  if (length(lab) == 0) return(paste0("DL", seq_len(m)))
  if (length(lab) != m)
    stop(sprintf("Please enter %d dose level labels (one per dose), or leave the box blank.", m))
  if (anyDuplicated(lab) > 0)
    stop("Dose level labels must be different.")
  lab
}

## ---- Combine dose-escalation and dose-expansion data per dose level ----
## design   : "SL" -> expansion at the LOWEST (DL1) and HIGHEST (DLm) dose levels
##            "UA" -> expansion at all m dose levels
## N_exp    : total expansion sample size, split EQUALLY across the DLs
##            in the expansion (so it must be divisible by 2 for SL, by m for UA)
## y_exp    : expansion responses, one per DL in the expansion
##            (SL: 2 numbers, lowest DL then highest DL; UA: m numbers, DL1..DLm)
## log_dose : if TRUE, the model dose scale is log(dose) (natural log)
## include_exp : if FALSE, expansion inputs are ignored (not validated) and the
##            totals equal the escalation data only
## Returns one row per DL with escalation, expansion and summed n / y.
## Stops with a user-readable message if any input is invalid.
build_combined_data <- function(x, n_esc, y_esc, design, N_exp, y_exp,
                                log_dose = FALSE, include_exp = TRUE) {
  
  m <- length(x)
  if (m < 2 || any(is.na(x)))
    stop("Please enter at least 2 valid doses.")
  if (isTRUE(log_dose) && any(x <= 0))
    stop("Doses must be positive to apply the log transformation.")
  if (length(n_esc) != m || any(is.na(n_esc)))
    stop("Dose escalation sample sizes must be valid numbers with the same length as doses.")
  if (length(y_esc) != m || any(is.na(y_esc)))
    stop("Dose escalation responses must be valid numbers with the same length as doses.")
  if (any(y_esc > n_esc))
    stop("Dose escalation responses cannot exceed dose escalation sample sizes.")
  if (any(n_esc <= 0))
    stop("Dose escalation sample sizes must be positive.")
  
  n_exp_full <- rep(0, m)
  y_exp_full <- rep(0, m)
  
  if (isTRUE(include_exp)) {
    
    ## ---- expansion total sample size ----
    if (length(N_exp) != 1 || is.na(N_exp) || N_exp < 1 || N_exp != round(N_exp))
      stop("Dose expansion total sample size must be a positive whole number.")
    
    ## ---- which DLs are in the expansion ----
    if (identical(design, "SL")) {
      dl_exp <- c(1, m)              # lowest and highest dose levels
    } else if (identical(design, "UA")) {
      dl_exp <- seq_len(m)
    } else {
      stop("Please select the dose expansion design (SL or UA).")
    }
    k <- length(dl_exp)
    
    ## ---- equal allocation of the expansion sample size ----
    if (N_exp %% k != 0)
      stop(sprintf("Dose expansion total sample size (%d) must be divisible by %d (the number of dose levels in the %s expansion).",
                   as.integer(N_exp), k, design))
    n_per <- N_exp / k
    
    ## ---- expansion responses ----
    if (length(y_exp) != k || any(is.na(y_exp)))
      stop(sprintf("%s design: please enter %d valid numbers for dose expansion responses (one per dose level in the expansion).",
                   design, k))
    if (any(y_exp != round(y_exp)) || any(y_exp < 0))
      stop("Dose expansion responses must be non-negative whole numbers.")
    if (any(y_exp > n_per))
      stop(sprintf("Dose expansion responses cannot exceed the expansion sample size per dose level (%g).",
                   n_per))
    
    ## ---- expansion data per DL (0 for DLs not in the expansion) ----
    n_exp_full[dl_exp] <- n_per
    y_exp_full[dl_exp] <- y_exp
  }
  
  ## ---- assemble per-DL data ----
  data.frame(DL      = seq_len(m),
             x       = x,                              # original dose
             x_model = if (isTRUE(log_dose)) log(x) else x,   # dose scale used in model
             n_esc   = n_esc,
             y_esc   = y_esc,
             n_exp   = n_exp_full,
             y_exp   = y_exp_full,
             n_tot   = n_esc + n_exp_full,       # sum of sample sizes per DL
             y_tot   = y_esc + y_exp_full)       # sum of responses per DL
}

run_sl_bayes_model <- function(x, n_cum, y_cum,
                               mu_prior    = c(0, 0),
                               sigma_prior = c(100, 100),
                               seed        = NULL) {
  
  ## fixing the R seed makes the MCMC run reproducible (the Stan seed is
  ## drawn from R's random number generator unless one is supplied explicitly)
  if (!is.null(seed)) set.seed(seed)
  
  data_dede_DL <- tibble(
    DL      = seq_along(x),
    x       = x,
    n_cum   = n_cum,
    y_E_cum = y_cum
  )
  
  CDR <- seq_along(x)   # use all provided dose levels
  
  fit <- sl_bayes_logistic_model_dede(
    data_dede_DL = data_dede_DL,
    CDR          = CDR,
    mu_prior     = mu_prior,
    sigma_prior  = sigma_prior
  )
  
  fit$x <- x
  fit
}


options(shiny.maxRequestSize = 200*1024^2)


shinyServer(function(input, output, session) {
  
  source("design_ui.R", local = TRUE)
  
  userLog <- reactiveValues()
  
  ##--------------------------------------
  ##---------main-------------------------
  ##--------------------------------------
  output$mainpage <- renderUI({
    tab_main()
  })
  
  ##--------------------------------------
  ##---------exit-------------------------
  ##--------------------------------------
  observeEvent(input$close, {
    stopApp()})
  
  ## ---- expansion responses / assumed rates when the design changes ----
  ## The values typed for each design are remembered, so switching SL <-> UA brings
  ## back what was entered for that design (matching its stored results).
  ## First visit to a design: SL = 2 numbers ("5, 5"); UA = one number per dose ("3, 3, 3").
  saved_inputs <- reactiveValues()
  last_design  <- reactiveVal("SL")
  
  observeEvent(input$design, {
    old_d <- last_design()
    new_d <- input$design
    saved_inputs[[old_d]] <- list(y = input$exp_y_vec, q = input$assume_rate)
    
    sv <- saved_inputs[[new_d]]
    if (!is.null(sv)) {
      updateTextInput(session, "exp_y_vec",  value = sv$y)
      updateTextInput(session, "assume_rate", value = sv$q)
    } else {
      m   <- length(parse_vec(input$x_vec))
      val <- if (identical(new_d, "SL")) {
        "5, 5"
      } else {
        paste(rep(3, max(m, 1)), collapse = ", ")
      }
      updateTextInput(session, "exp_y_vec", value = val)
      
      ## default assumed response rates (last tab): 2 numbers for SL, m for UA
      val_q <- if (identical(new_d, "SL")) {
        "0.5, 0.7"
      } else {
        paste(round(seq(0.5, 0.7, length.out = max(m, 2)), 2)[seq_len(max(m, 1))], collapse = ", ")
      }
      updateTextInput(session, "assume_rate", value = val_q)
    }
    last_design(new_d)
  }, ignoreInit = TRUE)
  
  ## ---- helper text (changes with the selected design) ----
  ## SL / UA description, sample size per DL, and how many responses to enter
  output$exp_alloc_help <- renderUI({
    req(input$x_vec, input$design)
    if (!identical(input$use_exp, "yes")) return(NULL)
    m <- length(parse_vec(input$x_vec))
    N <- input$exp_n
    
    if (identical(input$design, "SL")) {
      k      <- 2
      what   <- "SL: expansion at the LOWEST and HIGHEST dose levels (2 doses)."
      enter  <- "Enter 2 response numbers (lowest dose, then highest dose)."
    } else {
      k      <- m
      what   <- sprintf("UA: expansion at all dose levels (%d doses).", m)
      enter  <- sprintf("Enter %d response numbers (one per dose, in dose order).", m)
    }
    if (k < 1 || is.null(N) || is.na(N)) return(NULL)
    
    if (N %% k == 0) {
      helpText(what,
               sprintf("The total expansion sample size (%g) is split equally across the %d dose levels: %g per dose level.",
                       N, k, N / k),
               enter)
    } else {
      helpText(style = "color:#c0392b;",
               what,
               sprintf("The total sample size %g is not divisible by %d.", N, k))
    }
  })
  
  ## ---- combined escalation + expansion data used for the model ----
  ## formatting: counts as integers, dose as entered, other numbers to 3 decimals
  f_int  <- function(v) sprintf("%d", as.integer(round(v)))
  f_dose <- function(v) as.character(v)
  f_3dec <- function(v) ifelse(is.na(v), "-", formatC(v, format = "f", digits = 3))
  
  ## range bounds shown to 3 decimals: the lower bound is rounded UP and the upper
  ## bound is rounded DOWN, so the displayed range lies inside the computed range
  ## (1e-9: guards against floating-point error, e.g. 0.2 * 1000 = 200.00000000000003)
  f_lower3 <- function(v) f_3dec(ceiling(v * 1000 - 1e-9) / 1000)
  f_upper3 <- function(v) f_3dec(floor(v * 1000 + 1e-9) / 1000)
  
  ## Table 1: dose escalation data
  output$data_esc <- renderTable({
    cd <- fit_result()$data_comb
    data.frame(
      `Dose Level`         = cd$DL,
      `Dose`               = f_dose(cd$x),
      `Dose (model scale)` = f_3dec(cd$x_model),
      `N`                  = f_int(cd$n_esc),
      `Responses`          = f_int(cd$y_esc),
      `Observed rate`      = f_3dec(cd$y_esc / cd$n_esc),
      check.names = FALSE
    )
  })
  
  ## Table 2: dose expansion data (only DLs in the expansion)
  output$data_exp <- renderTable({
    fit <- fit_result()
    validate(need(isTRUE(fit$include_exp), "Dose expansion data not included."))
    cd  <- fit$data_comb
    cd  <- cd[cd$n_exp > 0, , drop = FALSE]
    data.frame(
      `Dose Level`    = cd$DL,
      `Dose`          = f_dose(cd$x),
      `Dose (model scale)` = f_3dec(cd$x_model),
      `N`             = f_int(cd$n_exp),
      `Responses`     = f_int(cd$y_exp),
      `Observed rate` = f_3dec(cd$y_exp / cd$n_exp),
      check.names = FALSE
    )
  })
  
  ## Table 3: escalation + expansion (data the model is fitted to)
  output$data_all <- renderTable({
    cd <- fit_result()$data_comb
    data.frame(
      `Dose Level`         = cd$DL,
      `Dose`               = f_dose(cd$x),
      `Dose (model scale)` = f_3dec(cd$x_model),
      `N`                  = f_int(cd$n_tot),
      `Responses`          = f_int(cd$y_tot),
      `Observed rate`      = f_3dec(cd$y_tot / cd$n_tot),
      check.names = FALSE
    )
  })
  
  
  ## ---- non-inferiority table: sl_noninf_table() in the sldesign package ----
  ## Created from the stored model fit (sl_logistic_expan_comb()); changing delta or tau
  ## recomputes the table without refitting any model.
  combo_dt <- reactive({
    fit <- combo_result()
    
    delta <- input$delta
    validate(need(length(delta) == 1 && !is.na(delta),
                  "Please enter a valid non-inferiority threshold (δ)."))
    tau <- input$tau
    validate(need(length(tau) == 1 && !is.na(tau) && tau >= 0 && tau <= 1,
                  "The posterior probability threshold (τ) must be between 0 and 1."))
    
    sl_noninf_table(fit, delta = delta, tau = tau)
  })
  
  ## ---- total probability of each lowest DL in AR ----
  ## probability of each combination of expansion responses under the assumed true
  ## response rates (product of binomial probabilities over the dose levels in the
  ## expansion), summed over the combinations with the same lowest DL in AR
  combo_dl_prob <- reactive({
    tab <- combo_dt()
    
    k      <- length(attr(tab, "fit")$dl_exp)
    q_true <- parse_vec(input$assume_rate)
    validate(need(length(q_true) == k && !any(is.na(q_true)) && all(q_true >= 0 & q_true <= 1),
                  sprintf("Please enter %d assumed response rates between 0 and 1 (one per dose level in the expansion).", k)))
    
    ## see sl_noninf_prob() in the sldesign package
    unlist(sl_noninf_prob(tab, q_true))
  })
  
  ## ---- Region of response rate: SL design more efficient than UA design ----
  ## sl_q_range(): sufficient (theorem) range of q for the lowest / highest dose
  ## sl_q_range_plot(): theorem, empirical and independent empirical regions
  ## x = doses on the model scale (log(dose) if selected), n = total expansion
  ## sample size, N = sample size per dose level before the expansion
  ## (= dose escalation sample sizes).
  ## validated inputs of the region calculations
  region_inputs <- reactive({
    validate(need(identical(input$use_exp, "yes"),
                  "Dose expansion data are not used, so the region of response rate is not shown."))
    
    x     <- parse_vec(input$x_vec)
    n_esc <- parse_vec(input$n_vec)
    validate(need(length(x) >= 3 && !any(is.na(x)),
                  "At least 3 valid dose levels are required."))
    use_log <- identical(input$dose_trans, "log")
    validate(need(!use_log || all(x > 0), "Doses must be positive to apply the log transformation."))
    x_model <- if (use_log) log(x) else x
    validate(need(all(diff(x_model) > 0), "Doses must be in increasing order with no ties."))
    validate(need(length(n_esc) == length(x) && !any(is.na(n_esc)) && all(n_esc >= 0),
                  "Dose escalation sample sizes must be valid numbers with the same length as doses."))
    
    n_exp <- input$exp_n
    validate(need(length(n_exp) == 1 && !is.na(n_exp) && n_exp >= 1,
                  "Dose expansion total sample size must be a positive number."))
    
    by <- input$region_by      # space between the grid search points
    validate(need(length(by) == 1 && !is.na(by) && by > 0 && by < 0.5,
                  "The space of grid search (by) must be a number between 0 and 0.5 (default 0.01)."))
    
    list(x = x_model, n = n_exp, N = n_esc, by = by)
  })
  
  ## ranges (fast): update automatically when the inputs change
  region_result <- reactive({
    ri <- region_inputs()
    
    ## sufficient (theorem) range
    q_range <- suppressWarnings(sl_q_range(ri$x, ri$n, ri$N))
    
    ## independent empirical range (symmetric at q = 0.5)
    q_ind <- tryCatch(suppressWarnings(sl_q_range_emp_ind(ri$x, ri$n, ri$N, ri$by)),
                      error = function(e) c(NA_real_, NA_real_))
    if (length(q_ind) != 2 || any(!is.finite(q_ind))) q_ind <- c(NA_real_, NA_real_)
    
    list(q_range = q_range, q_ind = q_ind)
  })
  
  ## empirical region plot (slower grid search): drawn only when the button is clicked
  region_plot_result <- eventReactive(input$run_region_plot, {
    ri  <- region_inputs()
    plt <- NULL
    withProgress(message = "Computing empirical region...", value = 0.5, {
      plt <- tryCatch(suppressWarnings(sl_q_range_plot(ri$x, ri$n, ri$N, ri$by)),
                      error = function(e) e)
    })
    plt
  })
  
  output$region_text <- renderUI({
    r <- region_result()
    if (any(is.na(r$q_range))) {
      helpText("The theorem gives no sufficient region for these inputs (no response rates are guaranteed to favor the SL design).")
    } else {
      helpText("The SL design is more efficient than the UA design if the true response rates lie within the following range.")
    }
  })
  
  ## sufficient (theorem) range: one-row table
  output$region_table <- renderTable({
    r <- region_result()
    data.frame(`Range`       = "Sufficient region (theorem)",
               `Lower bound` = f_lower3(r$q_range[1]),
               `Upper bound` = f_upper3(r$q_range[2]),
               check.names = FALSE)
  })
  
  ## independent empirical range: one-row table (grid search with the same
  ## space of grid search, "by", as the empirical region plot)
  output$region_ind_table <- renderTable({
    r <- region_result()
    data.frame(`Range`       = "Independent empirical range",
               `Lower bound` = f_lower3(r$q_ind[1]),
               `Upper bound` = f_upper3(r$q_ind[2]),
               check.names = FALSE)
  })
  
  output$region_plot <- renderPlot({
    plt <- region_plot_result()
    validate(need(!inherits(plt, "error"),
                  if (inherits(plt, "error")) conditionMessage(plt) else ""))
    plt
  })
  
  ## probability of being the lower bound of the AR (sums over all combinations)
  output$combo_sel_prob <- renderTable({
    ## one-row table: one column per dose level (the lowest DL in AR)
    q   <- combo_dl_prob()
    out <- as.data.frame(matrix(f_3dec(q), nrow = 1), stringsAsFactors = FALSE)
    colnames(out) <- names(q)
    out
  })
  
  ## interactive table: each column has its own filter box / slider
  ## ---- scatter plot of the non-inferiority table (SL design only) ----
  ## x / y = expansion responses at the lowest / highest dose level;
  ## colour = lowest DL in AR (grey = no lower DL is non-inferior to the highest DL).
  output$combo_plot <- renderPlot({
    req(identical(input$design, "SL"))
    sl_noninf_plot(combo_dt())
  })
  
  output$combo_table <- DT::renderDataTable({
    tab <- combo_dt()
    attr(tab, "fit") <- NULL      # model results are not part of the displayed table
    
    ## rates and probabilities are shown with 3 decimals (filters use the full values)
    ## (all columns except the first, the responses tuple, and the last, lowest DL in AR)
    dec_cols <- colnames(tab)[-c(1, ncol(tab))]
    
    DT::datatable(tab,
                  rownames = FALSE,
                  filter   = "top",
                  options  = list(pageLength = 10,
                                  lengthMenu = c(10, 25, 50, 100),
                                  scrollX    = TRUE)) %>%
      DT::formatRound(columns = dec_cols, digits = 3)
  })
  
  ## ---- Alpha & Beta density plots ----
  output$alpha_beta_plot <- renderPlot({
    fit <- fit_result()
    
    df <- data.frame(
      value = c(fit$post_alpha, fit$post_beta),
      param = rep(c("Intercept (α)", "Slope (β)"), each = length(fit$post_alpha))
    )
    
    ggplot(df, aes(x = value, fill = param)) +
      geom_density(alpha = 0.6) +
      facet_wrap(~ param, scales = "free") +
      theme_minimal(base_size = 18) +
      labs(x = "Value", y = "Posterior Density", fill = NULL) +
      theme(legend.position = "none")
  })
  
  output$alpha_beta_summary <- renderTable({
    fit <- fit_result()
    data.frame(
      Parameter = c("Intercept (α)", "Slope (β)"),
      Mean      = c(mean(fit$post_alpha), mean(fit$post_beta)),
      SD        = c(sd(fit$post_alpha), sd(fit$post_beta)),
      `2.5%`    = c(quantile(fit$post_alpha, 0.025), quantile(fit$post_beta, 0.025)),
      `97.5%`   = c(quantile(fit$post_alpha, 0.975), quantile(fit$post_beta, 0.975)),
      check.names = FALSE
    )
  }, digits = 3)
  
  ## ---- q (response rate) density plots per dose level ----
  output$q_plot <- renderPlot({
    fit <- fit_result()
    
    lab <- fit$data_comb$DL                   # dose level labels
    
    q_df <- as.data.frame(fit$post_q_all)
    colnames(q_df) <- lab
    
    q_df <- q_df %>%
      pivot_longer(everything(), names_to = "DL", values_to = "q") %>%
      mutate(DL = factor(DL, levels = lab))
    
    dose_labels <- setNames(lab, lab)
    
    ggplot(q_df, aes(x = q, fill = DL)) +
      geom_density(alpha = 0.6) +
      facet_wrap(~ DL, scales = "free_y",
                 labeller = as_labeller(dose_labels)) +
      theme_minimal(base_size = 18) +
      labs(x = "Response Rate", y = "Posterior Density") +
      theme(legend.position = "none") +
      xlim(0, 1)
  })
  
  output$q_summary <- renderTable({
    fit <- fit_result()
    cd  <- fit$data_comb       # escalation + expansion (data used in the model)
    data.frame(
      `Dose Level` = cd$DL,
      `Dose`       = f_dose(fit$x),
      `Dose (model scale)` = f_3dec(cd$x_model),
      `N`          = f_int(cd$n_tot),
      `Responses`  = f_int(cd$y_tot),
      `Observed RR` = cd$y_tot / cd$n_tot,
      Mean         = apply(fit$post_q_all, 2, mean),
      SD           = apply(fit$post_q_all, 2, sd),
      `2.5%`       = apply(fit$post_q_all, 2, quantile, probs = 0.025),
      `97.5%`      = apply(fit$post_q_all, 2, quantile, probs = 0.975),
      check.names  = FALSE
    )
  }, digits = 3)
  
  ## ---- Ratio q̃i / q̃m density plots ----
  output$ratio_plot <- renderPlot({
    rr <- ratio_result()
    
    r_df <- as.data.frame(rr$ratio_mat) %>%
      pivot_longer(everything(), names_to = "Comparison", values_to = "ratio") %>%
      mutate(Comparison = factor(Comparison, levels = colnames(rr$ratio_mat)))
    
    ggplot(r_df, aes(x = ratio, fill = Comparison)) +
      geom_density(alpha = 0.6) +
      geom_vline(xintercept = input$delta, linetype = "dashed", color = "red") +
      facet_wrap(~ Comparison, scales = "free") +
      theme_minimal(base_size = 18) +
      labs(x = expression(tilde(q)[i] / tilde(q)[m]), y = "Posterior Density",
           caption = paste0("Dashed red line: delta = ", input$delta)) +
      theme(legend.position = "none")
  })
  
  output$ratio_summary <- renderTable({
    rr <- ratio_result()
    cd <- fit_result()$data_comb      # escalation + expansion (data used in the model)
    m  <- nrow(cd)
    j  <- 1:(m - 1)                   # lower dose levels
    H  <- rr$H_idx                    # highest dose level
    
    obs_rate <- cd$y_tot / cd$n_tot
    
    data.frame(
      Comparison   = colnames(rr$ratio_mat),
      `Dose`       = paste0(f_dose(rr$x[j]), " / ", f_dose(rr$x[H])),
      `N`          = paste0(f_int(cd$n_tot[j]), " / ", f_int(cd$n_tot[H])),
      `Responses`  = paste0(f_int(cd$y_tot[j]), " / ", f_int(cd$y_tot[H])),
      `Observed Ratio` = obs_rate[j] / obs_rate[H],     # ratio of observed rates
      Mean         = apply(rr$ratio_mat, 2, mean),
      SD           = apply(rr$ratio_mat, 2, sd),
      `2.5%`       = apply(rr$ratio_mat, 2, quantile, probs = 0.025),
      `97.5%`      = apply(rr$ratio_mat, 2, quantile, probs = 0.975),
      `P(ratio > δ)` = apply(rr$ratio_mat, 2, function(r) mean(r > input$delta)),
      check.names  = FALSE
    )
  }, digits = 3)
  
  
})