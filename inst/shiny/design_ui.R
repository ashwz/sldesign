##-------------------------------------------------------------
##           UI FUNCTIONS
##-------------------------------------------------------------

##define the main tabset for beans
tab_main <- function() {
    tabsetPanel(type = "pills",
                id   = "mainpanel",
                tab_post_density(),
    )
}

tab_post_density <- function() {
    tabPanel("Bayesian Logistic Dose-Response Model",
             sidebarLayout(
                 sidebarPanel(
                     ## ------------------------------------------------
                     ## Dose escalation data
                     ## ------------------------------------------------
                     h4("Dose Escalation Data"),
                     textInput("x_vec", "Doses (comma-separated)",
                               value = "48, 96, 200"),
                     textInput("dl_vec", "Dose level labels (comma-separated text)",
                               value = "48mg, 96mg, 200mg"),
                     radioButtons("dose_trans", "Dose transformation",
                                  choices  = c("Log (natural log)" = "log",
                                               "None" = "none"),
                                  selected = "log",
                                  inline   = TRUE),
                     textInput("n_vec", "Dose escalation: sample size per dose (comma-separated)",
                               value = "7, 8, 12"),
                     textInput("y_vec", "Dose escalation: responses per dose (comma-separated)",
                               value = "2, 5, 10"),
                     
                     ## ------------------------------------------------
                     ## Dose expansion data
                     ## ------------------------------------------------
                     hr(),
                     h4("Dose Expansion Data"),
                     radioButtons("use_exp", "Include dose expansion data?",
                                  choices  = c("Yes" = "yes", "No" = "no"),
                                  selected = "yes",
                                  inline   = TRUE),
                     
                     conditionalPanel(
                         condition = "input.use_exp == 'yes'",
                         radioButtons("design", "Dose expansion design",
                                      choices  = c("SL" = "SL", "UA" = "UA"),
                                      selected = "SL",
                                      inline   = TRUE),
                         numericInput("exp_n", "Dose expansion: total sample size",
                                      value = 12, min = 1, step = 1),
                         textInput("exp_y_vec",
                                   "Dose expansion: responses per dose level (comma-separated)",
                                   value = "5, 5"),
                         uiOutput("exp_alloc_help")
                     ),
                     
                     hr(),
                     h4("Prior Hyperparameters"),
                     fluidRow(
                         column(6, numericInput("mu_alpha", "Mean prior for intercept (α)", value = 0)),
                         column(6, numericInput("mu_beta",  "Mean prior for slope (β)",  value = 0))
                     ),
                     fluidRow(
                         column(6, numericInput("sigma_alpha", "SD prior for intercept (α)", value = 100, min = 0.001)),
                         column(6, numericInput("sigma_beta",  "SD prior for slope (β)",  value = 100, min = 0.001))
                     ),
                     numericInput("seed", "Random seed (Bayesian model)",
                                  value = 123, min = 1, step = 1),
                     
                     hr(),
                     actionButton("run", "Fit Model", class = "btn-primary"),
                     hr(),
                     h4("Decision Criterion"),
                     numericInput("delta", paste0("Comparability threshold (δ) (no need to refit)"), value = 0.9,
                                  min = 0, max = 1, step = 0.01),
                     helpText("Ratio compared is q̃i / q̃m, where q̃m is the response",
                              "rate of the HIGHEST dose (the last dose entered above)",
                              "and q̃i is the response rate of each other dose level i.",
                              "Changing δ updates the results automatically without refitting."),
                 ),
                 
                 mainPanel(
                     tabsetPanel(
                         tabPanel("Response Rate",
                                  br(),
                                  h4("Ratio of Response Rates"),
                                  plotOutput("ratio_plot", height = "450px"),
                                  br(),
                                  tableOutput("ratio_summary"),
                                  hr(),
                                  h4("Response Rate per Dose Level"),
                                  plotOutput("q_plot", height = "450px"),
                                  br(),
                                  tableOutput("q_summary")),
                         tabPanel("Parameters",
                                  br(),
                                  plotOutput("alpha_beta_plot", height = "400px"),
                                  br(),
                                  tableOutput("alpha_beta_summary")),
                         tabPanel("Observed Data",
                                  br(),
                                  h4("Dose Escalation"),
                                  tableOutput("data_esc"),
                                  br(),
                                  h4("Dose Expansion"),
                                  tableOutput("data_exp"),
                                  br(),
                                  h4("Dose Escalation + Expansion"),
                                  tableOutput("data_all")),
                         tabPanel("Region of Response Rate",
                                  br(),
                                  helpText("Sufficient region of response rates for which the SL design is more",
                                           "efficient than the UA design (larger determinant of the Fisher",
                                           "information of the logistic regression).",
                                           "Uses the doses (on the model scale), the dose escalation sample sizes",
                                           "and the total dose expansion sample size entered on the left.",
                                           "Requires 'Include dose expansion data?' = Yes and at least 3 dose levels.",
                                           "The ranges update automatically when the inputs change;",
                                           "the empirical region plot is drawn when you click the button."),
                                  hr(),
                                  
                                  ## shown only when dose expansion data are used
                                  conditionalPanel(
                                      condition = "input.use_exp == 'yes'",
                                      h4("Range of Response Rates"),
                                      uiOutput("region_text"),
                                      tableOutput("region_table"),
                                      hr(),
                                      h4("Empirical region"),
                                      ## space of the grid search: used for the table below and the plot
                                      numericInput("region_by", "Space of grid search",
                                                   value = 0.001, min = 0.001, max = 0.1, step = 0.005,
                                                   width = "260px"),
                                      tableOutput("region_ind_table"),
                                      div(class = "help-block",
                                          p("Each point (q̃1, q̃m) is a pair of response rates at the lowest and",
                                            "highest dose levels; the pair defines the logistic dose-response curve.",
                                            "The plot shows the ranges of response rates for which the SL design is more",
                                            "efficient than the UA design (larger determinant of the Fisher information),",
                                            "given the doses and sample sizes entered on the left."),
                                          tags$ul(
                                              tags$li(tags$b("Theorem (blue):"),
                                                      "the sufficient range of response rates shown in the table above.",
                                                      "If q̃1 and q̃m both lie within this range, the SL design is more efficient."),
                                              tags$li(tags$b("Empirical (green):"),
                                                      "the region found by grid search over all (q̃1, q̃m) pairs for which",
                                                      "the SL design is at least as efficient as the UA design.",
                                                      "It contains the theorem region, since the theorem is sufficient but not necessary."),
                                              tags$li(tags$b("Empirical (Independent) (yellow):"),
                                                      "the range of response rates from the grid search when q̃m = 1 - q̃1",
                                                      "(symmetric around 0.5), shown as a square like the theorem range.")
                                          )),
                                      actionButton("run_region_plot", "Run empirical region plot",
                                                   class = "btn-primary"),
                                      br(), br(),
                                      plotOutput("region_plot", height = "600px")
                                  ),
                                  conditionalPanel(
                                      condition = "input.use_exp == 'no'",
                                      helpText("Dose expansion data are not used, so the region of response rate is not shown.")
                                  )),
                         tabPanel("Decision Table and Plot",
                                  br(),
                                  actionButton("run_combo", "Find response combinations",
                                               class = "btn-primary"),
                                  br(), br(),
                                  numericInput("tau", "Posterior probability threshold (τ) (no need to refit)",
                                               value = 0.7, min = 0, max = 1, step = 0.01,
                                               width = "230px"),
                                  helpText("Lowest DL in AR (Admissible Region) = the LOWEST dose level i with",
                                           "posterior probability P(q̃i / q̃m > δ) > τ. If no lower dose level meets",
                                           "the threshold, the highest dose level is selected.",
                                           "Highest dose level in AR is the highest dose level. ", 
                                           "Changing τ updates the table without refitting."),
                                  textInput("assume_rate",
                                            "Assumption of response rate per dose level (comma-separated) (no need to refit)",
                                            value = "0.5, 0.7", width = "320px"),
                                  helpText("True response rates assumed for the dose levels in the expansion:",
                                           "SL: 2 numbers (lowest dose, then highest dose); UA: m numbers",
                                           "(one per dose, in dose order).",
                                           "The probability of each combination of expansion responses is the product",
                                           "of the binomial probabilities at these assumed rates."),
                                  
                                  hr(),
                                  h4("Probability of being the lowest dose that is comparable to that of the high dose"),
                                  tableOutput("combo_sel_prob"),
                                  hr(),
                                  ## SL design only
                                  conditionalPanel(
                                      condition = "input.design == 'SL'",
                                      h4("Decision Plot"),
                                      helpText("Each block is one combination of expansion responses (x: lowest dose in the AR,",
                                               "y: highest dose in the AR), labeled and colored by the comparable dose in the AR.",
                                               "Grey block: no lower dose is comparable to the high dose",
                                               "(the high dose is selected).",
                                               "Other colors: the lowest dose that is comparable to the high dose."),
                                      plotOutput("combo_plot", height = "600px"),
                                      hr()
                                  ),
                                  h4("Decision Table"),
                                  DT::dataTableOutput("combo_table"))
                     )
                 )
             )
    )
}

# ---------------------- calculation --------------------
## Results are cached per expansion design, so switching SL <-> UA only
## re-displays the stored results (no model refit). The model is refitted only
## when "Fit Model" is clicked, and the new result replaces the one stored for the
## design selected at that moment.
## key = "SL" / "UA" (or "none" when dose expansion is not included)
cur_key <- reactive({
    if (identical(input$use_exp, "yes")) input$design else "none"
})

fit_store <- reactiveValues()    # fitted model per key
fit_msg   <- reactiveValues()    # error message per key (if the last run failed)

## Dose level labels used in all tables and plots, e.g. "48mg, 96mg, 200mg"
## (instead of DL1, DL2, DL3)
dl_labels <- reactive({
    m  <- length(parse_vec(input$x_vec))
    dl <- tryCatch(parse_dl(input$dl_vec, m), error = function(e) e)
    msg <- if (inherits(dl, "error")) conditionMessage(dl) else ""
    validate(need(!inherits(dl, "error"), msg))
    dl
})

fit_compute <- function() {
    
    x     <- parse_vec(input$x_vec)
    n_esc <- parse_vec(input$n_vec)
    y_esc <- parse_vec(input$y_vec)
    
    ## ---- combine dose-escalation and dose-expansion data ----
    ## (all input checks live in build_combined_data(), see server.R)
    cd <- tryCatch(
        build_combined_data(x       = x,
                            n_esc   = n_esc,
                            y_esc   = y_esc,
                            design  = input$design,
                            N_exp   = input$exp_n,
                            y_exp   = parse_vec(input$exp_y_vec),
                            log_dose    = identical(input$dose_trans, "log"),
                            include_exp = identical(input$use_exp, "yes")),
        error = function(e) e
    )
    
    err_msg <- if (inherits(cd, "error")) conditionMessage(cd) else ""
    validate(need(!inherits(cd, "error"), err_msg))
    
    cd$DL <- dl_labels()      # dose level labels used in tables / plots
    
    mu_prior    <- c(input$mu_alpha, input$mu_beta)
    sigma_prior <- c(input$sigma_alpha, input$sigma_beta)
    seed        <- input$seed
    validate(need(length(seed) == 1 && !is.na(seed) && seed == round(seed),
                  "The random seed must be a whole number."))
    
    withProgress(message = "Running Stan MCMC sampling...", value = 0.5, {
        
        ## model is fitted to the SUM of escalation + expansion data per DL
        ## (escalation data only if "Include dose expansion data?" = No)
        ## (cd$x_model = dose, or log(dose) if the log transformation is selected)
        fit <- run_sl_bayes_model(cd$x_model, cd$n_tot, cd$y_tot,
                                  mu_prior    = mu_prior,
                                  sigma_prior = sigma_prior,
                                  seed        = seed)
    })
    
    fit$x_model <- cd$x_model   # dose scale used in the model
    fit$x       <- cd$x         # original doses, used for labeling in tables
    
    fit$data_comb   <- cd                                  # for the "Data Used in Model" tab
    fit$include_exp <- identical(input$use_exp, "yes")
    fit
}

observeEvent(input$run, {
    key <- cur_key()
    res <- tryCatch(fit_compute(), error = function(e) e)
    if (inherits(res, "error")) {
        fit_msg[[key]]   <- trimws(conditionMessage(res))
        fit_store[[key]] <- NULL
    } else {
        fit_store[[key]] <- res
        fit_msg[[key]]   <- NULL
    }
})

## result shown in all tabs = the stored result of the currently selected design
fit_result <- reactive({
    key <- cur_key()
    msg <- fit_msg[[key]]
    validate(need(is.null(msg), msg))
    f <- fit_store[[key]]
    validate(need(!is.null(f), "Click 'Fit Model' to run the model for this design."))
    f
})

## ---- Ratio q̃i / q̃m, where H = highest dose level (last column) ----
## Matches the convention used in sl_get_ARLD(): ratio = q̃i / q̃m
ratio_result <- reactive({
    fit <- fit_result()
    m   <- length(fit$x)
    
    validate(
        need(m >= 2, "Need at least 2 dose levels to compute a ratio.")
    )
    
    H_idx <- m                      # highest dose = last dose entered
    q_H   <- fit$post_q_all[, H_idx]
    
    ratio_mat <- sapply(1:(m - 1), function(j) fit$post_q_all[, j] / q_H)
    if (is.null(dim(ratio_mat))) ratio_mat <- matrix(ratio_mat, ncol = m - 1)
    dl <- fit$data_comb$DL           # dose level labels
    colnames(ratio_mat) <- paste0(dl[1:(m - 1)], " / ", dl[H_idx])
    
    list(ratio_mat = ratio_mat,
         x         = fit$x,
         H_idx     = H_idx)
})

## ---- Decision table for all expansion-response combinations ----
## "Find response combinations" fits the Bayesian logistic model to every combination
## of expansion responses with sl_logistic_expan_comb() (sldesign package) and keeps the fit,
## one per design (SL / UA). The decision table is then created from the stored fit
## with sl_decision_table(), so changing delta or tau does not refit any model.
max_combo <- 1000    # safety cap on the number of model fits

combo_store <- reactiveValues()   # response-combination results per design
combo_msg   <- reactiveValues()

combo_compute <- function() {
    
    validate(need(identical(input$use_exp, "yes"),
                  "Set 'Include dose expansion data?' to Yes to use this tab."))
    
    x     <- parse_vec(input$x_vec)
    n_esc <- parse_vec(input$n_vec)
    y_esc <- parse_vec(input$y_vec)
    k     <- if (identical(input$design, "SL")) 2 else length(x)
    
    ## validate all inputs with a dummy (all-zero) set of expansion responses
    cd0 <- tryCatch(
        build_combined_data(x = x, n_esc = n_esc, y_esc = y_esc,
                            design = input$design, N_exp = input$exp_n,
                            y_exp = rep(0, k),
                            log_dose = identical(input$dose_trans, "log"),
                            include_exp = TRUE),
        error = function(e) e
    )
    err_msg <- if (inherits(cd0, "error")) conditionMessage(cd0) else ""
    validate(need(!inherits(cd0, "error"), err_msg))
    
    dl_exp <- which(cd0$n_exp > 0)          # DLs in the expansion
    n_per  <- cd0$n_exp[dl_exp[1]]          # expansion sample size per DL
    
    ## number of possible expansion response combinations
    n_combo <- (n_per + 1)^length(dl_exp)
    validate(need(n_combo <= max_combo,
                  sprintf("%d response combinations would need %d model fits (maximum allowed: %d). Please reduce the expansion sample size.",
                          n_combo, n_combo, max_combo)))
    
    seed <- input$seed
    validate(need(length(seed) == 1 && !is.na(seed) && seed == round(seed),
                  "The random seed must be a whole number."))
    
    ## fit the Bayesian model for all combinations of expansion responses with
    ## sl_logistic_expan_comb(): dose escalation data (N, Y), total dose expansion sample
    ## size (n), design (method), doses on the model scale (x). The progress
    ## (k / total combinations) is shown in the app and in the R console; errors of
    ## sl_logistic_expan_comb() are shown in the tab.
    set.seed(seed)
    fit <- NULL
    withProgress(message = "Evaluating response combinations...",
                 value = 0, {
                     fit <- sl_logistic_expan_comb(N           = cd0$n_esc,
                                                   Y           = cd0$y_esc,
                                                   n           = input$exp_n,
                                                   method      = input$design,
                                                   x           = cd0$x_model,
                                                   mu_prior    = c(input$mu_alpha, input$mu_beta),
                                                   sigma_prior = c(input$sigma_alpha, input$sigma_beta),
                                                   DL          = dl_labels())
                 })
    
    fit
}

observeEvent(input$run_combo, {
    key <- cur_key()
    res <- tryCatch(combo_compute(), error = function(e) e)
    if (inherits(res, "error")) {
        combo_msg[[key]]   <- trimws(conditionMessage(res))
        combo_store[[key]] <- NULL
    } else {
        combo_store[[key]] <- res
        combo_msg[[key]]   <- NULL
    }
})

## result shown in the last tab = the stored result of the currently selected design
combo_result <- reactive({
    validate(need(identical(input$use_exp, "yes"),
                  "Set 'Include dose expansion data?' to Yes to use this tab."))
    key <- cur_key()
    msg <- combo_msg[[key]]
    validate(need(is.null(msg), msg))
    r <- combo_store[[key]]
    validate(need(!is.null(r), "Click 'Find response combinations' to run this design."))
    r
})