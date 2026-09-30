
rm(list = ls())
options(error = recover)
require(devtools)

# for local package
load_all("sldesign")
document("sldesign")

# for online package
# library(devtools)
# install_github("ashwz/sldesign")

require(sldesign)
ls(getNamespace("sldesign"))

# simulate a single trial
trial = sl_dede_full_sim(
    p           = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55), # toxicity rates
    q           = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61), # response rates
    gamma       = 0.5,         # odds ratio
    dose        = (1:6)*0.1,   # doses
    x           = (1:6)*0.1,   # transformed doses
    p_DLT       = 0.3,         # target toxicity rate
    n_size      = 3,           # cohort size
    n_max       = 15,          # maximum sample size per dose
    n           = 36,          # total number of subjects for expansion
    skip        = 2,           # number of DLs to be skipped
    method      = "SL",        # SL: the SL design; UA: the UA design
    mu_prior    = c(0, 0),     # mean prior for parameters of logistic regression
    sigma_prior = c(100, 100), # sd prior for parameters of logistic regression
    delta       = 0.9,         # relative efficiency for non-inferiority
    cc          = c(0.7, 0.8)  # posterior probability cut
)


