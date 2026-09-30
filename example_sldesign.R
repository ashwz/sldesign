
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
    p           = c(0.03, 0.07, 0.10, 0.13, 0.25, 0.55),
    q           = c(0.20, 0.30, 0.54, 0.57, 0.60, 0.61),
    gamma       = 0.5,
    dose        = (1:6)*0.1,
    x           = (1:6)*0.1,
    p_DLT       = 0.3,
    n_size      = 3,
    n_max       = 15,
    n           = 36,
    skip        = 2,
    method      = "SL",
    mu_prior    = c(0, 0),
    sigma_prior = c(100, 100),
    delta       = 0.9,
    cc          = c(0.7, 0.8)
)


