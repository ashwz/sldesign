
rm(list = ls())
options(error = recover)
require(devtools)

# for local package
load_all("sldesign")
document("sldesign")

# for online package
# install_github("ashwz/sldesign")

require(sldesign)
ls(getNamespace("sldesign"))

set.seed(1000)

# set design parameters
gamma       = 0.5
p           = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
q           = c(0.05, 0.2, 0.35, 0.5, 0.45, 0.4)
dose        = (1:6)*0.1
x           = dose
p_DLT       = 0.3
n_size      = 3
n_max       = 15
n           = 30
skip        = 1
method      = "SL"
mu_prior    = c(0, 0)
sigma_prior = c(100, 100)
delta       = 0.9
cc          = c(0.7, 0.8)

# simulation
rst = sl_dede_full_sim(p, q, gamma, dose, x, 
                       p_DLT, n_size, n_max, 
                       n, skip, method, 
                       mu_prior, sigma_prior, 
                       delta, cc)

rst


