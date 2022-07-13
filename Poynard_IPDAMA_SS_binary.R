# Binary covariate

logit <- function(p) log(p/(1-p))

# Demidenko paper, appendix specified the information matrix in terms of M1-M4
# so start by setting this up

M1 <- matrix(c(1,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0), ncol = 4, byrow = T)
M2 <- matrix(c(1,1,0,0, 1,1,0,0, 0,0,0,0, 0,0,0,0), ncol = 4, byrow = T)
M3 <- matrix(c(1,0,1,0, 0,0,0,0, 1,0,1,0, 0,0,0,0), ncol = 4, byrow = T)
M4 <- matrix(c(1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1), ncol = 4, byrow = T)

# let x be the treatment and z be the covariate 
# define our model as logit-p = alpha + beta*x + gamma*z + interaction*x*z
# and we center our covariate by its mean, to help interpret alpha and beta

# Poynard application - 4 trials, sex covariate
DAT <- data.frame(n_C            = c(112, 89, 49, 53),
                  n_T            = c(118, 85, 30, 53),
                  age_mean_C     = c(54, 53, 55, 57),
                  age_sd_C       = c(11, 11, 9, 12),
                  age_mean_T     = c(54, 55, 53, 55),
                  age_sd_T       = c(9, 11, 7, 11),
                  percent_male_C = c(71, 73, 71, 76),
                  percent_male_T = c(71, 67, 73, 74),
                  events_C       = c(30, 28, 11, 13),
                  events_T       = c(19, 16, 1, 13))

# define vector alpha = log-odds in the control group
DAT$ln_odds_C <- logit(DAT$events_C / DAT$n_C)
alpha         <- DAT$ln_odds_C

# define vector beta = overall treatment effect
DAT$non_events_C <- DAT$n_C - DAT$events_C
DAT$non_events_T <- DAT$n_T - DAT$events_T
DAT$OR           <- (DAT$events_T / DAT$non_events_T) / (DAT$events_C / DAT$non_events_C)
DAT$lnOR         <- log(DAT$OR)
beta             <- rep(DAT$lnOR, 4)

# specify gamma = assumed prognostic effect of Z - assumed common for each trial here
gamma <-c(log(1), log(1) , log(1) , log(1))
# specify assumed interaction - assumed common for each trial here
# 1.3 is the value assumed by Kovalchick
interaction <- c(log(1.3), log(1.3), log(1.3), log(1.3))
# Kovalchick say 3.65 gets 80% power
#interaction = c(log(3.65), log(3.65), log(3.65), log(3.65))

# total sample size of each trial
DAT$total <- DAT$n_C + DAT$n_T

# define number of patients by X and Z categories
DAT$n_z0_C <- DAT$n_C * (1 - (DAT$percent_male_C / 100))
DAT$n_z0_T <- DAT$n_T * (1 - (DAT$percent_male_T / 100))
DAT$n_z1_C <- DAT$n_C * (DAT$percent_male_C / 100)
DAT$n_z1_T <- DAT$n_T * (DAT$percent_male_T / 100)

DAT$prop_z0_C <- DAT$n_z0_C / DAT$total
DAT$prop_z0_T <- DAT$n_z0_T / DAT$total
DAT$prop_z1_C <- DAT$n_z1_C / DAT$total
DAT$prop_z1_T <- DAT$n_z1_T / DAT$total

prop_z0_C <- DAT$prop_z0_C
prop_z0_T <- DAT$prop_z0_T
prop_z1_C <- DAT$prop_z1_C
prop_z1_T <- DAT$prop_z1_T

eqM1 <- array(dim = c(4, 4, nrow(DAT)))
eqM2 <- array(dim = c(4, 4, nrow(DAT)))
eqM3 <- array(dim = c(4, 4, nrow(DAT)))
eqM4 <- array(dim = c(4, 4, nrow(DAT)))

for(i in 1:nrow(DAT)){
  eqM1[,,i] <- (exp(alpha[i]) / (1 + exp(alpha[i]))^2) * M1 * prop_z0_C[i]
  eqM2[,,i] <- (exp(alpha[i] + beta[i]) / (1 + exp(alpha[i] + beta[i]))^2) * M2 * prop_z0_T[i]
  eqM3[,,i] <- (exp(alpha[i] + gamma[i]) / (1 + exp(alpha[i] + gamma[i]))^2) * M3 * prop_z1_C[i]
  eqM4[,,i] <- (exp(alpha[i] + beta[i] + gamma[i] + interaction[i]) / (1 + exp(alpha[i] + beta[i] + gamma[i] + interaction[i]))^2) * M4 * prop_z1_T[i]
  
}
Imat <- eqM1 + eqM2 + eqM3 + eqM4

vars <- vector(length = nrow(DAT), mode = 'numeric')
for(i in 1:nrow(DAT)){
  if(any(is.na(Imat[,,i])) | any(Imat[,,i]==0)){    
    vars[i] <- NA
  } else {
    vars[i] <- solve(Imat[,,i])[4,4] / DAT$total[i]
  }
}
var.summary.int <- 1/sum(vars^-1)

# Power (%) # Table 2
study.power <- vector(mode = 'numeric', length = nrow(DAT))
for(i in 1:nrow(DAT)){
  study.power[i] <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[i] / sqrt(vars[i])) + pnorm(-qnorm(0.025, lower = F) - interaction[i] / sqrt(vars[i])))
}

overall.power <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[1] / sqrt(var.summary.int)) + pnorm(-qnorm(0.025, lower = F) - interaction[1] / sqrt(var.summary.int)))

inv_var <- 1/vars
weights <- inv_var/sum(inv_var)*100

data.frame(study = c("trial_1", "trial_2", "trial_3", "trial_4"), 
           variance = round(vars, 5), 
           power = round(study.power, 2), 
           weight = round(weights, 2))

overall.power
var.summary.int

