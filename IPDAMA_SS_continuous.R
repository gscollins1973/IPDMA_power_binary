# Continuous covariate

# Demidenko paper, appendix specified the information matrix
# this is the expected value of a 4 by 4 matrix
# The expected value of a matrix is define ｩas
# the matrix of expected values
# so we will generate a large dataset, and derive the values of each element 
# for each patient, and then take their mean.

# Poynard application - 4 trials, age covariate
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

# first define the logistic equation parameters for each study

# define vector alpha = log-odds in the control group
DAT$ln_odds_C <- log(DAT$events_C/DAT$n_C)
alpha         <- DAT$ln_odds_C

# define vector beta = overall treatment effect
DAT$non_events_C <- DAT$n_C - DAT$events_C
DAT$non_events_T <- DAT$n_T - DAT$events_T
DAT$OR           <- (DAT$events_T / DAT$non_events_T) / (DAT$events_C / DAT$non_events_C)
DAT$lnOR         <- log(DAT$OR)

beta <- rep(DAT$lnOR, 4)

# specify gamma = assumed prognostic effect of Z - assumed common for each trial here
gamma <- c(log(1), log(1) , log(1) , log(1)) 
# gamma <- c(log(1.025), log(1.025), log(1.025), log(1.025)) 

# specify assumed interaction - assumed common for each trial here
# 1.3 is the OR value assumed by Kovalchick for a 10-unit increase in age
# corresponds to a lnOR of 0.262364
# - for a 1-year increase in age, this corresponds to 0.0262364

interaction <- c(0.0262364, 0.0262364, 0.0262364, 0.0262364) 
# Kovalchick says to gets 80% power the OR is 1.74 for a 10 year increase
# this corresponds to a logOR of 0.55388511
# - and a lnOR of 0.055388511 for a 1-year increase
# interaction <- c(0.055388511, 0.055388511, 0.055388511, 0.055388511) 

# total sample size of each trial
DAT$total = DAT$n_C + DAT$n_T

# generate the covariate Z of interest;  assume normal for now & that we 
# know the distribution for treat and control groups separately
obs     <- 1000000
treat   <- matrix(0, ncol = nrow(DAT), nrow = obs)
control <- matrix(0, ncol = nrow(DAT), nrow = obs)
id      <- seq(1, obs, 1)
for(i in 1:nrow(DAT)){
  x1 <- treat[,i]
  x1[id > obs * DAT$n_C[i] / DAT$total[i]] <- 1
  treat[,i] <- x1
  control[x1 == 0, i] <- 1
}

z <- matrix(0, ncol = nrow(DAT), nrow = obs)
for(i in 1:nrow(DAT)){
  x1 <- treat[,i]
  x2 <- rep(0, obs)
  x2[x1 == 0] <- rnorm(obs - sum(treat[,i]), DAT$age_mean_C[i], DAT$age_sd_C[i])
  x2[x1 == 1] <- rnorm(sum(treat[,i]), DAT$age_mean_T[i], DAT$age_sd_T[i])
  z[,i] <- x2
}

z_cent <- matrix(ncol = 4, nrow = obs)
for(i in 1:nrow(DAT)){
  z_cent[,i] <- scale(z[,i], scale = F)  
}

# now generate the 4 by 4 matrix entries corresponding to logit-p = alpha + beta*x + gamma*z + interaction*z*x
# for each study separately 

M_11 <- matrix(ncol = nrow(DAT), nrow = obs)
M_12 <- matrix(ncol = nrow(DAT), nrow = obs)
M_13 <- matrix(ncol = nrow(DAT), nrow = obs)
M_14 <- matrix(ncol = nrow(DAT), nrow = obs)
M_21 <- matrix(ncol = nrow(DAT), nrow = obs)
M_22 <- matrix(ncol = nrow(DAT), nrow = obs)
M_23 <- matrix(ncol = nrow(DAT), nrow = obs)
M_24 <- matrix(ncol = nrow(DAT), nrow = obs)
M_31 <- matrix(ncol = nrow(DAT), nrow = obs)
M_32 <- matrix(ncol = nrow(DAT), nrow = obs)
M_33 <- matrix(ncol = nrow(DAT), nrow = obs)
M_34 <- matrix(ncol = nrow(DAT), nrow = obs)
M_41 <- matrix(ncol = nrow(DAT), nrow = obs)
M_42 <- matrix(ncol = nrow(DAT), nrow = obs)
M_43 <- matrix(ncol = nrow(DAT), nrow = obs)
M_44 <- matrix(ncol = nrow(DAT), nrow = obs)
LP   <- matrix(ncol = nrow(DAT), nrow = obs)

for(i in 1:nrow(DAT)){
  LP[, i]  <- alpha[i] + (beta[i]*treat[,i]) + (gamma[i] * z_cent[,i]) + (interaction[i] * z_cent[,i] * treat[,i])
  M_11[,i] <- exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_12[,i] <- treat[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_13[,i] <- z_cent[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_14[,i] <- treat[,i] * z_cent[,i] * exp(LP[,i])/((1 + exp(LP[,i]))^2)
  M_21[,i] <- M_12[,i]
  M_22[,i] <- treat[,i] * treat[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_23[,i] <- treat[,i] * z_cent[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_24[,i] <- treat[,i] * treat[,i] * z_cent[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_31[,i] <- M_13[,i]
  M_32[,i] <- M_23[,i]
  M_33[,i] <- z_cent[,i] * z_cent[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_34[,i] <- z_cent[,i] * z_cent[,i] * treat[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
  M_41[,i] <- M_14[,i]
  M_42[,i] <- M_24[,i]
  M_43[,i] <- M_34[,i]
  M_44[,i] <- z_cent[,i] * treat[,i] * z_cent[,i] * treat[,i] * exp(LP[,i]) / ((1 + exp(LP[,i]))^2)
}


Imat <- array(dim = c(4, 4, nrow(DAT)))
for(i in 1:nrow(DAT)){
  Imat[1,1,] <- colMeans(M_11)
  Imat[1,2,] <- colMeans(M_12)
  Imat[1,3,] <- colMeans(M_13)
  Imat[1,4,] <- colMeans(M_14)
  Imat[2,1,] <- colMeans(M_21)
  Imat[2,2,] <- colMeans(M_22)
  Imat[2,3,] <- colMeans(M_23)
  Imat[2,4,] <- colMeans(M_24)
  Imat[3,1,] <- colMeans(M_31)
  Imat[3,2,] <- colMeans(M_32)
  Imat[3,3,] <- colMeans(M_33)
  Imat[3,4,] <- colMeans(M_34)
  Imat[4,1,] <- colMeans(M_41)
  Imat[4,2,] <- colMeans(M_42)
  Imat[4,3,] <- colMeans(M_43)
  Imat[4,4,] <- colMeans(M_44)
}

Imat_inv <- array(dim=c(4, 4, nrow(DAT)))
for(i in 1:nrow(DAT)){
  Imat_inv[,,i] <- solve(Imat[,,i])
}

vars <- vector(mode = 'numeric', length = nrow(DAT))
X    <- array(dim = c(4, 4, nrow(DAT)))
for(i in 1:nrow(DAT)){
  X[,,i]  <- Imat_inv[,,i] / DAT$total[i]
  vars[i] <- X[4,4,i]
}

overall.var <- 1 / sum((c(vars))^-1)

# Power (%) Table 2 in the paper
study.power <- c(100 * (pnorm(-qnorm(0.025, lower = F) + interaction[1] / sqrt(vars[1])) + pnorm(-qnorm(0.025, lower = F) - interaction[1] / sqrt(vars[1]))),
                 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[2] / sqrt(vars[2])) + pnorm(-qnorm(0.025, lower = F) - interaction[2] / sqrt(vars[2]))),
                 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[3] / sqrt(vars[3])) + pnorm(-qnorm(0.025, lower = F) - interaction[3] / sqrt(vars[3]))),
                 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[4] / sqrt(vars[4])) + pnorm(-qnorm(0.025, lower = F) - interaction[4] / sqrt(vars[4]))))

overall.power <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[1] / sqrt(overall.var)) + pnorm(-qnorm(0.025, lower = F) - interaction[1] / sqrt(overall.var)))

study.power
overall.power

