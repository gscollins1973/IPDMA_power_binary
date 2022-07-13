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

# STEER-OR example
DAT <- data.frame(n_C = c(68, 24, 44, 99, 74, 159, 63, 41, 43, 16, 108, 23, 35, 49, 140, 147, 54, 78, 40, 44, 9, 11, 54, 54, 20, 17, 102, 156, 27, 102, 23), 
                  n_T = c(142, 24, 45, 100, 148, 153, 63, 111, 45, 23, 109, 25, 36, 56, 278, 71, 53, 80, 40, 43, 19, 21, 55, 55, 20, 17, 101, 235, 28, 98, 23), 
                  percent_male_C = 100*c(0.220588237047195, 0.20833332836628, 0.545454561710358, 0.30303031206131, 0.364864856004715, 0.377358496189117, 0.269841283559799, 0.170731708407402, 0.465116292238235, 0.25, 0.351851850748062, 0.260869562625885, 0.314285725355148, 0.122448980808258, 0.314285725355148, 0.591836750507355, 0.462962955236435, NA, NA, 0, 0.444444447755814, 0.0909090936183929, 0.481481492519379, 0.277777791023254, 0.349999994039536, 0.235294118523598, 0.450980395078659, 0.358974367380142, 0.333333343267441, 0.205882355570793, 0.52173912525177), 
                  percent_male_T = 100*c(0.309859156608582, 0.375, 0.488888889551163, 0.400000005960464, 0.283783793449402, 0.366013079881668, 0.222222223877907, 0.297297298908234, 0.244444444775581, 0.217391297221184, 0.348623842000961, 0.119999997317791, 0.333333343267441, 0.214285716414452, 0.287769794464111, 0.591549277305603, 0.433962255716324, NA, NA, 0, 0.473684221506119, 0.142857149243355, 0.436363637447357, 0.345454543828964, 0.0500000007450581, 0.235294118523598, 0.376237630844116, 0.361702114343643, 0.214285716414452, 0.224489793181419, 0.60869562625885), 
                  age_mean_C = c(64.2544097900391, 67.1666641235352, 64.5926132202148, 62.5353546142578, 62.2837829589844, 69.626579284668, 64.9354858398438, 69.6097564697266, 60.3720932006836, 75.1875, 68.2314834594727, 61.3025016784668, 61.5309524536133, 65.2040786743164, 66.7642822265625, 59.1088447570801, 63.5644454956055, 68.6059341430664, 58.2000007629395, 63.9090919494629, 70.4444427490234, 71.1818161010742, 57.2055549621582, 68.6037750244141, 67.1031951904297, 70.7647094726562, 66.5882339477539, 61.8910255432129, 78.9259262084961, 67.7843170166016, 66.7826080322266), 
                  age_sd_C = c(12.2120180130005, 8.13295364379883, 7.55268573760986, 5.36107110977173, 6.77327585220337, 6.26229667663574, 9.42941188812256, 6.09868049621582, 9.91397094726562, 4.57848215103149, 7.97614002227783, 7.06486034393311, 7.79656314849854, 5.73359823226929, 8.7174768447876, 9.98500919342041, 8.72769832611084, 6.13285255432129, 4.25591611862183, 2.36070728302002, 7.82801246643066, 5.25010824203491, 9.81989574432373, 6.87341499328613, 5.35772657394409, 4.71075105667114, 9.56204605102539, 9.58733940124512, 8.30164909362793, 9.22290229797363, 7.27359199523926), 
                  age_mean_T = c(65.2908477783203, 65.1666641235352, 64.51611328125, 60.9599990844727, 63.9324340820312, 69.8562088012695, 63.2063484191895, 70.3873901367188, 61.7555541992188, 73.5652160644531, 67.9357833862305, 65.0339202880859, 63.3449058532715, 65.4285736083984, 66.5215835571289, 57.7887306213379, 65.6188659667969, 69.0233383178711, 57.25, 63.8139533996582, 66.3684234619141, 72.4285736083984, 58.403636932373, 66.8545455932617, 66.1433334350586, 69.5882339477539, 64.1782150268555, 61.5446815490723, 78.8928604125977, 68.2857131958008, 67.5652160644531), 
                  age_sd_T = c(11.4627475738525, 6.72223234176636, 9.05132865905762, 5.92021036148071, 9.39363098144531, 6.82045888900757, 8.38423728942871, 6.28297996520996, 9.48720550537109, 7.29773378372192, 8.54213428497314, 8.91479015350342, 9.54590320587158, 5.27380180358887, 8.24651145935059, 10.4182748794556, 8.2117919921875, 6.54862785339355, 3.97911214828491, 2.41282105445862, 5.549880027771, 6.52358341217041, 9.99824047088623, 7.40974855422974, 8.73675537109375, 6.73664236068726, 8.52454853057861, 9.57522678375244, 6.90860605239868, 8.44997406005859, 7.86144828796387), 
                  events_C = c(39, 9, 18, 40, 40, 55, 22, 19, 16, 9, 60, 10, 13, 20, 71, 78, 16, 44, 17, 29, 6, 5, 27, 21, 10, 11, 66, 78, 16, 55, 4), 
                  events_T = c(88, 15, 30, 64, 75, 56, 40, 83, 28, 8, 73, 7, 29, 32, 193, 50, 33, 40, 26, 30, 13, 16, 38, 23, 15, 10, 64, 110, 20, 77, 7))

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
gamma <-c(rep(log(1), nrow(DAT)))
# specify assumed interaction - assumed common for each trial here
interaction <- c(rep(log(1.3), nrow(DAT)))

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

var.summary.int <- 1/sum(vars^-1, na.rm=T)

# Power (%) # Table 2
study.power <- vector(mode = 'numeric', length = nrow(DAT))
for(i in 1:nrow(DAT)){
  study.power[i] <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[i] / sqrt(vars[i])) + pnorm(-qnorm(0.025, lower = F) - interaction[i] / sqrt(vars[i])))
}

overall.power <- 100 * (pnorm(-qnorm(0.025, lower = F) + interaction[1] / sqrt(var.summary.int)) + pnorm(-qnorm(0.025, lower = F) - interaction[1] / sqrt(var.summary.int)))

inv_var <- 1/vars
weights <- inv_var/sum(inv_var, na.rm=T)*100

data.frame(study = paste("trial_", seq(1, nrow(DAT), 1), sep=''), 
           variance = round(vars, 5), 
           power = round(study.power, 2), 
           weight = round(weights, 2))

overall.power
var.summary.int

