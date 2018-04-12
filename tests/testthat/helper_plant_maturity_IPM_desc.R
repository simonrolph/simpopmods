################ IPM_desc

states <- c("immature", "mature")
states_z <- c(T,T)

# $par_info
par_info <- data.frame(t(data.frame(
  ## survival
  surv.int  =  c(par = "surv.int",min = -10,max = 4,samp = T),
  surv.z    =   c(par = "surv.z",min = 0,max = 10,samp = T),

  ## flowering
  flow.int  = c(par = "flow.int",min = -10,max = 10,samp = T),
  flow.z    =   c(par = "flow.z",min = 0,max = 10,samp = T),

  ## growth
  grow.int  =   c(par = "grow.int",min = 0,max = 1,samp = T),
  grow.z    =   c(par = "grow.z",min = 0,max = 1,samp = F), # constrained from 0 to 1
  grow.sd   =   c(par = "grow.sd",min = 0.01,max = 0.5,samp = T),

  ## recruit size
  rcsz.int  =   c(par = "rcsz.int",min = 0,max = 0,samp = F), # doesn't change
  rcsz.sd   =   c(par = "rcsz.sd",min = 0.1,max = 0.5,samp = T),

  ## seed size
  seed.int  =   c(par = "seed.int",min = 0,max = 100,samp = T),
  seed.z    =   c(par = "seed.z",min = 0,max = 10,samp = T),

  # maturity
  matu.int =   c(par = "matu.int",min = 0,max = 1,samp = T),
  matu.z   =   c(par = "matu.z",min = 0,max = 1,samp = T),

  ## recruitment probability
  p.r       =   c(par = "p.r",min = 0,max = 0.99,samp = T) # constrained from 0 to 1
)),stringsAsFactors=F)

par_info$min <- as.numeric(par_info$min)
par_info$max <- as.numeric(par_info$max)

# $par_fns
par_fns <- list(
  rcsz.int = function(params){0},
  grow.z = function(params){1 - params["grow.int"]}
)

# $demo_fns
demo_fns <- list(
  ## G_z1z -  Growth function, given you are size z now returns the pdf of size z1 next time
  grow = function(z1, z, params){
    mu <- params["grow.int"] + params["grow.z"] * z           # mean size next year
    sig <- params["grow.sd"]                                 # sd about mean
    p.den.grow <- dnorm(z1, mean = mu, sd = sig)             # pdf that you are size z1 given you were size z
    return(p.den.grow)
  },

  ## s_z Survival function, logistic regression
  surv = function(z, params){

    linear.p <- params["surv.int"] + params["surv.z"] * z  # linear predictor
    p <- pmin(1/(1+exp(-linear.p)),rep(0.99,length(z))) # logistic transformation to probability with constant cap
    p <- 1/(1+exp(-linear.p))
    return(p)
  },

  ## p_bz Probability of flowering function, logistic regression
  flow = function(z, params){
    linear.p <- params["flow.int"] + params["flow.z"] * z      # linear predictor
    p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
    return(p)
  },

  ## b_z Seed production function
  seed = function(z, params){
    N <- exp(params["seed.int"] + params["seed.z"] * z)   # seed production of a size z plant
    return(N)
  },

  ## c_0z1 Recruit size pdf
  rcsz = function(z1, params){
    mu <- params["rcsz.int"]
    sig <- params["rcsz.sd"]
    p.deRecr <- dnorm(z1, mean = mu, sd = sig)              # pdf of a size z1 recruit
    return(p.deRecr)
  },

  matu = function(z, params){
    linear.p <- params["matu.int"] + params["matu.z"] * z      # linear predictor
    p <- 1/(1+exp(-linear.p))                                # logistic transformation to probability
    return(p)
  }
)

kernels <- c("P","Fec")

# initalise kernels nested list
kernel_fns <- init_nested_list(kernels,states)

# kernel_fns$P$   TO   $   FROM

# SURVIVAL/GROWTH -------------------------

# immature to immature
kernel_fns$P$immature$immature <- function (z1, z, params) {
  # probability of survival * probability of not maturing * growth pdf
  return(IPM_desc_maturity$demo_fns$surv(z, params) * (1-IPM_desc_maturity$demo_fns$matu(z, params)) * IPM_desc_maturity$demo_fns$grow(z1, z, params))
}

# immature to mature survival (recruitment)
kernel_fns$P$mature$immature <- function (z1, z, params) {
  # probability of survival * probability of maturing * growth pdf
  return(IPM_desc_maturity$demo_fns$surv(z, params) * IPM_desc_maturity$demo_fns$matu(z, params) * IPM_desc_maturity$demo_fns$grow(z1, z, params))
}

# mature to mature survival
kernel_fns$P$mature$mature <- function (z1, z, params) {
  # survival * growth pdf
  return(IPM_desc_maturity$demo_fns$surv(z, params) * IPM_desc_maturity$demo_fns$grow(z1, z, params))
}

# FECUNDITY ------------------------------

kernel_fns$Fec$immature$mature <- function (z1, z, params) {
  #flowering * number of seeds * recruit svival * recruit size
  return( IPM_desc_maturity$demo_fns$flow(z, params) * IPM_desc_maturity$demo_fns$seed(z, params) * params["p.r"] * IPM_desc_maturity$demo_fns$rcsz(z1, params))
}



# KERNEL LIMIT FUNCTIONS ------------------

# lower size limit to prevent eviction
limit_lower <- function(params, IPM_desc_maturity){
  max_sd <- params["rcsz.sd"]+params["grow.sd"]
  # returns the lower size limit of the kernel
  # a value of 0.05 means that 5% of new recruits are evicted

  qnorm(0.01,mean = 0, sd = max_sd)
}

#upper size limit to prevent eviction
limit_upper <- function(params,IPM_desc_maturity){

  upper_lims <- seq(from = 2, to = 10, by = 1)

  #Calculate the proportion of individuals wrongly evicted at a size with a size limit
  fac1 <- IPM_desc_maturity$demo_fns$surv(upper_lims, params) # survival probability ~ z
  #integrate(function(x) IPM_desc_maturity$demo_fns$grow(x, z, params), U, Inf)$value

  inter <- function(x) {
    integrate(function(x) IPM_desc_maturity$demo_fns$grow(x, x, params), x, Inf)$value
  }

  fac2 <- sapply(upper_lims,inter)

  props <- fac1 * fac2
  props[props>0.01] <- 0
  upper_lim <- upper_lims[which.max(props)]

  return(upper_lim)
}

kernel_res <- function(params,IPM_desc_maturity,units_per_sd = 25){
  min_sd <- min(params["rcsz.sd"],params["grow.sd"])
  return(units_per_sd/min_sd)
}

# IPM_desc_maturity object
IPM_desc_maturity <- list(
  states = states,
  states_z = states_z,
  kernels = kernels,
  par_info = par_info,
  par_fns = par_fns,
  demo_fns = demo_fns,
  kernel_fns = kernel_fns,
  limit_lower = limit_lower,
  limit_upper = limit_upper,
  kernel_res = kernel_res
)
