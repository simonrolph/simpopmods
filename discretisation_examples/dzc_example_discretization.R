### IPM descriptor object


# $par_info
par_info <- data.frame(t(
  data.frame(
    ## survival
    surv.int  =  c(min = -3, max = 1, samp = T),
    surv.z    =   c(min = 0, max = 5, samp = T),
    #surv.cap = c(min = 0.01,max = 0.99),
    
    ## flowering
    flow.int  = c(min = -10, max = 1, samp = T),
    flow.z    =   c(min = 0, max = 20, samp = T),
    
    ## growth
    grow.int  =   c(min = 0, max = 1, samp = T),
    grow.z    =   c(min = 0, max = 1, samp = F),
    # constrained from 0 to 1
    grow.sd   =   c(min = 0.01, max = 0.5, samp = T),
    
    ## recruit size
    rcsz.int  =   c(min = 0, max = 0, samp = F),
    rcsz.sd   =   c(min = 0.1, max = 0.5, samp = T),
    
    ## seed size
    seed.int  =   c(min = -10, max = 100, samp = T),
    seed.z    =   c(min = 0, max = 15, samp = T),
    
    ## recruitment probability
    p.r       =   c(min = 0, max = 1, samp = T) # constrained from 0 to 1
  )
))

# $par_fns
par_fns <- list(
  rcsz.int = function(params) {
    0
  },
  grow.z = function(params) {
    1 - params["grow.int"]
  }
)

# $demo_fns
demo_fns <- list(
  ## G_z1z -  Growth function, given you are size z now returns the pdf of size z1 next time
  grow = function(z1, z, params) {
    mu <-
      params["grow.int"] + params["grow.z"] * z           # mean size next year
    sig <-
      params["grow.sd"]                                 # sd about mean
    p.den.grow <-
      dnorm(z1, mean = mu, sd = sig)             # pdf that you are size z1 given you were size z
    return(p.den.grow)
  },
  ## s_z Survival function, logistic regression
  surv = function(z, params) {
    linear.p <-
      params["surv.int"] + params["surv.z"] * z  # linear predictor
    #p <- (1/(1+exp(-linear.p))) # logistic transformation to probability with no cap
    p <-
      pmin(1 / (1 + exp(-linear.p)), rep(0.99, length(z))) # logistic transformation to probability with constant cap
    #p <- min((1/(1+exp(-linear.p))),params["surv.cap"]) # logistic transformation to probability with cut-off cap
    #p <- (1/(1+exp(-linear.p)))*params["surv.cap"] # logistic transformation to probability with scalar cap
    return(p)
  },
  ## p_bz Probability of flowering function, logistic regression
  flow = function(z, params) {
    linear.p <-
      params["flow.int"] + params["flow.z"] * z      # linear predictor
    p <-
      1 / (1 + exp(-linear.p))                                # logistic transformation to probability
    return(p)
  },
  ## b_z Seed production function
  seed = function(z, params) {
    #if (params["seed.int"] < 0) {params["seed.int"] <- params["seed.int"] * params["seed.z"]}
    N <-
      exp(params["seed.int"] + params["seed.z"] * z)   # seed production of a size z plant
    return(N)
  },
  ## c_0z1 Recruit size pdf
  rcsz = function(z1, params) {
    mu <- params["rcsz.int"]
    sig <- params["rcsz.sd"]
    p.deRecr <-
      dnorm(z1, mean = mu, sd = sig)              # pdf of a size z1 recruit
    return(p.deRecr)
  }
)

P_z1z = function (z1, z, params) {
  # survival * growth
  return(IPM_desc$demo_fns$surv(z, params) * IPM_desc$demo_fns$grow(z1, z, params))
}
## Define the fecundity kernel
F_z1z = function (z1, z, params) {
  # survival * flowering * number of seeds * recruit survival * recruit size
  return(
    IPM_desc$demo_fns$flow(z, params) * IPM_desc$demo_fns$seed(z, params) * params["p.r"] * IPM_desc$demo_fns$rcsz(z1, params)
  )
}

# IPM_desc object
IPM_desc <- list(
  par_info = par_info,
  par_fns = par_fns,
  demo_fns = demo_fns
)

params <- c(
  surv.int  =  0.1,
  surv.z    =   0.3,
  flow.int  = 0.1,
  flow.z    =   0.2,
  grow.int  =   0.27,
  grow.z    =   1 - 0.27,
  grow.sd   =   0.38,
  rcsz.int  =   0,
  rcsz.sd   =   0.2,
  seed.int  =   2,
  seed.z    =   1,
  p.r       =   0.1
)


#### DZC example
meshpts <- seq(from = -1, to = 4, by = 0.005)

dlta <- meshpts[2] - meshpts[1]
lwr <- min(meshpts) - dlta / 2
upr <- max(meshpts) + dlta / 2


P_kern <- outer(meshpts, meshpts, P_z1z,params) * dlta
F_kern <- outer(meshpts, meshpts, F_z1z,params) * dlta
K_kern <- P_kern+F_kern


# produce ssd from parameter set and IPM_desc
calc_dom_eig <- function(A, tol = 1e-8) {
  qmax <- 10 * tol
  lam <- 1
  x <- rep(1, nrow(A))
  
  while (qmax > tol) {
    x1 <- A %*% x
    qmax <- sum(abs(x1 - lam * x))
    lam <- sum(x1)
    x <- x1 / lam
  }
  
  # Having found w (within tol), get lambda
  x1 <- A %*% x
  lam <- sum(x1)
  x <- x1 / lam
  return(list(lambda = lam, w = x / sum(x)))
}

ssd <- calc_dom_eig(K_kern)$w
lambda.IPM <- calc_dom_eig(K_kern)$lambda

# approximation to the size distribition function
n_t <- approxfun(meshpts, ssd, rule = 2)
# approximation to the cumulative density function
cum_dist_n_t <- approxfun(meshpts, cumsum(ssd) / sum(ssd))

n_state_mpm <- 6 # how many classes do I want to make my MPM with?
qtiles <- seq(0, 1, length.out = n_state_mpm + 1)
qtiles <- qtiles[-c(1, length(qtiles))]

fun <- function(x, which) {
  cum_dist_n_t(x) - which
}

stage_boundaries <- qtiles
for (i in seq_along(qtiles)) {
  stage_boundaries[i] <-
    uniroot(fun, interval = range(meshpts),  which =  qtiles[i])$root
}
stage_boundaries <- c(lwr, stage_boundaries, upr)

# make a function to compute integration stuff
mk_meshpts <- function (l, u, m) {
  h <- (u - l) / m
  meshpts <- l + ((1:m) - 1 / 2) * h
  list(
    m = m,
    h = h,
    l = l,
    u = u,
    meshpts = meshpts
  )
}

# extract the IPM function we'll use
s_z <- IPM_desc$demo_fns$surv
g_z <- IPM_desc$demo_fns$grow
f_z <- F_z1z

# do the survival rates
s_stage <- numeric(n_state_mpm)
for (i in seq_len(n_state_mpm)) {
  sub_mesh <- mk_meshpts(stage_boundaries[i], stage_boundaries[i + 1], 100)
  mpts_now <- sub_mesh$meshpts
  s_stage[i] <- sum(s_z(z = mpts_now, params = params) * n_t(mpts_now)) / sum(n_t(mpts_now))
}

# do the growth transitions
g_stage <- matrix(nrow = n_state_mpm, ncol = n_state_mpm)
for (i in seq_len(n_state_mpm)) {
  sub_mesh_i <-
    mk_meshpts(stage_boundaries[i], stage_boundaries[i + 1], 250)
  mpts_now_i <- sub_mesh_i$meshpts
  for (j in seq_len(n_state_mpm)) {
    sub_mesh_j <- mk_meshpts(stage_boundaries[j], stage_boundaries[j + 1], 250)
    mpts_now_j <- sub_mesh_j$meshpts
    
    g_now <- outer(mpts_now_i, mpts_now_j, g_z, params = params)
    s_now <- s_z(mpts_now_j, params = params)
    
    g_stage[i, j] <- sum(g_now %*% (s_now * n_t(mpts_now_j))) / sum(s_now * n_t(mpts_now_j))
    g_stage[i, j] <- g_stage[i, j] * sub_mesh_i$h
  }
}

# fecundity
f_stage <- matrix(nrow = n_state_mpm, ncol = n_state_mpm)
for (i in seq_len(n_state_mpm)) {
  
  sub_mesh_i <- mk_meshpts(stage_boundaries[i], stage_boundaries[i + 1], 250)
  mpts_now_i <- sub_mesh_i$meshpts
  
  for (j in seq_len(n_state_mpm)) {
    
    sub_mesh_j <- mk_meshpts(stage_boundaries[j], stage_boundaries[j + 1], 250)
    mpts_now_j <- sub_mesh_j$meshpts
    
    f_now <- outer(mpts_now_i, mpts_now_j, f_z, params = params)
    f_stage[i, j] <- sum(f_now %*% (n_t(mpts_now_j))) / sum(n_t(mpts_now_j))
    f_stage[i, j] <- f_stage[i, j] * sub_mesh_i$h
    
  }
}

# put matrices together
P_mat <- g_stage * outer(rep(1,n_state_mpm),s_stage)
F_mat <- f_stage
K_mat <- P_mat + F_mat

# calculate lambda  
lambda.MPM <- calc_dom_eig(K_mat)$lambda
  
  
  
print(lambda.IPM) 
print(lambda.MPM)



