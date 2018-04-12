# Building IPMs

params <- c(
  surv.int  =  0,
  surv.z    =   2,
  flow.int  = 0.3,
  flow.z    =   0.1,
  grow.int  =   0.2,
  grow.sd   =   0.25,
  rcsz.sd   =   0.3,
  seed.int  =   2.35,
  seed.z    =   2.37,
  p.r       =   0.4
)

context("Make a basic IPM")

# create a complete parameter set
params_c <- make_full_params(params_ic = params, IPM_desc = IPM_desc)

test_that("make_full_params completes parameter set", {
  expect_equal(params_c["grow.z"][[1]], 0.8)
  expect_equal(params_c["rcsz.int"][[1]], 0)
})

# create an IPM
IPM <- make_kernels(params_c,IPM_desc = IPM_desc)
IPM_eigen <- calc_dom_eig(MIPM=IPM,kernel = "K")

test_that("calc_dom_eig calculated the dominant eigenvalue by iteration", {
  expect_equal(IPM_eigen$lambda, Re(eigen(IPM[["megakernels"]]$K)$values[1]))
})

context("Discritisation to MPM")

# create an MPM
qtiles <- c(0,0.2,0.4,0.7,1)

MPM <- make_MPM(params_c = params_c,IPM_desc = IPM_desc,qtiles = qtiles)
MPM_eigen <- calc_dom_eig(MIPM=MPM,kernel = "K")




test_that("MPM should have same lambda value as the IPM", {
  expect_equal(MPM_eigen$lambda, IPM_eigen$lambda, tolerance=0.001)
})

test_that("MPM should have same stabe stage distribution as quantiles specified", {
  expect_equal(c(0,cumsum(MPM_eigen$w)), qtiles, tolerance=0.01)
})



###-------------------------------------


context("Make an IPM with maturity")

params <- c(
  surv.int  =  0,
  surv.z    =   2,
  flow.int  = 0.3,
  flow.z    =   0.1,
  grow.int  =   0.2,
  grow.sd   =   0.25,
  rcsz.sd   =   0.3,
  seed.int  =   2.35,
  seed.z    =   2.37,
  p.r       =   0.5,
  matu.int  =  0,
  matu.z    =   2
)


# create a complete parameter set
params_c <- make_full_params(params_ic = params, IPM_desc = IPM_desc_maturity)

test_that("make_full_params completes parameter set", {
  expect_equal(params_c["grow.z"][[1]], 0.8)
  expect_equal(params_c["rcsz.int"][[1]], 0)
})

# create an IPM
IPM <- make_kernels(params_c,IPM_desc = IPM_desc_maturity)
IPM_eigen <- calc_dom_eig(MIPM=IPM,kernel = "K")

image(t(IPM$megakernels$Fec))

test_that("calc_dom_eig calculated the dominant eigenvalue by iteration", {
  expect_equal(IPM_eigen$lambda, Re(eigen(IPM[["megakernels"]]$K)$values[1]))
})

context("Discritisation to MPM")

# create an MPM
qtiles <- c(0,0.2,0.4,0.7,1)

MPM <- make_MPM(params_c = params_c,IPM_desc = IPM_desc_maturity,qtiles = qtiles)
MPM_eigen <- calc_dom_eig(MIPM=MPM,kernel = "K")

ssd_overall <- cumsum(c(0,MPM_eigen$w[1:4]+MPM_eigen$w[5:8]))

image(t(MPM$megakernels$Fec))

test_that("MPM should have same lambda value as the IPM", {
  expect_equal(MPM_eigen$lambda, IPM_eigen$lambda, tolerance=0.001)
})

test_that("MPM should have same stabe stage distribution as quantiles specified", {
  expect_equal(ssd_overall, qtiles, tolerance=0.01)
})

###-------------------------------------

context("Make an IPM with seeds")

params <- c(
  surv.int  =  0,
  surv.z    =   2,
  flow.int  = 0.3,
  flow.z    =   0.1,
  grow.int  =   0.2,
  grow.sd   =   0.25,
  rcsz.sd   =   0.3,
  seed.int  =   2.35,
  seed.z    =   2.37,
  p.r       =   0.5,
  sdbk.surv = 0.5,
  sdbk.exit = 0.5,
  sdbk.entr = 0.5
)

# create a complete parameter set
params_c <- make_full_params(params_ic = params, IPM_desc = IPM_desc_seed)

# create an IPM
IPM <- make_kernels(params_c,IPM_desc = IPM_desc_seed)
IPM_eigen <- calc_dom_eig(MIPM=IPM,kernel = "K")

context("Make an IPM with seeds - Discritisation to MPM")

qtiles <- c(0,0.2,0.4,0.7,1)

MPM <- make_MPM(params_c = params_c,IPM_desc = IPM_desc_seed,qtiles = qtiles)
MPM_eigen <- calc_dom_eig(MIPM=MPM,kernel = "K")

MPM_eigen$lambda

test_that("MPM should have same lambda value as the IPM", {
  expect_equal(MPM_eigen$lambda, IPM_eigen$lambda, tolerance=0.001)
})

ssd_z <- MPM_eigen$w[2:5]
ssd_z <- ssd_z/sum(ssd_z)

ssd_overall <- cumsum(c(0,ssd_z))


test_that("MPM should have same stable stage distribution as quantiles specified", {
  expect_equal(ssd_overall, qtiles, tolerance=0.01)
})


# ---------------------------------------------------

context("Make an IPM descriptor template")

name = "egg"
states = c("egg","not_an_egg")
states_z = c(T,T)
kernels = c("Fec","P")
params = c("par1","par2")
n_demo_fns = 3

IPM_descriptor_template <- init_IPM_desc(name="test_IPM",states=states,states_z = states_z,kernels = kernels,params=2,n_demo_fns = 2)

