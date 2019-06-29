make_MPM <- function(params_c, IPM_desc, qtiles, submesh_res = 200) {
  # CHECKS
  if (qtiles[1] != 0) {
    stop("First value of quantiles must be 0")
  }
  if (qtiles[length(qtiles)] != 1) {
    stop("Last value of quantiles must be 1")
  }
  if (!all(qtiles == cummax(qtiles))) {
    stop("Quantiles must be monotonically increasing")
  }

  # FUNCTIONALITY

  # make the IPM kernel to get the SSD for finding quantiles
  IPM <- make_kernels(params_c, IPM_desc)
  meshpts <- IPM$meshpts
  dlta <- meshpts[2] - meshpts[1]

  # produce ssd from parameter set and IPM_desc
  eigen_stuff <- calc_dom_eig(IPM)
  IPM$ssd <- ssd <- eigen_stuff$w

  n_t <- function(state, z, ssd, stage_ref) {
    if (sum(stage_ref == state) == 1) {
      return(rep(1, length(z)))
    } else {
      ssd <- ssd[stage_ref == state]
      ssd_fun <- approxfun(meshpts, ssd, rule = 2)
      return(ssd_fun(z))
    }
  }

  # average ssd across all continious states

  # get only the SSD from stages with continious size
  states_with_z <- IPM$IPM_desc$states[IPM$IPM_desc$states_z]
  ssd_avg <- ssd[IPM$stage_ref %in% states_with_z]

  ssd_avg2 <- rep(0, length(ssd_avg) / length(states_with_z))

  for (i in states_with_z) {
    ssd_avg2 <- ssd_avg2 + ssd[IPM$stage_ref == i]
  }

  ssd <- ssd_avg2
  ssd <- ssd / sum(ssd)

  # approximation to the cumulative density function
  cum_dist_n_t <- approxfun(meshpts, cumsum(ssd) / sum(ssd))

  # how many classes do I want to make my MPM with?
  n_state_mpm <- length(qtiles) - 1

  # little function
  fun <- function(x, which) {
    cum_dist_n_t(x) - which
  }

  #remove the first and last number from the quantiles
  qtiles <- qtiles[-c(1, length(qtiles))]

  #find stage boundaries
  stage_boundaries <- qtiles
  for (i in seq_along(qtiles)) {
    stage_boundaries[i] <-
      uniroot(fun, interval = range(meshpts),  which =  qtiles[i])$root
  }

  # make some useful objects/variables about the stage boundaries
  lwr <- min(meshpts) - dlta / 2
  upr <- max(meshpts) + dlta / 2
  stage_boundaries <- c(lwr, stage_boundaries, upr)

  #other stuff about the class boundaries
  llims <- stage_boundaries[1:n_state_mpm]
  ulims <- stage_boundaries[2:(n_state_mpm + 1)]
  bin_sizes <- ulims - llims
  mids <- llims + bin_sizes / 2

  # function to create a grid to integrate over
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

  # initialise a nested list to store the final matrices
  kernels <- init_nested_list(IPM_desc$kernels, IPM_desc$states)

  # loop through each kernel (survival, fecundity, clonal etc.)
  for (krnl in IPM_desc$kernels) {
    # loop through each of the DESTINATION states
    for (dst in IPM_desc$states) {
      # loop through all of the ORIGIN states
      for (ori in IPM_desc$states) {
        nrows <- ncols <- 1

        if (dst %in% states_with_z) {
          nrows <- n_state_mpm
        }

        if (ori %in% states_with_z) {
          ncols <- n_state_mpm
        }

        matrix_values <- matrix(nrow = nrows, ncol = ncols)

        # get the kernel function from the IPM descriptor
        kernel_fn <- IPM_desc$kernel_fns[[krnl]][[dst]][[ori]]

        # if we actually are dealing with a function
        if (is.function(kernel_fn)) {
          # loop through matrix i
          for (i in seq_len(nrows)) {
            if (dst %in% states_with_z) {
              # create the submesh to integrate over
              sub_mesh_i <-
                mk_meshpts(stage_boundaries[i],
                           stage_boundaries[i + 1],
                           submesh_res)
              mpts_now_i <- sub_mesh_i$meshpts
              h <- sub_mesh_i$h
            } else {
              mpts_now_i <- 1
              h <- 1
            }

            # loop through matrix j
            for (j in seq_len(ncols)) {
              if (ori %in% states_with_z) {
                # create submesh in 2nd dimensions
                sub_mesh_j <-
                  mk_meshpts(stage_boundaries[j],
                             stage_boundaries[j + 1],
                             submesh_res)
                mpts_now_j <- sub_mesh_j$meshpts
              } else {
                mpts_now_j <- 1
              }

              # integrate
              fn_now <-
                outer(mpts_now_i, mpts_now_j, kernel_fn, params = params_c)

              # ssd
              n_t_now <- n_t(ori, mpts_now_j, IPM$ssd, IPM$stage_ref)

              matrix_values[i, j] <-
                sum(fn_now %*% n_t_now) / sum(n_t_now) * h
            }
          }



          # put the discritised kernel in the nested list structure
          kernels[[krnl]][[dst]][[ori]] <- matrix_values
        }

      }
    }


  }

  # make the stage ref
  states_z <- IPM_desc$states_z
  mat_dim <- n_state_mpm
  howmany <- states_z * mat_dim
  howmany[howmany == 0] <- 1
  stage_ref <- rep(IPM_desc$states, howmany)

  # Make the MPM object
  MPM <-
    list(
      kernels = kernels,
      meshpts = mids,
      params = params_c,
      IPM_desc = IPM_desc,
      stage_ref = stage_ref
    )

  # make megakernels
  megakernels <- make_mega_kernel(MPM)

  # return all the bits and bobs
  return(
    list(
      megakernels = megakernels,
      kernels = kernels,
      meshpts = mids,
      stage_boundaries = stage_boundaries,
      params = params_c,
      IPM_desc = IPM_desc,
      stage_ref = stage_ref
    )
  )

}
