#' Build IPM kernels
#'
#' @param m Named numeric vector of the IPM parameters
#' @param params Named numeric vector of the IPM parameters
#' @param IPM_desc IPM descriptor object
#' @return List containing full kernel (matrix), meshpoints (vector), survival kernel (matrix), fecuntity kernel (matrix) and params (named numeric vector).
#' @seealso \code{\link{params_to_IPM}} which uses this function
#' @export
make_kernels <- function(params,IPM_desc) {

  # 1.
  # Make things about the IPM kernel
  res <- IPM_desc$kernel_res(params = params,IPM_desc = IPM_desc) # find the correct resolution to where
  L <- IPM_desc$limit_lower(params = params,IPM_desc = IPM_desc) #find the lower size limit
  U <- IPM_desc$limit_upper(params = params,IPM_desc = IPM_desc) # find the upper size limit

  m <- min(1000,max(100,ceiling((U-L)*res))) #find the dimensions of the kernel
  h <- (U - L)/m # kernel pixel width
  meshpts <- L + ((1:m) - 1/2) * h # mesh points for z

  # initialise the nested list structure for storing all the kernels
  kernels <- init_nested_list(IPM_desc$kernels,IPM_desc$states)

  # 3.
  #make the kernels from the kernel functions in the kernel_fns nested list object

  # loop through each kernel (survival, fecundity, clonal etc.)
  for (krnl in IPM_desc$kernels){
    # loop through each of the DESTINATION states
    for (dst in IPM_desc$states){
      # loop through all of the ORIGIN states
      for (ori in IPM_desc$states){

        meshpts_s <- meshpts_ts <- 0
        if (IPM_desc$states_z[IPM_desc$states == dst]){ meshpts_ts <- meshpts }
        if (IPM_desc$states_z[IPM_desc$states == ori]){ meshpts_s <- meshpts }

        # get the kernel function from the IPM descriptor
        kernel_fn <- IPM_desc$kernel_fns[[krnl]][[dst]][[ori]]

        if (is.function(kernel_fn)){
          kern <- (outer(meshpts_ts, meshpts_s,kernel_fn, params = params))

          # multiply by width of meshpts
          if (IPM_desc$states_z[IPM_desc$states == dst]){kern <- kern * h}


          # put the discritised kernel in the nested list structure
          kernels[[krnl]][[dst]][[ori]] <- kern

        } else {
          kern <- matrix(0,  length(meshpts_ts),length(meshpts_s))
        }
      }
    }
  }

  states_z <- IPM_desc$states_z
  mat_dim <- length(meshpts)
  howmany <- states_z*mat_dim
  howmany[howmany==0]<- 1
  stage_ref <- rep(IPM_desc$states,howmany)

  IPM <- list(kernels = kernels,
              meshpts = meshpts,
              params = params,
              IPM_desc = IPM_desc,
              stage_ref = stage_ref)

  # make megakernels
  megakernels <- make_mega_kernel(IPM)

  #return it all:
  return(list(megakernels = megakernels,
              kernels = kernels,
              meshpts = meshpts,
              params = params,
              IPM_desc = IPM_desc,
              stage_ref = stage_ref))

}
