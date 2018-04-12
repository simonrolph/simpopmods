#' Calculate other parameters
#'
#' @param params_ic incomplete parameter set
#' @param IPM_desc IPM descriptor object
#' @return full param set, named list
#' @export
make_full_params <- function(params_ic,IPM_desc){
  params <- params_ic
  for (param in IPM_desc$par_info$par[IPM_desc$par_info$samp == F]) {
    params[param] <- IPM_desc$par_fns[param][[1]](params_ic)
  }
  params <- params[IPM_desc$par_info$par]

  params
}
