#' initialise nested lsits for making IPM kernels
#'
#' @param kernels string vector of names of kernels
#' @param params string vector of names of discrete stages
#' @return nested list
#' @export
init_nested_list <- function(kernels,states){
  list0 <- list1 <- list2 <- list()
  for (s in states){
    list0[s] <- 0
  }

  for (s in states){
    list1[[s]] <- list0
  }

  for (k in kernels){
    list2[[k]] <- list1
  }

  return(list2)
}
