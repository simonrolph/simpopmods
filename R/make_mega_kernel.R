make_mega_kernel <- function(IPM){
  # get the info needed

  stage_ref <- IPM$stage_ref
  mega_matrix <- matrix(NA,nrow = length(stage_ref), ncol = length(stage_ref))

  mega_matrices <- replicate(length(IPM$kernels),mega_matrix,simplify = F)
  names(mega_matrices) <- names(IPM$kernels)

  # loop through the nested list kernel function structure and generate kernels

  # loop through each kernel (survival, fecundity, clonal etc.)
  for (krnl in 1:length(IPM$kernels)){
    # loop through each of the DESTINATION states
    for (dst in IPM$IPM_desc$states) {
      # loop through all of the ORIGIN states
      for (ori in IPM$IPM_desc$states) {

        mega_matrices[[krnl]][stage_ref == dst,stage_ref== ori] <- IPM$kernels[[krnl]][[dst]][[ori]]

      }
    }
  }

  # make total matrix
  mega_matrices[["K"]] <- Reduce('+', mega_matrices)

  mega_matrices
}
