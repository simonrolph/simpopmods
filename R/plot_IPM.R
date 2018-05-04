plot_IPM <- function(MIPM){
  par(mfrow=(c(1,1)))
  plot.new()
  IPM_desc <- MIPM$IPM_desc

  par(mfrow=(c(length(IPM_desc$states),length(IPM_desc$states))))


  for (krnl in IPM_desc$kernels) {
    title(krnl, side = 3, line = -21, outer = TRUE)


    # loop through each of the DESTINATION states
    for (dst in IPM_desc$states) {
      # loop through all of the ORIGIN states
      for (ori in IPM_desc$states) {

        if (sum(MIPM$kernels[[krnl]][[dst]][[ori]]) != 0){
          image(x = MIPM$meshpts,
                y = MIPM$meshpts,
                z = t(MIPM$kernels[[krnl]][[dst]][[ori]]),
                main = paste(krnl,ori,dst,sep="-"),
                xlab = "size",
                ylab = "size next")
        } else {
          frame()
        }
      }
    }
  }





  par(mfrow=(c(1,1)))
}
