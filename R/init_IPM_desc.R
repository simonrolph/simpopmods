init_IPM_desc <- function(name = "untitled_IPM_descriptor",states,states_z,kernels,params,n_demo_fns){
  if (length(states) != length(states_z)){
    stop(paste0("Number of defined states (n=",length(states),") and logical vector of whether they have continious size (n=",length(states_z),") are not the same length."))
  }


  if (is.numeric(params)){
    n_params <- params
    params <- paste0("param_",1:n_params)
  } else {
    n_params <- length(params)
  }


  o <- c("")
  o <- c(o,paste0("IPM descriptor object - ",name))
  o <- c(o,"")
  o <- c(o,"")
  o <- c(o,"```{r}")
  o <- c(o,"")
  o <- c(o,paste0('states <- c(',paste(shQuote(states), collapse=", "),')'))
  o <- c(o,paste0('states_z <- c(',paste(shQuote(states_z), collapse=", "),')'))
  o <- c(o,paste0('kernels <- c(',paste(shQuote(kernels), collapse=", "),')'))
  o <- c(o,"")
  o <- c(o,"par_info <- data.frame(t(data.frame(")
  for (i in 1:n_params){
    line <- paste0("   ",params[i],'  =  c(par = "',params[i],'",min = -Inf,max = Inf,samp = T)')
    if (i != n_params) {
      line <- paste0(line,",")
    }
    o <- c(o,line)
  }
  o <- c(o,')),stringsAsFactors=F)')
  o <- c(o,"")
  o <- c(o,'par_info$min <- as.numeric(par_info$min);par_info$max <- as.numeric(par_info$max)')
  o <- c(o,"")
  o <- c(o,'par_fns <- list(')
  for (i in 1:n_params){
    line <- paste0("#    ",params[i],'  = function(params){0}')
    if (i != n_params) {
      line <- paste0(line,",")
    }
    o <- c(o,line)
  }
  o <- c(o,")")
  o <- c(o,"```")
  o <- c(o,"")
  o <- c(o,"")
  o <- c(o,"```{r}")
  o <- c(o,'demo_fns <- list(')
  for (i in 1:n_demo_fns){
    line <- paste0('  demo_fn',i,'= function(z1, z, params){')
    line <- c(line,'    #insert function here')
    line <- c(line,'    return()')

    if (i != n_demo_fns) {
      line <- c(line,'  },')
    } else {
      line <- c(line,'  }')
    }
    o <- c(o,line)
  }
  o <- c(o,")")

  o <- c(o,"```")
  o <- c(o,"")
  o <- c(o,"```{r}")
  o <- c(o,"kernel_fns <- init_nested_list(kernels,states)")
  for(kern in kernels){
    o <- c(o,"")
    o <- c(o,paste0("### KERNEL: ",kern))
    for (ori in states) {
      for (dst in states){
        line <- paste0("# ",ori," -> ",dst)
        line <- c(line,paste0("kernel_fns$",kern,"$",dst,"$",ori," <- function (z1, z, params) {"))
        line <- c(line,"   return()")
        line <- c(line,"}")
        o <- c(o,line)
        o <- c(o,"")
      }
    }
  }
  o <- c(o,"```")

  o <- c(o,"Kernel size limit functions")
  o <- c(o,"")
  o <- c(o,"```{r}")
  o <- c(o,"# lower size limit")
  o <- c(o,"limit_lower <- function(params,",name,"){")
  o <- c(o,"  ")
  o <- c(o,"  return()")
  o <- c(o,"}")
  o <- c(o,"")
  o <- c(o,"#upper size limit")
  o <- c(o,"limit_upper <- function(params,",name,"){")
  o <- c(o,"  ")
  o <- c(o,"  return()")
  o <- c(o,"}")
  o <- c(o,"")
  o <- c(o,"# function to determine the resolution of the kernel")
  o <- c(o,"kernel_res <- function(params,",name,"){")
  o <- c(o,"  ")
  o <- c(o,"  return()")
  o <- c(o,"}")
  o <- c(o,"```")
  o <- c(o,"")
  o <- c(o,"```{r}")
  o <- c(o,paste0(name,"<- list("))
  o <- c(o,"  states = states,")
  o <- c(o,"  states_z = states_z,")
  o <- c(o,"  kernels = kernels,")
  o <- c(o,"  par_info = par_info,")
  o <- c(o,"  par_fns = par_fns,")
  o <- c(o,"  demo_fns = demo_fns,")
  o <- c(o,"  kernel_fns = kernel_fns,")
  o <- c(o,"  limit_lower = limit_lower,")
  o <- c(o,"  limit_upper = limit_upper,")
  o <- c(o,"  kernel_res = kernel_res")
  o <- c(o,")")
  o <- c(o,"")
  o <- c(o,paste0("class(",name,") <- append(class(",name,"),'IPM_descriptor')"))
  o <- c(o,"")
  o <- c(o,"# Save the IPM descriptor object")
  o <- c(o,paste0("save(",name,",file = stop(\"'file' must be specified\"))"))
  o <- c(o,"```")

  return(writeLines(o))
}









