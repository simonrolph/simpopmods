

get_IPM_desc <- function(){
  get(paste0(gsub( "@.*$", "",sys.call(-1))))
}

get_dm_fns <- function(){
  print(sys.function(1))
  get(gsub( "@.*$", "",sys.call(-1)))@demo_fns
}
