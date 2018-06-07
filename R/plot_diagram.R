
plot_diagram <- function(IPM_desc){
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  IPM_desc_name <- deparse(substitute(IPM_desc))


  krnl_colours <- data.frame(krnl = unique(IPM_desc$kernels))
  krnl_colours$col <- gg_color_hue(nrow(krnl_colours))

  MIPM_gra <- create_graph()

  # add nodes for each life stage
  for(i in 1:length(IPM_desc$states)){
    if (IPM_desc$states_z[i]){
      MIPM_gra <- add_node(MIPM_gra,
                           label = IPM_desc$states[i],
                           node_aes = node_aes(fillcolor = "grey",
                                               fontcolor = "black"))
    } else {
      MIPM_gra <- add_node(MIPM_gra,
                           label = IPM_desc$states[i],
                           node_aes = node_aes(fillcolor = "grey",
                                               fontcolor = "black",
                                               shape = "diamond"))
    }
  }

  for (krnl in IPM_desc$kernels){
    # loop through each of the DESTINATION states
    for (dst in IPM_desc$states){
      # loop through all of the ORIGIN states
      for (ori in IPM_desc$states){

        kernel_fn <- IPM_desc$kernel_fns[[krnl]][[dst]][[ori]]

        if (is.function(kernel_fn)){
          labe <- paste(deparse(kernel_fn),collapse = "")
          labe <- gsub("\"","",labe)
          labe <- gsub("\\{","START",labe)
          labe <- gsub("\\}","END",labe)

          labe <- sub(".*START *(.*?) *END.*", "\\1", labe)

          labe <- gsub("return\\(","",labe)
          labe <- gsub("params\\[","",labe)
          labe <- gsub(paste0(IPM_desc_name,"\\$demo_fns\\$"),"",labe)
          labe <- gsub(" params","",labe)
          labe <- gsub("\\]","",labe)
          labe <- substr(labe, 1, nchar(labe)-1)
          labe <-gsub("\\s+"," ",labe)
          labe <-gsub(",\\)","\\)",labe)
          labe <-gsub(", \\)","\\)",labe)

          MIPM_gra <- add_node(MIPM_gra,
                               label = labe,
                               node_aes = node_aes(shape = "rectangle",
                                                   fontsize = 5,
                                                   height = 0.3,
                                                   width = 2,
                                                   fillcolor = krnl_colours[krnl_colours$krnl == krnl,"col"],
                                                   fontcolor = "black")
          )



          MIPM_gra <- add_edge(MIPM_gra,
            from = ori,
            to = labe,
            rel = krnl,
            edge_aes = edge_aes(color = krnl_colours[krnl_colours$krnl == krnl,"col"],arrowhead = "none",fontcolor = "black")

            )

          MIPM_gra <- add_edge(MIPM_gra,
                               from = labe,
                               to = dst,
                               rel = krnl,
                               edge_aes = edge_aes(color = krnl_colours[krnl_colours$krnl == krnl,"col"],fontcolor = "black")

          )



        }


      }
    }
  }

  MIPM_gra
}

