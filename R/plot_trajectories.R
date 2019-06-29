#alpha
plot_trajectories <- function(MIPM,pop_n=100,t_steps = 50){

  #basic dataframe
  trajectories <- data.frame(t = 0,id = 0,z = 0,fec = 0,state = "",stringsAsFactors = F)


  #initialise organisms (from the fecundity kernel at ssd?)

  for (i in 1:pop_n){
    trajectories[i,] <- list(t = 0,id = i,z = seq(from=-0.25,to = 0.25,length.out = pop_n)[i],fec = 0,state = "immature")
  }

  for (t in 0:t_steps){
    if( nrow(trajectories[trajectories$t==t,])>0){

      for (id in 1:nrow(trajectories[trajectories$t==t,]) ) {
        new_row <- trajectories[trajectories$t==t,][id,]

        new_row$t <- t+1

        #find closes position on meshpts
        mesh_pos <- which(abs(MIPM$meshpts-new_row$z)==min(abs(MIPM$meshpts-new_row$z)))


        if (runif(1) < sum(MIPM$megakernels$P[,mesh_pos])){ #mortality

          # growth
          mesh_pos_next <- sample(1:dim(MIPM$megakernels$P)[1], 1, prob = MIPM$megakernels$P[,mesh_pos])
          if (mesh_pos_next > length(MIPM$meshpts)){
            new_row$state <- "mature"
            mesh_pos_next <- mesh_pos_next - length(MIPM$meshpts)
          }
          new_row$z <- MIPM$meshpts[mesh_pos_next]

          #fecundity
          fec_post_next <- 0
          if (new_row$state == "mature" && sum(MIPM$megakernels$Fec[,mesh_pos+length(MIPM$meshpts)]) != 0){
            new_row$fec <- sum(MIPM$megakernels$Fec[,mesh_pos_next+length(MIPM$meshpts)])
          } else {
            new_row$fec <- 0
          }

          #add new row to trajectories df
          row.names(new_row) <- nrow(trajectories)+1
          trajectories  <- rbind(trajectories,new_row)
        }


      }
    }

  }





  ggplot(trajectories,aes(x = t,y = z,colour = state,group = id))+
    geom_point(data = trajectories[trajectories$fec>0,],aes(size = fec),alpha = 0.3)+
    geom_line()

}
