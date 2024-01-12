noise_test <- function(data, outcome, run_var, nsims, nonparametric, plot_title){

  packs <- c('sf', 'sp', 'dplyr', 'purrr', 'rdrobust', 'ggplot2', 
             'ggthemes')
  suppressMessages(lapply(packs, require, character.only = T))
  
  nsims <- nsims 
  
  df_spat <- st_as_sf(data) %>% 
    st_transform(st_crs("+proj=longlat +datum=WGS84"))
  
  y <- st_coordinates(df_spat$geometry)[,1]
  x <- st_coordinates(df_spat$geometry)[,2]
  xy <- as.data.frame(cbind(x, y))
  colnames(xy) <- c('x', 'y')
  
  
  df_spat$ycoor <- st_coordinates(df_spat$geometry)[,1]
  df_spat$xcoor <- st_coordinates(df_spat$geometry)[,2]
  
  df_spat <- cbind.data.frame("outcome" = outcome, "run_var" = run_var, df_spat)
  
  df_spat <- as_Spatial(st_as_sf(df_spat))
  
  set.seed(168637)

  v = gstat::variogram(outcome ~ 1 + ycoor + xcoor, data = df_spat)
  v.fit = fit.variogram(v, vgm(model = 'Mat'), fit.kappa = TRUE)

  g.dummy <- gstat::gstat(formula= z ~ 1, locations=~x+y, dummy=TRUE, beta=mean(outcome), 
                   model=vgm(range=v.fit[,3],psill=v.fit[,2], model='Mat', maxdist =v.fit[,3]))
  library(gstat)
  yy <- predict(g.dummy, newdata=xy, nsim=nsims)

  colnames(yy)[1:2] <- c('xcoor','ycoor')
  df_spat <- st_as_sf(df_spat) %>% 
    st_transform(st_crs("+proj=longlat +datum=WGS84"))
  
  df_spat$geometry <- NULL 
  df_spat <- df_spat[ ,c('run_var', 'xcoor','ycoor', "outcome")]
  df_spat <- merge(yy, df_spat, by = c('xcoor', 'ycoor'))
  
  t_list <- as.data.frame(matrix(NA, ncol = 1, nrow = nsims))
  
  if(nonparametric == FALSE){
    
    bw <- rdbwselect(df_spat$outcome, df_spat$run_var, kernel = 'uni')$bws[1]
    df_spat <- df_spat %>% 
      filter(abs(run_var)<= bw)

  for(i in 1:nsims){
    t_list[i,] = lm(df_spat[,i] ~ I(as.numeric(run_var>=0, 1,0))*run_var,
                                               data = df_spat)$statistic[2]
      }
  t_list$sig_emp <- as.numeric(abs(t_list$V1) >=  1.96, 1,0)
  share_rej <- round(sum(t_list$sig_emp)/nsims, 2)
 
  real_t <- lm(df_spat$outcome ~ I(as.numeric(run_var>=0, 1,0))*run_var,
               data = df_spat)$statistic[2]
  
  p1 <- ggplot(t_list, aes(x=V1)) + 
    geom_histogram() + 
    geom_vline(data=t_list, aes(xintercept=real_t), col = 'black') +
    geom_vline(data=t_list, aes(xintercept=-1.96), linetype = "dashed", label = 'test', col = 'red') +
    geom_vline(data=t_list, aes(xintercept=1.96), linetype = "dashed", col = 'red') +
    labs(x = "Z-statistic", y = 'Frequency', title =  plot_title) + 
    annotate("label", x = -Inf, y = Inf, label = paste('Share Rejected:', share_rej), 
             hjust=0, vjust=1, family = "Times") + 
    theme_tufte()
  return(p1)
  }
  
  
  else{
 
    real_z <- rdrobust(df_spat$outcome, df_spat$run_var)$z[3]
    
    for(i in 1:nsims){
      t_list[i,] = rdrobust(df_spat[,i], df_spat$run_var)$z[3]
    }
    
    t_list$sig_emp <- as.numeric(abs(t_list$V1) >=  abs(real_z), 1,0)
    share_rej <- round(sum(t_list$sig_emp)/nsims, 2)
    
    p1 <- ggplot(t_list, aes(x=V1)) + 
      geom_histogram() + 
      geom_vline(data=t_list, aes(xintercept=real_z), col = 'black') +
      geom_vline(data=t_list, aes(xintercept=-1.96), linetype = "dashed", label = 'test', col = 'red') +
      geom_vline(data=t_list, aes(xintercept=1.96), linetype = "dashed", col = 'red') +
      labs(x = "Z-statistic", y = 'Frequency', title = plot_title) + 
      annotate("label", x = -Inf, y = Inf, label = paste('Share Rejected:', share_rej), 
               hjust=0, vjust=1, family = "Times") + 
      theme_tufte()
    return(p1)
  }
}
