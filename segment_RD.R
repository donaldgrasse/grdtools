segment_RD <- function(data, outcome, run_var, segment, h, att_only){
  
  data = data %>% 
    select(outcome, run_var, segment)
  
  colnames(data)[1:3] <- c('outcome', 'run_var', 'border_point')
  
  coef_list <- list()
  pv_list <- list()
  nt_list <- list()
  nc_list <- list()
  
  
  for(i in 1:length(unique(data$border_point))){ 
    
    coef_list[i] <- rdrobust(data$outcome[data$border_point==i], 
                             data$run_var[data$border_point==i], kernel = 'tri', h = h)$coef[1]
    pv_list[i] <- rdrobust(data$outcome[data$border_point==i], 
                           data$run_var[data$border_point==i], kernel = 'tri', h = h)$pv[1]
    nt_list[i] <- rdrobust(data$outcome[data$border_point==i], 
                           data$run_var[data$border_point==i], kernel = 'tri', h = h)$N_b[1]
    nc_list[i] <- rdrobust(data$outcome[data$border_point==i], 
                           data$run_var[data$border_point==i], kernel = 'tri', h = h)$N_b[1]
    
  }
  
  result_df <- data.frame(data.frame(matrix(unlist(coef_list), nrow=length(coef_list), byrow=TRUE)), 
                          data.frame(matrix(unlist(pv_list), nrow=length(pv_list), byrow=TRUE)),
                          data.frame(matrix(unlist(nt_list), nrow=length(pv_list), byrow=TRUE)),
                          data.frame(matrix(unlist(nc_list), nrow=length(pv_list), byrow=TRUE))
  )
  colnames(result_df)[1:4] <- c('est', 'pv', 'n_c', 'n_t')
  result_df$nbin <- result_df$n_c+result_df$n_t
  result_df$ntot <- sum(result_df$nbin)
  
  result_df$att <- sum(result_df$est*(result_df$nbin)/(result_df$ntot))
  
  if(att_only == TRUE){
    
    return(result_df$att[1])
  }
  
  else{
    return(result_df)
  }
  
  
  
}
