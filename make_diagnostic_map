make_diagnostic_map <- function(data, polygon, study_border, outcome, treat,
                                legend_title, n_breaks, plot_labs, 
                                lab_position, treat_poly, different_slope, 
                                projection){
  
  packs <- c('automap','BAMMtools', 'gstat', 'purrr', 'sf', 'sp', 'tmap')
  suppressMessages(lapply(packs, require, character.only = T))
  
  outcome <- outcome
  treat <- treat
  
  treat_poly <- treat_poly
  
  data <- data[ ,c(outcome, treat)]
  
  projection <- "+proj=longlat +ellps=WGS84 +no_defs"
  
  polygon <- polygon %>% 
    st_transform(st_crs(projection))
  
  study_border <- study_border %>% 
    st_transform(st_crs(projection))
  
  polygon_sp <- as_Spatial(polygon)
  grid <- makegrid(polygon_sp, cellsize = .003)
  grid <- SpatialPoints(grid, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))
  grid <- grid[polygon_sp, ]
  
  
  grid <- st_as_sf(grid)
  
  grid <- grid %>% 
    st_transform(st_crs("+proj=longlat +ellps=WGS84 +no_defs"))
  
  grid$treat <- as.numeric(as.numeric(st_distance(grid$geometry, treat_poly$geometry))==0,1,0)
  grid$y <- as.numeric(map_chr(grid$geometry, 1))
  grid$x <- as.numeric(map_chr(grid$geometry, 2))
  
  grid <- as_Spatial(grid)
  
  df_sp <- data %>% 
    st_transform(st_crs("+proj=longlat +ellps=WGS84 +no_defs"))
  
  treat_poly <- treat_poly %>% 
    st_transform(st_crs("+proj=longlat +ellps=WGS84 +no_defs"))
  
  colnames(df_sp)[1] <- 'outcome'
  df_sp$y <- as.numeric(map_chr(df_sp$geometry, 1))
  df_sp$x <- as.numeric(map_chr(df_sp$geometry, 2))
  
  
  df_sp <- as_Spatial(df_sp)
  
  if(different_slope == FALSE){
  
  
  vario.fit = autofitVariogram(outcome~1+x+y,
                               df_sp,
                               model = c("Exp", "Sph", "Mat"),
                               kappa = c(0, 0.01, 0.05, seq(0.2, 2, 0.1), 5, 10),
                               fix.values = c(NA),
                               start_vals = c(NA),
                               verbose = T)
  fit.v = vario.fit$var_model
  
  
  heat = krige(outcome ~ 1, locations = df_sp, beta = mean(df_sp$outcome), newdata = grid, model = fit.v)
  
  heat_2 <- st_as_sf(heat)
  heat_2 <- st_buffer(heat_2, .003, endCapStyle="SQUARE")
  
  colnames(data)[1] <- 'outcome'
  
  jenks_breaks = getJenksBreaks(data$outcome, n_breaks)
  
  
  map <- tm_shape(heat_2) + tm_fill('var1.pred', title = legend_title,
                                    breaks = dput(jenks_breaks),        
                                    palette = 'BuPu'
  ) + 
    tm_shape(data) + tm_bubbles('outcome', palette = 'BuPu', 
                                size = .1, 
                                breaks = dput(jenks_breaks), 
                                legend.size.show = F, legend.col.show = F) + 
    tm_shape(line_as_line) + tm_lines("red") + 
    tm_add_legend("line", col = 'red', lty = 1, labels = 'Study Border') + 
    tm_credits(plot_labs[1], size = 1.5, position = c(lab_position[1],lab_position[2])) + 
    tm_layout(legend.outside = T) + 
    tm_credits(plot_labs[2], size = 1.5, position = c(lab_position[3], lab_position[4]))
  
  return(map)
  }
  
  else{
    
    vario.fit = autofitVariogram(outcome ~ treat + x+y+x:y,
                                 df_sp,
                                 model = c("Exp", "Sph", "Mat"),
                                 kappa = c(0, 0.01, 0.05, seq(0.2, 2, 0.1), 5, 10),
                                 fix.values = c(NA),
                                 start_vals = c(NA),
                                 verbose = T)
    fit.v = vario.fit$var_model
    g <- gstat(NULL, 'outcome', outcome~treat*(x+y+x:y), df_sp, model = fit.v)
    
    heat = predict(g, 
                 locations = df_sp, newdata = grid, model = fit.v)
    
    heat_2 <- st_as_sf(heat)
    heat_2 <- st_buffer(heat_2, .004, endCapStyle="SQUARE")
    data <- st_as_sf(df_sp)
    #return(heat_2)
    map <- tm_shape(heat_2) + tm_fill('outcome.pred', title = legend_title,
                                      style = "cont", palette = 'BuPu'
    ) + 
      tm_shape(data) + tm_bubbles('outcome', 
                                  style = "cont",
                                 shape = 1
                                 ) + 
      tm_shape(line_as_line) + tm_lines("red") +
      tm_add_legend("line", col = 'red', lty = 1, labels = 'Study Border') + 
      tm_credits(plot_labs[1], size = 1.5, position = c(lab_position[1],lab_position[2])) + 
      tm_layout(legend.outside = T) + 
      tm_credits(plot_labs[2], size = 1.5, position = c(lab_position[3], lab_position[4]))
    
    return(map)
    
  }
}
