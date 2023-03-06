#### About #### 
#This function `defines' the RD - it takes as arguments the data (an sf object with coordinates), 
#the polygon that defines the `action' (or treatment), the polygon that contains control units, 
#a polygon that defines the study area (most likely the union of action and control), 
#the coordinate reference system, and a boolean argument for whether fixed effects 
#should be constructed for difference border points (see Dell 2010). 


define_RD <- function(data, action_polygon, control_polygon, study_area, crs, segment_fe){
  
  #read required packages 
  packs <- c('sf',  'sp', 'dplyr', 'purrr', 'nngeo')
  suppressMessages(lapply(packs, require, character.only = T))
  
  #prep data to all be in the same 
  #crs 
  action_polygon <- action_polygon
  control_polygon <- control_polygon
  crs <- crs 
  data <- data 
  
  study_area <- study_area %>% 
    st_transform(crs = crs)
  
  action_polygon <- action_polygon %>% 
    st_transform(crs = crs)
  
  control_polygon <- control_polygon %>% 
    st_transform(crs = crs)
  
  data <- data %>% 
    st_transform(crs = crs)
  
  
  #function to make the border for the study
  #the cutoff for the running variable 
  makes_study_border <- function(action_polygon, control_polygon, crs){
    
    x_and_y <- st_intersection(action_polygon$geometry, control_polygon$geometry) %>% 
      st_transform(crs)
    
    line_as_point <- st_cast(x_and_y, "MULTIPOINT")
    line_as_line <- st_cast(line_as_point, "MULTILINESTRING")
    
    line_as_line <- line_as_line %>% 
      st_set_crs(crs)
    line_as_line <- st_cast(line_as_point, "LINESTRING")
    return(line_as_line)
    
  }
  
  line_as_line <- makes_study_border(action_polygon, control_polygon, crs)
  
  #define treatment and the running variable 
  data$treat <- as.numeric(as.numeric(st_distance(data$geometry, action_polygon$geometry))==0,1,0)
  data$dist_border <- ifelse(data$treat == 1, as.numeric(st_distance(data$geometry, line_as_line))*1, 
                             as.numeric(st_distance(data$geometry, line_as_line))*-1)
  if(segment_fe == FALSE){
    data <- st_as_sf(data)
    return(data)
    
    
  }
  
  else{
    
    x_and_y <- st_intersection(action_polygon$geometry, control_polygon$geometry) %>% 
      st_transform(crs)
    line_as_point <- st_cast(x_and_y, "POINT") %>% 
      st_set_crs(crs)
    
    point <- sf:::as_Spatial(line_as_point)
    point <- st_as_sf(point)
    point <- point %>% 
      st_set_crs(st_crs("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs"))
    
    study_area$polyID <- row.names(study_area)
    
    line_as_point_2 <- st_join(point, st_as_sf(study_area), st_within)
    line_as_point_2 <- line_as_point_2 %>% 
      filter(is.na(polyID) == F) %>% 
      distinct()    
    data$dist.segment <- as.factor(as.numeric(st_nn(data$geometry, line_as_point_2$geometry)))
    
    dummies <- model.matrix( ~ dist.segment - 1, data=data)
    data <- cbind.data.frame(data, dummies)
    data <- st_as_sf(data)
    
    return(data)
    
  }
  
}
