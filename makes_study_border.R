makes_study_border <- function(action_polygon, control_polygon, crs){
  action_polygon <- action_polygon
  control_polygon <- control_polygon
  crs <- crs 
  
  x_and_y <- st_intersection(action_polygon$geometry, control_polygon$geometry) %>% 
    st_transform(crs)
  
  line_as_point <- st_cast(x_and_y, "MULTIPOINT")
  line_as_line <- st_cast(line_as_point, "MULTILINESTRING")
  
  line_as_line <- line_as_line %>% 
    st_set_crs(crs)
  line_as_line <- st_cast(line_as_point, "LINESTRING")
  return(line_as_line)
  
}
