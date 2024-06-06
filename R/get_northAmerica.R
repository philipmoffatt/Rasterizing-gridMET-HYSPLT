#' `get_northAmerica` provides the spatial data for the coastline necessary for delineating between the whole trajectory indices and those indices that reflect airmass characteristics overland. 
#' 
#' 

get_northAmerica <- function(){
  usa = ne_countries(country = "United States of America") 
  canada = ne_countries(country = "Canada") 
  mexico = ne_countries(country = "Mexico") 
  
  canada <- canada %>% st_as_sf
  usa <- usa %>% st_as_sf
  mexico <- mexico %>% st_as_sf
  
  north_america <- st_union(usa, canada) %>% st_union(., mexico) %>% 
    mutate(LC = factor("land")) %>% dplyr::select(LC)
}