####################################################################################
# Functions to calculate spatial distribution of cells on the image
####################################################################################
## API
## New method of Image coordinate method
CalculateSpatialDistribution <- function(cyto_object, method = c('all','Circle','Sphere'), radius = 70) {
    # 0. extract relevant objects from inputs
    loc_use <- cyto_object@misc$cytospace_result$assigned_locations # Location files 
    # st_img_use <- st_use@images[[1]] # Image
    cyto_coord <- cyto_object@images[[1]]@coordinates # Image Coordinate

    # 1. Get coordinate table with imagerow and image col
    coor_df_new_all <- cyto_coord %>% 
        rownames_to_column('UniqueCID') %>%
        left_join(
            loc_use, by = c('UniqueCID', 'row', 'col')
        ) 

    #@@ 4 Rearrange spots on the image
    # A. Get Cell Count information Per Spot 
    # 4.1. calcultate number of cell per spot
    ncount_spot = coor_df_new_all %>% count(SpotID) 
    # 4.2. Get center location
    distinct_df = coor_df_new_all[,c('SpotID','imagerow','imagecol')] %>% distinct
    # 4.3 Combine
    spot_df = left_join(ncount_spot, distinct_df, by = 'SpotID')

    # 4.4 Get New Locations
    # 4.5. Replace image coordinate with new cell coordinate
    method = match.arg(method)
    # 4.5.1. Circle
    if(method == 'Circle' | method == 'all'){
      message(str_c('Calculating new spot location using circle_divide_points method, radius = ', radius))
      subspot_circle_df = CalculatePointCoordinate(spot_df, radius = radius, method = 'circle_divide_points')
      coor_df_new_clean_circle = cbind(
          subspot_circle_df %>% arrange(SpotID),
          coor_df_new_all %>% arrange(SpotID)
      ) %>% 
          column_to_rownames('UniqueCID') %>%
          # Clean up columns
          select(c('tissue','row','col','imagerow_new','imagecol_new')) %>% 
          rename(imagerow = imagerow_new, imagecol = imagecol_new)
      
      # Add result to misc slot
      cyto_object@misc$cytospace_coordinate$Circle = coor_df_new_clean_circle
      cyto_object@misc$cytospace_coordinate$Circle_radius = radius
    }
    
    # 4.5.2 Sphere
    if(method == 'Sphere' | method == 'all'){
      message(str_c('Calculating new spot location using generatePointsInCircle method, radius = ', radius))
      subspot_sphere_df = CalculatePointCoordinate(spot_df, radius = radius, method = 'generatePointsInCircle')
      coor_df_new_clean_sphere = cbind(
          subspot_sphere_df %>% arrange(SpotID),
          coor_df_new_all %>% arrange(SpotID)
      ) %>% 
          column_to_rownames('UniqueCID') %>% 
          # Clean up columns
          select(c('tissue','row','col','imagerow_new','imagecol_new')) %>% 
          rename(imagerow = imagerow_new, imagecol = imagecol_new)
      
      # Add result to misc slot
      cyto_object@misc$cytospace_coordinate$Sphere = coor_df_new_clean_sphere
      cyto_object@misc$cytospace_coordinate$Sphere_radius = radius
    }
    
    
    return(cyto_object)
}

####################################################################################
### Helper function
# Spot df need to have columns: SpotID, imagerow, imagecol, n
CalculatePointCoordinate = function(spot_df, radius = 70, method = c('circle_divide_points', 'generatePointsInCircle')){
  calculate_func = if(match.arg(method) == 'circle_divide_points'){
    circle_divide_points
  } else if(match.arg(method) == 'generatePointsInCircle'){
    generatePointsInCircle
  }
  pmap(spot_df, function(SpotID, imagerow, imagecol, n){
        df = calculate_func(
          center_x = imagerow, 
          center_y = imagecol, 
          radius = radius, 
          n = n # number of points
        ) %>%
          as.data.frame %>% 
          setNames(c('imagerow_new', 'imagecol_new')) %>%
          mutate(SpotID = SpotID)
        return(df)
    }) %>% bind_rows()
}


########################################################
# Funcion to calculate spatial distributions
# Both with help of ChatGPT
########################################################
## A. divide a circle into n parts
circle_divide_points <- function(center_x, center_y, radius, n, start_angle = 0) {
  theta <- seq(start_angle, start_angle+2*pi, length.out=n+1)[-1] # angles between points
  x <- center_x + radius*cos(theta) # x-coordinates of points
  y <- center_y + radius*sin(theta) # y-coordinates of points
  coords <- cbind(x,y) # combine x and y coordinates into a matrix
  return(coords)
}

## B. Generate non-overlap points in a circle with a minimum distance between them
generatePointsInCircle <- function(n, radius, center_x=0, center_y=0, min_dist_ratio = 3, start_from_center = TRUE) {
  # Currently set d to be 1/10 of the radius
  # d is the minimum distance between the centers of two spheres/points
  d <- radius/min_dist_ratio

  # Calculate the maximum number of iterations
  max_iter <- 10000
  
  # Generate the first point at a random position within the circle
  points <- matrix(0, nrow=n, ncol=2)
  points[1,] <- if(start_from_center) c(0,0) else runif(2, -radius, radius)
  
  # Check if n > 2, if not return result
  if (n == 1) {
    return(as.data.frame(points + c(center_x, center_y)))
  }

  # Generate the remaining points
  for (i in 2:n) {
    # Initialize the distance to the closest neighbor to be greater than the sum of the sphere diameters
    dist_min <- d
    
    # Perform a Monte Carlo simulation to find a valid position for the new sphere
    iter <- 0
    while (iter < max_iter) {
      # Generate a random position within the circle
      x <- runif(1, -radius, radius)
      y <- runif(1, -radius, radius)
      point <- c(x, y)
      
      # Check if the new sphere overlaps with any of the existing spheres
      overlap <- FALSE
      for (j in 1:(i-1)) {
        dist <- sqrt((point[1] - points[j,1])^2 + (point[2] - points[j,2])^2)
        if (dist < dist_min) {
          overlap <- TRUE
          break
        }
      }
      
      # If the new sphere does not overlap with any existing spheres, add it to the list of points and break out of the loop
      if (!overlap) {
        points[i,] <- point
        break
      }
      
      # If the maximum number of iterations is reached, print a warning message and 
      # New: Just return the current random position
      if (iter == max_iter - 1) {
        warning("Warning: Maximum number of iterations reached. Could not generate valid configuration.")
        warning("Still assign as the random position")
	points[i, ] <- point
	break
      }
      
      # If the new sphere overlaps with an existing sphere, update the distance to the closest neighbor and continue the simulation
      #dist_min <- dist - d
      iter <- iter + 1
    }
  }
  
  # Shift the points so that the center of the circle is at (center_x, center_y)
  #points <- points + c(center_x, center_y)
  points <- as.data.frame(points) %>% setNames(c('x', 'y'))
  points$x <- points$x + center_x
  points$y <- points$y + center_y

  # Return the list of points
  return(points)
}

########################################################
# TESTING SCRIPT 
########################################################

# ## check if generatePointsInCircle is working fine
# ## test case : a uniform box 1:10 x 1:10
# ## 1. generate 100 points in box, x and y should be between 1 and 10
# box_array = expand.grid(1:3,1:3) %>% setNames(c('x','y'))
# ## 2. Run generatePointsInCircle on all points in boxr = 1, n_points = 20
# points_df = pmap(box_array, function(x,y, ...){
#     generatePointsInCircle(center_x = x,center_y = y, radius = 0.1, n = 10)
# }) %>% bind_rows %>% setNames(c('x','y'))

# # Plot
# pdf('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/PKD_local/Batch1-4/5_Cytospace/5_createCytospaceSeurat/out/STDim_TEST.pdf', w=12, h =6)
# ggplot() + 
#     geom_point(data = box_array, aes(x = x, y = y), color = 'black', size = 0.5) +
#     geom_point(data = points_df, aes(x = x, y = y), color = 'red', size = 0.5) +
#     theme_bw() +
#     theme(legend.position = 'none')
# dev.off()

