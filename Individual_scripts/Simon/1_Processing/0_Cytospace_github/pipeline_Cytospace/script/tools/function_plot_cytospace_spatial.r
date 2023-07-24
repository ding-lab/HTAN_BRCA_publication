
########################################################
# The modified SpatialPlot function
# Replace image coordinate with new cell coordinate in @misc$cytospace_coordinate$Circle or @misc$cytospace_coordinate$Sphere
########################################################
SpatialCytoPlot = function(cyto_object, 
  method = c('Circle', 'Sphere'),
  pt.size.factor = 0.6,
  stroke = NA,
  image.alpha = 0.2,
  features = NULL,
  order_by_value = F,
  stroke_color = 'gray20',
  zero_expression_color = 'gray70',
  ...
){
  # 1. Get the coordinate data
  coor_df = cyto_object@misc$cytospace_coordinate[[match.arg(method)]]

  # 2. Replace coordinate in the image
  cyto_object@images[[1]]@coordinates <- coor_df # replace the coordinate
  
  # 3. Reoder object by feature value
  if(!is.null(features) & order_by_value){
    message("Order by value, currently only support using regular ggdotplot")
    print(features)
    # Fetch datas 
    exp_df = FetchData(cyto_object, vars = features)
    print(head(exp_df))
    exp_coord_df = bind_cols(exp_df, coor_df[rownames(exp_df),])
    exp_coord_df = exp_coord_df %>% arrange(across(all_of(features)))
    # Use geom_point to plot
    ggplot(exp_coord_df, aes(x = imagecol, y = imagerow, fill = !!sym(features))) + 
      geom_point( shape = 21, size = 2 * pt.size.factor, stroke = stroke, color = stroke_color) +
      scale_fill_gradientn(colours = c(zero_expression_color, RColorBrewer::brewer.pal("YlOrRd", n = 9)[4:8])) +
      theme_void() 

  }else{

    # 3. Plot the image using SpatialPlot
    SpatialPlot(cyto_object, pt.size.factor = pt.size.factor, stroke = stroke, image.alpha=image.alpha, features=features,...)
  }
}
