# THIS SCRIPT CONTAINS THE FUNCTIONS THAT ARE USED TO ANALYSE MOBILE BICYCLE MEASUREMENTS AND TO CREATE TEMPERATURE CHARTS BASED ON THESE MEASUREMENTS

#--------------------------------------------------------------------------INITIALIZATION FUNCTION-------------------------------------------------------------------------------#

# Function that checks if the names in a vector are appearing in a dataframe, if not the program stops
names_check = function(dataframe, name_vector){
  
  number_wrong_names = 0
  for (name in name_vector){
    if (!(name %in% dataframe$Name)){
      print(paste0(name, " is not a valid station"))
      number_wrong_names = number_wrong_names + 1
    }
  }
  if (number_wrong_names != 0){
    stop(call. = FALSE)
  }
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------READ OUT FUNCTION---------------------------------------------------------------------------------#

# Function that creates a dataframe from the measured data stored in a csv-file and that converts the coordinates from (D)DDMM.MMMM format to numerical format
# The coordinates are also transformed from WGS84 (= GPS coordinate system) to Lambert-72 system (used in the BBK land cover map of Flanders)
reading_data = function(datafile){
  
  # Read the data from a csv-datafile and give every data point a name (number)
  dataframe_measurements = fread(datafile, sep=",", dec=".")
  
  names(dataframe_measurements) = c("date", "time", "Y_WGS_GPS", "X_WGS_GPS", "temperature", "pressure", "humidity")
  dataframe_measurements$Name = seq(1, NROW(dataframe_measurements), 1)
  
  # Convert the (D)DDMM.MMMM format to numerical format
  degreeslist_Y = as.numeric(substr(dataframe_measurements$Y_WGS_GPS, 1, 2))
  degreeslist_X = as.numeric(substr(dataframe_measurements$X_WGS_GPS, 1, 1))
  minutelist_Y= as.numeric(substr(dataframe_measurements$Y_WGS_GPS, 3, 8))
  minutelist_X = as.numeric(substr(dataframe_measurements$X_WGS_GPS, 2, 7))
  dataframe_measurements$Y_WGS_GPS = degreeslist_Y + (minutelist_Y/60)
  dataframe_measurements$X_WGS_GPS = degreeslist_X + (minutelist_X/60)
  
  # Create a column with date + time
  dataframe_measurements$fulldate = with(dataframe_measurements, paste(dataframe_measurements$date, dataframe_measurements$time))
  dataframe_measurements$fulldate = as.POSIXct(dataframe_measurements$fulldate, format= "%d/%m/%Y %H:%M:%S", tz=Sys.timezone(), origin)
  
  # Transform the coordinates from WGS84 system to Lambert-72 system (EPSG:31370)
  dataframe_measurements = coordinate_projection_lambert(dataframe_measurements)
  return(dataframe_measurements)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------COORDINATE FUNCTIONS-------------------------------------------------------------------------------#

# Function that transforms the coordinates from WGS84 (= GPS coordinate system) to Lambert-72 system (used in the BBK land cover map of Flanders)
coordinate_projection_lambert = function(dataframe){
  
  # Setting existing coordinate as lat-long system
  cord.dec = SpatialPoints(cbind(dataframe$X_WGS_GPS, dataframe$Y_WGS_GPS), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  # Transform the coordinates to Lambert_72 system using EPSG:31370
  cord.UTM = spTransform(cord.dec, CRS("+init=epsg:31370"))
  
  # Extract the Lambert-72 coordinates
  dataframe$X_lambert = cord.UTM$coords.x1
  dataframe$Y_lambert = cord.UTM$coords.x2
  
  return(dataframe)
}

# Function that transforms a dataframe of WGS84 coordinates to a list of Lambert-72 coordinates (used in the BBK land cover map of Flanders)
WGS_to_lambert = function(dataframe_WGS){
  
  # Transform the WGS84 coordinates in the dataframe to Lambert-72 coordinates 
  dataframe_lambert = coordinate_projection_lambert(dataframe_WGS)
  
  # Initialize a list to store the Lambert-72 coordinates of the dataframe
  list_lambert_coords = vector(mode = "list", length = length(dataframe_lambert$X_lambert))
  
  # Fill the list with the Lambert-72 coordinates of the dataframe
  for (row in 1:length(dataframe_lambert$X_lambert)){
    list_lambert_coords[[row]] = c(dataframe_lambert$X_lambert[row], dataframe_lambert$Y_lambert[row])
  }
  return(list_lambert_coords)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------LAND COVER FUNCTIONS-------------------------------------------------------------------------------#

# Function that creates circular buffers for different distances (based on the BBK land cover map of Flanders)
creating_buffers_BBK = function(dataframe_measurements, distance_vector){
  
  # Create a geometry with the right projection (EPSG:31370)
  points.sf.lambert = st_as_sf(dataframe_measurements, coords = c("X_lambert", "Y_lambert"), crs = 31370, agr = "constant")
  
  # Create buffers for each distance in the vector of distances specified in the general settings
  for (index in 1:length(distance_vector)) {
    
    # Create the buffer and add a column with the distance of the buffer
    buffer = st_buffer(points.sf.lambert, distance_vector[index], nQuadSegs = 200)
    buffer$distance = distance_vector[index]
    
    # Add Lambert-72 coordinates to the buffer 
    buffer$X = dataframe_measurements$X_lambert
    buffer$Y = dataframe_measurements$Y_lambert	
    
    # Create a dataframe with all buffers
    if (index == 1) {
      buffers = buffer
    } 
    else {
      buffers = rbind(buffers, buffer)
    }
    print(paste0("Buffers created for radius: ", distance_vector[index], " m (based on BBK land cover map)"))
  }
  return(buffers)
}

# Function that creates a small land cover map based on the BBK land cover map of Flanders (used to calculate the land cover fractions)
creating_small_landcover_map_BBK = function(BBK_path, buffer_dataframe, distance_vector){
  
  # Create a RasterLayer object from the BBK land cover map
  BBK_map = raster(BBK_path)   
  
  # Define the limits of the small land cover map
  X_minimum = min(buffer_dataframe$X) - max(distance_vector) - 100.
  X_maximum = max(buffer_dataframe$X) + max(distance_vector) + 100.
  Y_minimum = min(buffer_dataframe$Y) - max(distance_vector) - 100.
  Y_maximum = max(buffer_dataframe$Y) + max(distance_vector) + 100.

  # Calculate the time it takes for the program to create the small land cover map
  start_time = Sys.time()
  chart = as(extent(X_minimum, X_maximum, Y_minimum, Y_maximum), 'SpatialPolygons')
  crs(chart) = "+init=epsg:31370"
  small_landcover_map = crop(BBK_map, chart)
  end_time = Sys.time()
  duration = end_time-start_time
  print(paste0("Time to make the BBK small land cover map: ", duration, " seconds"))
  return(small_landcover_map)
}

# Function that creates a small land cover map based on the BBK land cover map of Flanders (used to calculate the land cover fractions)
# This will be used to create the visualisation of the detailed maps around 'Rivierenhof'
creating_small_landcover_map_BBK_detailed = function(BBK_path, buffer_dataframe, distance_vector){
  
  # Create a RasterLayer object from the BBK land cover map
  BBK_map = raster(BBK_path)   
  
  # Define the limits of the small land cover map
  X_minimum = min(buffer_dataframe$X) - max(distance_vector) - 500.
  X_maximum = max(buffer_dataframe$X) + max(distance_vector) + 500.
  Y_minimum = min(buffer_dataframe$Y) - max(distance_vector) - 100.
  Y_maximum = max(buffer_dataframe$Y) + max(distance_vector) + 100.
  
  # Calculate the time it takes for the program to create the small land cover map
  start_time = Sys.time()
  chart = as(extent(X_minimum, X_maximum, Y_minimum, Y_maximum), 'SpatialPolygons')
  crs(chart) = "+init=epsg:31370"
  small_landcover_map = crop(BBK_map, chart)
  end_time = Sys.time()
  duration = end_time-start_time
  print(paste0("Time to make the BBK small land cover map: ", duration, " seconds"))
  return(small_landcover_map)
}

# Function that changes the resolution of a RasterLayer object by a specified factor 
# This is done by splitting the data in blocks for memory consuming reasons
aggregate_blocks = function(raster_object, res_factor, fun=mean, filename=''){
  
  # Split the input raster in blocks
  blocks_data = blockSize(raster_object)                   # Find automatic blocks
  blocks_data$startrow = seq(1, nrow(raster_object), res_factor)        # Calculate the starting row of each block
  blocks_data$nrows = diff(c(blocks_data$startrow,(nrow(raster_object)+1)))   # Calculate the number of rows of each block
  blocks_data$endrow = blocks_data$startrow + blocks_data$nrows - 1          # Calculate the ending row of each block
  blocks_data$n = length(blocks_data$startrow)                 # number of blocks
  
  # Write the output raster
  raster_out = raster_object                            # Get raster basic structure
  res(raster_out) = res(raster_out)*res_factor                 # Set the resolution of the output raster
  raster_out = writeStart(raster_out, filename=filename)
  
  # Split the output raster in blocks
  blocks_data_out = blockSize(raster_out)                          # Find automatic blocks
  blocks_data_out$row = seq(1, nrow(raster_out), by=1)  
  
  # Loop the blocks defined above for aggregating
  for (j in 1:nrow(raster_out)) {
    
    # Crop the input raster using the block
    ext = extent(raster_object, r1=blocks_data$startrow[j], r2=blocks_data$endrow[j], c1=1, c2=ncol(raster_object))
    sub_raster = crop(raster_object, ext)
    
    # Aggregate pixel values for the cropped input raster
    agg_values = aggregate(sub_raster, res_factor, fun=fun)
    agg_values = getValues(agg_values)
    
    # Write the values of the pixel aggregation to the output raster
    writeValues(raster_out, agg_values, blocks_data_out$row[j])
    
  }
  raster_out = writeStop(raster_out)           # Stop writing to output raster
}

# Function that calculates the land cover fractions
# The land cover fractions of fourteen land cover types are calculated based on the BBK land cover map of Flanders
calculate_landcover_BBK = function(buffer_dataframe, distance_vector, small_landcover_map, cores){
  
  # Create a list of the observation numbers
  pointlist = unique(buffer_dataframe$Name)

  # Sort the vector of buffer distances in increasing order
  distance_vector = sort(distance_vector, decreasing = FALSE)
  
  # Make a cluster of cores used to calculate the land cover fractions
  cl = makeCluster(cores)
  clusterExport(cl=cl, list("small_landcover_map", "buffer_dataframe", "distance_vector", "pointlist"), envir = environment())
  
  # Calculate the land cover fractions for each data point (for all buffers)
  landcoverlist = foreach (j = 1:length(pointlist), .export=c("small_landcover_map", "buffer_dataframe", "distance_vector", "pointlist"), .combine = rbind, .packages = c("raster", "fasterize", "sf")) %dopar% {
    datapoint_number = pointlist[j]
    
    # Collection of all buffers belonging to a certain data point
    buffer_collection = buffer_dataframe[which(buffer_dataframe$Name == datapoint_number),] 
    classtat.results = data.frame()
    
    # Crop the small land cover map to the extent of the largest buffer of the collection of all buffers (largest buffer map)
    biggest_geometry_point = buffer_collection[which(buffer_collection$distance == max(buffer_collection$distance)),]
    largest_buffer_map = raster::crop(small_landcover_map, extent(biggest_geometry_point))
    
    # Initialize a dataframe for the land cover fractions of the data point
    landcover_datapoint = data.frame(impervious=numeric(), green=numeric(), water=numeric())
    
    # Calculate the land cover fractions of all buffers belonging to the data point
    for (distance in distance_vector){
      
      # Only for the smallest buffer of the collection of all buffers
      if (distance == min(distance_vector)){
        smallest_buffer = buffer_collection[which(buffer_collection$distance == distance),]   # Extract the smallest buffer from the collection of all buffers 
        smallest_buffer_map = raster::crop(largest_buffer_map, extent(smallest_buffer))   # Crop the largest buffer map to the extent of the smallest buffer (smallest buffer map)
        cliptest = raster::mask(smallest_buffer_map, smallest_buffer)
        fr = fasterize(smallest_buffer, smallest_buffer_map)   # This makes a raster of the sf geometry type   
        area_raster = raster::mask(x=smallest_buffer_map, mask=fr)
        small_geometry = smallest_buffer
      } 
      
      # For all buffers of the collection of buffers except the smallest buffer
      else {
        current_buffer = buffer_collection[which(buffer_collection$distance == distance),]   # Extract the current buffer from the collection of all buffers 
        
        # Subtract the area of the previous buffer (which is smaller than the current buffer) from the current buffer area (a thorus shape is left)
        thorus_geometry = st_difference(current_buffer, small_geometry)
        thorus_map = raster::crop(largest_buffer_map, extent(thorus_geometry))   # Crop the largest buffer map to the extent of the thorus shape
        cliptest = raster::mask(thorus_map, thorus_geometry) 
        fr = fasterize(thorus_geometry, thorus_map)   # This makes a raster of the sf geometry type   
        area_raster = raster::mask(x=thorus_map, mask=fr)
        small_geometry = current_buffer
      }
      
      # Calculate the fourteen different land cover fractions
      # This is done by building a frequency table containing the number of occurences for each of the fourteen different land cover types in the geometry and normalizing this table
      raster_vector_format = table(as.vector(area_raster))   # Create the frequency table   
      
      frequency_table = prop.table(raster_vector_format)   # Normalize the frequency table
      
      # Extract the fourteen different land cover fractions based on the frequency table 
      building = as.numeric(frequency_table["1"])
      building[is.na(building)] = 0
      road = as.numeric(frequency_table["2"])
      road[is.na(road)] = 0
      rest_impervious = as.numeric(frequency_table["3"])
      rest_impervious[is.na(rest_impervious)] = 0
      rail_road = as.numeric(frequency_table["4"])
      rail_road[is.na(rail_road)] = 0
      water = as.numeric(frequency_table["5"])
      water[is.na(water)] = 0
      rest_non_impervious = as.numeric(frequency_table["6"])
      rest_non_impervious[is.na(rest_non_impervious)] = 0
      crop_land = as.numeric(frequency_table["7"])
      crop_land[is.na(crop_land)] = 0
      grass_shrub = as.numeric(frequency_table["8"])
      grass_shrub[is.na(grass_shrub)] = 0
      
      tree = as.numeric(frequency_table["9"])
      tree[is.na(tree)] = 0
      grass_shrub_agriculture = as.numeric(frequency_table["10"])
      grass_shrub_agriculture[is.na(grass_shrub_agriculture)] = 0
      grass_shrub_road = as.numeric(frequency_table["11"])
      grass_shrub_road[is.na(grass_shrub_road)] = 0
      trees_road = as.numeric(frequency_table["12"])
      trees_road[is.na(trees_road)] = 0
      grass_shrub_water = as.numeric(frequency_table["13"])
      grass_shrub_water[is.na(grass_shrub_water)] = 0
      trees_water = as.numeric(frequency_table["14"])
      trees_water[is.na(trees_water)] = 0
      
      watertest = sum(water, grass_shrub_water)
      water = watertest - grass_shrub_water
      
      if (distance == min(distance_vector)){
        # Create a vector with the fourteen land cover fractions that are calculated (for the smallest buffer)
        landcover_buffer = cbind(building, road, rest_impervious, rail_road, water, rest_non_impervious, crop_land, grass_shrub, tree, grass_shrub_agriculture, grass_shrub_road, trees_road, grass_shrub_water, trees_water)
      }
      
      else{
        # Create a vector with the fourteen land cover fractions calculated (for all buffers except the smallest buffer)
        thorus_landcover_buffer = cbind(building, road, rest_impervious, rail_road, water, rest_non_impervious, crop_land, grass_shrub, tree, grass_shrub_agriculture, grass_shrub_road, trees_road, grass_shrub_water, trees_water)   # The land cover fractions for the thorus shape
        circle_landcover_buffer = (thorus_landcover_buffer * distance^2) + (previous_buffer_landcover * previous_distance^2)   # Calculate the land cover fractions for the whole buffer (multiply by the square of the distance to scale the land cover fractions by area)
        landcover_buffer = prop.table(circle_landcover_buffer)   # Normalize the land cover fractions of the whole buffer  
      }
      
      landcover_buffer = as.data.frame(landcover_buffer)
      previous_buffer_landcover = landcover_buffer
      previous_distance = distance
      
      # Store the land cover fractions for all buffers of a certain data point in a dataframe 
      # Each row represents the land cover fractions for a certain buffer distance
      landcover_datapoint = rbind(landcover_datapoint, landcover_buffer)
    }
    
    return(landcover_datapoint)
  }
  
  buffer_dataframe = arrange(buffer_dataframe, buffer_dataframe$Name)
  landcover_dataframe = cbind(as.data.frame(buffer_dataframe), as.data.frame(landcoverlist)) 
  
  # The format of the land cover dataframe is such that the land cover fractions of a certain data point are displayed in multiple rows (one row for each buffer distance)
  
  # Create a vector with the names of the fourteen land cover types
  landcover_classes = c('building', 'road', 'rest_impervious', 'rail_road', 'water', 'rest_non_impervious', 'crop_land', 'grass_shrub', 'tree', 'grass_shrub_agriculture', 'grass_shrub_road', 'trees_road', 'grass_shrub_water', 'trees_water')
  
  # Change the format of the land cover dataframe such that the land cover fractions of a certain data point are displayed in one row with multiple columns indicating the land cover fractions for different buffer distances
  not_landcover_data = landcover_dataframe[which(landcover_dataframe$distance == distance_vector[1]), !(names(landcover_dataframe) %in% landcover_classes)]
  final_landcover_dataframe = data.frame(not_landcover_data)
  for (distance in distance_vector) {
    landcover_data = landcover_dataframe[which(landcover_dataframe$distance == distance), names(landcover_dataframe) %in% landcover_classes]
    distance_string = as.character(distance)
    landcover_colnames = paste(landcover_classes, distance_string, sep="")
    colnames(landcover_data) = landcover_colnames
    final_landcover_dataframe = cbind(final_landcover_dataframe, landcover_data)
  }
  
  return(final_landcover_dataframe)
}

# Function that calculates the green and impervious land cover fractions
calculate_combined_landcover = function(landcover_dataframe, distance_vector){
  
  # Define which land cover types belong to the 'impervious' class and the 'green' class
  impervious_variables = c('building', 'road', 'rest_impervious', 'rail_road')
  green_variables = c('tree', 'rest_non_impervious', 'grass_shrub', 'crop_land', 'grass_shrub_agriculture', 'grass_shrub_road', 'grass_shrub_water', 'trees_water', 'trees_road')
  
  # Add extra columns to the land cover dataframe with the impervious and green land cover fractions for each buffer distance
  for (distance in distance_vector){
    distance_string = as.character(distance)
    impervious_colname = paste0("impervious", distance_string)
    green_colname = paste0("green", distance_string)
    
    impervious_columns = paste(impervious_variables, distance_string, sep="")
    green_columns = paste(green_variables, distance_string, sep="")
    
    # Calculate the impervious and green land cover fractions as the sum of the land cover fractions for the specific land cover types belonging to the impervious and green classes
    landcover_dataframe[impervious_colname] = rowSums(landcover_dataframe[, impervious_columns])
    landcover_dataframe[green_colname] = rowSums(landcover_dataframe[, green_columns])
  }
  
  return(landcover_dataframe)
}

# Function that calculates the land cover fractions for the set of weather stations
# The land cover fractions are calculated based on the BBK land cover map of Flanders
landcover_weather_stations_BBK = function(distance_vector, station_locations, BBK_directory, cores){
  
  # Transform the coordinates from WGS84 system to Lambert-72 system and create buffers for the weather stations
  station_locations = coordinate_projection_lambert(station_locations)
  stationbuffers = creating_buffers_BBK(station_locations, distance_vector)
  
  # Build a small land cover map for the weather stations (based on the BBK land cover map)
  print("Building the BBK small land cover map of the environment of the set of weather stations")
  
  station_landcover_map = creating_small_landcover_map_BBK(BBK_directory, stationbuffers, distance_vector)
  
  print("The small land cover map is built")
  
  # Calculate the land cover fractions for the weather stations (based on the BBK land cover map)
  station_landcover_dataframe = calculate_landcover_BBK(stationbuffers, distance_vector, station_landcover_map, cores)
  
  return(station_landcover_dataframe)
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------TEMPERATURE CORRECTION FUNCTIONS-------------------------------------------------------------------------#

# Function that corrects the measured temperature for thermal inertia
thermal_inertia_correction = function(landcover_dataframe, inertia_files_directory, degrees_of_freedom){
  
  print(paste0("Degrees of freedom used for smoothing of the temperature derivative: ", degrees_of_freedom))
  
  # Create a vector of datafiles used for the thermal inertia correction
  # Each file gives a value of tau which is a measure for the speed of change of the temperature of the sensor when there is a considerable temperature difference
  inertia_files = c("temp_inertia_1", "temp_inertia_2", "temp_inertia_3", "temp_inertia_4", "temp_inertia_5", "temp_inertia_6", "temp_inertia_7")
  
  # Function that calculates the derivative and smoothed derivative of temperature
  temperature_derivative = function(dataframe, degrees_of_freedom){
    
    # Calculate the derivative of temperature and the smoothed derivative of temperature 
    temp_derivative_data = predict(sm.spline(dataframe$fulldate, dataframe$temperature), dataframe$fulldate, 1)
    temp_derivative = temp_derivative_data[,1]
    derivative_dataframe = data.frame("date" = dataframe$fulldate, "temp_derivative" = temp_derivative)
    smoothed_temp_derivative_data = smooth.spline(derivative_dataframe$date, derivative_dataframe$temp_derivative, df = degrees_of_freedom)
    derivative_dataframe$smooth = smoothed_temp_derivative_data$y
    return(derivative_dataframe) 
  }
  
  # Calculate the derivative and smoothed derivative of temperature for the measured data points
  temp_derivative_dataframe = temperature_derivative(landcover_dataframe, degrees_of_freedom)
  smoothed_derivative = temp_derivative_dataframe$smooth
 
  # Plot the derivative of temperature and the smoothed derivative of temperature vs. time
  png("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/temperature_derivative.png")
  plot = ggplot() + geom_line(data = temp_derivative_dataframe, aes(x = date, y = temp_derivative, linetype='Derivative of temperature')) 
  plot2 = plot + geom_line(data = temp_derivative_dataframe, aes(x = date, y = smooth, linetype='Smoothed derivative of temperature')) + labs(x='Time (in UTC)', y="dT/dt (°C/s)") + scale_linetype_manual(values = c("Derivative of temperature" = "dashed", "Smoothed derivative of temperature" = "solid")) + theme(plot.title = element_text(size = 15, hjust = 0), axis.title=element_text(size=17), axis.text=element_text(size=17), legend.title=element_blank(), legend.text=element_text(size=13), legend.position = "top")
  print(plot2)
  dev.off()
  
  
  # Function that calculates tau via a fit of a non-linear model T ~ a + b*(1 - exp(-t*c)) to the data in the inertia datafiles
  calculation_tau = function(inertia_datafile, inertia_files_directory){
    
    # Read the data from the inertia datafile 
    path_inertia_data = paste0(inertia_files_directory, inertia_datafile, ".txt")
    inertia_data = fread(path_inertia_data, sep=",", dec=".")
    names(inertia_data) = c("date", "time", "XWGS", "YWGS", "temperature", "pressure", "humidity")
    inertia_data$fulldate = paste(inertia_data$date, inertia_data$time)
    inertia_data$fulldate = as.POSIXct(inertia_data$fulldate, format="%d/%m/%Y %H:%M:%S", tz=Sys.timezone())
    
    # Create a dataframe with the measured temperature and seconds passed since begin of the signal
    inertia_data$seconds = as.numeric(difftime(inertia_data$fulldate, inertia_data$fulldate[1], units = "secs"))
    temp = inertia_data$temperature
    secs = inertia_data$seconds
    temp_evolution_dataframe = data.frame(secs, temp)
    
    # Use a non-linear model T ~ a + b*(1 - exp(-t*c)) to fit the parameters a, b and c
    fit = nlsfit(data = temp_evolution_dataframe, model = 12, start = c(a = max(inertia_data$temperature), b = min(inertia_data$temperature), c = 0.01))
    
    # Extract the best fit coefficients a, b and c
    coefficient_list = as.numeric(unlist(fit[2])) 
    a = coefficient_list[1]
    b = coefficient_list[2]
    c = coefficient_list[3]
  
    # Calculate tau which is a measure for the speed of change of the temperature of the sensor when there is a considerable temperature difference
    tau = 1/c
    
    return(tau)
  }
  
  # Create a vector used for different values of tau (for each inertia datafile)
  tauvector = c()
  
  # Calculate tau for each inertia datafile and put the values in a vector 
  for (inertia_file in inertia_files){
    tau = calculation_tau(inertia_file, inertia_files_directory)
    tauvector = c(tauvector, tau)
  }
  print(c("Values of tau: ", paste0(tauvector, " s")))
  
  # Create one tau by taking the mean of all the elements in the vector of tau values
  finaltau = mean(tauvector)
  print(paste0("Final tau value: ", finaltau, " s"))
  
  # Calculate the temperature corrected for thermal inertia and store it as a column in the landcover dataframe
  landcover_dataframe$thermal_inertia_corrected = landcover_dataframe$temperature + (finaltau * smoothed_derivative)
  
  return(landcover_dataframe$thermal_inertia_corrected)
}

# Function that calculates the temperature residuals and smoothed temperature residuals for a number of weather stations 
station_smoothed_residuals = function(landcover_dataframe, station_directory_list, station_name_list, number_reference){
  
  # Calculate the start time and end time of the measurements
  starttime = min(landcover_dataframe$fulldate)
  endtime = max(landcover_dataframe$fulldate)

  # Define columnnames for the data of the weather stations
  column_names_station = c("fulldate", "temperature")
  
  # Create a dataframe for the residuals and smoothed residuals of the weather stations. The time of the measurements is already stored in the dataframe
  station_residual_dataframe = data.frame("date" = landcover_dataframe$fulldate)
  
  # Calculate the temperature residuals and smoothed temperature residuals for all specified weather stations
  index=1
  while (index <= length(station_name_list)){
    
    # Read the data from the datafile of the weather station
    stationdata = read.table(station_directory_list[index], header = FALSE, sep = ',', skip = 1)
    names(stationdata) = column_names_station
    stationdata$fulldate = as.POSIXct(stationdata$fulldate, format="%Y-%m-%d %H:%M:%S", tz=Sys.timezone(), origin)
    stationdata$fulldate = as.POSIXct(stationdata$fulldate, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    
    # Only retain the data points that are relevant (same time interval as the measured data points)
    station_subset = subset(stationdata, stationdata$fulldate < (endtime+900) & stationdata$fulldate > (starttime-900))
    # Interpolation of the temperatures of the station data points to the exact same times as the measured data points of the route
    station_interpolated = approx(x = station_subset$fulldate, y = station_subset$temperature, xout = unique(landcover_dataframe$fulldate))
    interpolated_temp = as.numeric(unlist(station_interpolated[2]))

    # Add the station interpolated temperatures to the dataframe together with the temperature residuals and smoothed temperature residuals
    station_residual_dataframe = add_column(station_residual_dataframe, !!(station_name_list[index]) := interpolated_temp )
    station_residual_dataframe = add_column(station_residual_dataframe, !!(paste0('residu_', station_name_list[index])) := station_residual_dataframe[number_reference,2+(index-1)*3]-station_residual_dataframe[,2+(index-1)*3])

    station_smoothed_residuals = smooth.spline(station_residual_dataframe$date, station_residual_dataframe[,index*3], df=5)
    smoothed_residuals = station_smoothed_residuals$y
    station_residual_dataframe = add_column(station_residual_dataframe, !!(paste0('residu_', station_name_list[index], '_smoothed')) := smoothed_residuals)
    
    index = index+1
  }
  
  return(station_residual_dataframe)
}

# Function that calculates the temperature residuals for the temperature decline
temperature_residu_cal = function(landcover_dataframe, stations_smoothed_residuals, stations_landcover, distance_vector, name_list){
  
  # Create a vector with the names of the land cover variables that are used to calculate the temperature residuals
  variable_list = c()
  for(distance in distance_vector){
    distance_string = as.character(distance)
    variables_distance = c(paste0('green', distance_string), paste0('impervious', distance_string), paste0('water', distance_string))
    variable_list = append(variable_list, variables_distance)
  }
  
  # Function that calculates the likelihood between the land cover of a certain location and a certain weather station
  likelihood = function(landcover, stations_landcover) {
    landcover = as.numeric(landcover)
    stations_landcover = as.numeric(stations_landcover)
    weight = 1/sqrt(sum((landcover-stations_landcover)^2))
    weight = replace_na(weight, 9999)
    weight[is.infinite(weight)] = 9999
    weight = sum(weight)
    return(weight)
  }
  
  # Extract the impervious, green and water land cover fractions for the measured data points 
  landcover_measurements = select(landcover_dataframe, variable_list)

  # Extract the impervious, green and water land cover fractions for the weather stations
  landcover_of_stations = select(stations_landcover, c('Name', variable_list))
  residual_dataframe = c()
  weights_dataframe = c()
  for (station in name_list){
    # Extract the smoothed temperature residuals for each station
    smoothed_residual = select(stations_smoothed_residuals, paste0("residu_", station, "_smoothed"))
    smoothed_residual = as.vector(unlist(smoothed_residual))
    residual_dataframe = cbind(residual_dataframe, smoothed_residual)
    
    # Calculate the weight factors for each station
    specific_landcover = subset(landcover_of_stations, Name == station) %>% select(-Name)
    weights = apply(landcover_measurements, 1, likelihood, stations_landcover = as.numeric(specific_landcover[1,]))
    weights_dataframe = cbind(weights_dataframe, weights)
  }
  residual_dataframe = as.data.frame(residual_dataframe)
  weights_dataframe = as.data.frame(weights_dataframe)
  
  # Normalize the weights
  for (row in 1:length(weights_dataframe[,1])){
    weights_sum = sum(weights_dataframe[row,])
    for (column in 1:length(weights_dataframe)){
      weights_dataframe[row, column] = weights_dataframe[row, column]/weights_sum
    }
  }
  weights_dataframe = cbind(weights_dataframe, landcover_dataframe$fulldate)
  
  # Replace infinite weight values by 99999
  weights_dataframe[weights_dataframe == Inf] = 99999 
  
  # Create a dataframe with the smoothed temperature residuals of the stations together with the likelihood coefficients
  residual_weight_dataframe = data.frame(residual_dataframe, weights_dataframe)
  
  # Function that calculates the temperature residuals as a weighted sum of the smoothed temperature residuals of the weather stations
  # The weights are the likelihood coefficients between the land cover of the measured data points and the weather stations
  calculate_temp_residual = function(dataframe){
    
    residuals = as.numeric(dataframe)[1:length(name_list)]   # The smoothed temperature residuals of the stations
    weights = as.numeric(dataframe)[(length(name_list)+1):(2*length(name_list))]   # The likelihood coefficients 
    temp_residu = weighted.mean(residuals, weights)
    return(temp_residu)
  }
  
  # Calculate the temperature residuals 
  temperature_residu = apply(residual_weight_dataframe, 1, calculate_temp_residual)
  
  return(temperature_residu)
}

# Function that corrects the temperatures for temperature decline based on a specific reference data point (reference time) and for a particular set of weather stations
temp_decline_corr = function(landcover_dataframe, landcover_weather_stations, station_directory_list, station_name_list, distance_vector, reference_num, smoothing){
  
  # Calculate the temperature residuals and smoothed temperature residuals for the set of weather stations
  smoothed_residual_dataframe = station_smoothed_residuals(landcover_dataframe, station_directory_list, station_name_list, reference_num)
  
  # Calculate the temperature residuals based on the smoothed temperature residuals of the specified weather stations and the similarity in land cover fractions between the measured data points and the weather stations
  temperature_residu = temperature_residu_cal(landcover_dataframe, smoothed_residual_dataframe, landcover_weather_stations, distance_vector, station_name_list)
  
  if (smoothing){
    # Calculate the smoothed temperature residuals of the measured data points
    smoothed_residuals_list = smooth.spline(landcover_dataframe$fulldate, temperature_residu, df=5)
    smoothed_residuals = smoothed_residuals_list$y
    
    # Correct the temperatures for temperature decline by adding the final smoothed temperature residuals
    corrected_temp = landcover_dataframe$temperature + smoothed_residuals
  }
  
  else{
    # Correct the temperatures for temperature decline by adding the final temperature residuals (without smoothing)
    corrected_temp = landcover_dataframe$temperature + temperature_residu
  }
    
  
  return(corrected_temp)
}

# Function that calculates the temperatures corrected for temperature decline at different times in the time range of the measurements
# The begin time of the measurements is always included and the time interval is specified in the general settings in the main script by the 'seconds_per_interval' variable
# These will also be the times for which temperature charts will be made. The decline correction is based on the set of specified weather stations
# The temperature decline correction is based on the BBK map of Flanders
temp_decline_timeseries = function(landcover_dataframe, seconds_per_interval, landcover_weather_stations, station_directory_list, station_name_list, distance_vector){
  
  # Create a vector of the numbers of the reference data points that are used for different temperature decline corrections
  # The time interval between two consecutive reference data points is specified by the 'seconds_per_interval' variable
  vector_ref_numbers = seq(1, length(landcover_dataframe$fulldate), by = seconds_per_interval/10)
  
  # Perform the temperature decline for each reference time 
  for (reference_number in vector_ref_numbers){
    print(paste0("Correcting for temperature decline at reference time: ", landcover_dataframe$fulldate[reference_number], " based on the BBK land cover map of Flanders"))
    corrected_temperature = temp_decline_corr(landcover_dataframe, landcover_weather_stations, station_directory_list, station_name_list, distance_vector, reference_number, FALSE)
    
    # Add the decline corrected temperatures to a dataframe
    if (reference_number == 1){
      dataframe_corr_temp_timeseries = data.frame(corrected_temperature)
      names(dataframe_corr_temp_timeseries) = paste0("corrected_", reference_number)
    }
    else{
      dataframe_corr_temp_timeseries = add_column(dataframe_corr_temp_timeseries, !!(paste0("corrected_", reference_number)) := corrected_temperature)
    }
  }
  
  return(dataframe_corr_temp_timeseries)
}

# Function that plots the evolution of: the temperatures corrected for temperature decline based on the set of weather stations and smoothing of the final temperature residuals,
# the temperatures corrected for temperature decline based on the set of weather stations and no smoothing of the final temperature residuals, 
# the temperatures corrected for thermal inertia but not for temperature decline and the measured uncorrected temperatures
plot_temp_evolutions = function(temperature_dataframe){
  
  png("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/temp_evolutions.png")
  plot = ggplot(data = temperature_dataframe) + geom_line(aes(x = date, y = decline_corrected_no_smoothing, color = "Decline corrected no smoothing"), size = 0.4) + geom_line(aes(x = date, y = inertia_corrected, color = "Inertia corrected (no decline)"), size = 0.4) + geom_line(aes(x = date, y = decline_corrected, color = "Decline corrected"), size = 0.4) + geom_line(aes(x = date, y = measured, color = "Measured"), size = 0.4)
  plot2 = plot + scale_color_manual(values = c("Measured" = "orange2", "Inertia corrected (no decline)" = "blue","Decline corrected no smoothing" = "green", "Decline corrected" = "black")) 
  plot3 = plot2 + labs(x = 'Time (in UTC)', y = 'Temperature (°C)') + theme(plot.title = element_text(hjust = 0.5), legend.position="top", legend.title=element_blank(), legend.text=element_text(size=13), axis.title=element_text(size=15), axis.text=element_text(size=15)) + guides(color=guide_legend(ncol=2))
  print(plot3)
  dev.off()
  
  return(temperature_dataframe)
}

# Function that creates a plot of the temperature residuals and smoothed temperature residuals of the specified weather stations during the period of measurements
residu_plot_stations = function(station_smoothed_residuals, name_list_station, simple_residu){
  
  # Create a dataframe with the temperature residuals and smoothed temperature residuals
  station_smoothed_residuals = station_smoothed_residuals[,!(names(station_smoothed_residuals) %in% name_list_station)]
  station_smoothed_residuals$residu_traditional = simple_residu
  station_smoothed_residuals$residu_traditional_smoothed = simple_residu
  residual_data = station_smoothed_residuals[, c(1, seq(2, length(names(station_smoothed_residuals)), by = 2))]
  smoothed_residual_data = station_smoothed_residuals[, c(1, seq(3, length(names(station_smoothed_residuals)), by = 2))]
  
  melted_residual = melt(residual_data, id.vars = "date")
  melted_smoothed_residual = melt(smoothed_residual_data, id.vars = "date")
  
  # For the Antwerp weather stations
  names_vector = c("AWV-Boechout (1)", "AWV-Kontich (2)", "AWV-Ranst (3)", "AWV-Stabroek (4)", "Infrabel-Berchem (5)", "Infrabel-Ekeren (6)", "Infrabel-Lier (7)", "VLINDER-Beveren (8)", "VLINDER-Eilandje (9)", "VLINDER-ITG (10)", "VLINDER-Kontich (11)", "VLINDER-Linkeroever (12)", "VLINDER-Zoo", "VMM-Lier", "VMM-Luchtbal (13)", "VMM-Melsele (14)", "VMM-Stabroek", "VMM-Vremde", "WOW-Borgerhout (15)", "WOW-Kruibeke", "traditional")
  color_vec = c("black", "gray50", "darkblue", "blue", "cyan", "aquamarine2", "darkgreen", "green", "darkolivegreen1", "burlywood1", "red", "darkkhaki", "darkorange3", "saddlebrown", "orange", "yellow", "plum1", "lightseagreen", "purple", "darkmagenta", "olivedrab")
  linetype_vec = rep(c("solid", "twodash"), times = length(names_vector))
  
  melted_smoothed_residual_used = melted_smoothed_residual[which(!melted_smoothed_residual$variable %in% c("residu_VLINDER-Zoo_smoothed", "residu_VMM-Lier_smoothed", "residu_VMM-Stabroek_smoothed", "residu_VMM-Vremde_smoothed", "residu_WOW-Kruibeke_smoothed")),]
  melted_smoothed_residual_not_used = melted_smoothed_residual[which(melted_smoothed_residual$variable %in% c("residu_VLINDER-Zoo_smoothed", "residu_VMM-Lier_smoothed", "residu_VMM-Stabroek_smoothed", "residu_VMM-Vremde_smoothed", "residu_WOW-Kruibeke_smoothed")),]
  
  # For the Ghent weather stations
  #names_vector = c("AWV-Deinze (1)", "AWV-Destelbergen (2)", "AWV-Evergem (3)", "AWV-Nevele", "AWV-Oosterzele (4)", "AWV-Terdonk (5)", "Infrabel-St-Pieters (6)", "MOCCA-Botanic garden (7)", "MOCCA-Harbour (8)", "MOCCA-Melle (9)", "MOCCA-St-Bavo (10)", "MOCCA-Wondelgem (11)", "VLINDER-Gentbrugge (12)", "VLINDER-Melle (13)", "VLINDER-Observatory (14)", "VLINDER-Oostakker (15)", "VLINDER-Ottogracht (16)", "VLINDER-UGent (17)", "WOW-Belfort (18)", "WOW-De Pinte (19)", "WOW-Destelbergen (20)", "WOW-Lovendegem (21)", "WOW-Merendree (22)", "WOW-Nazareth (23)", "WOW-Wetteren (24)", "WOW-Zaffelare (25)", "WOW-Zevergem (26)", "WOW-Zwijnaarde (27)", "traditional")
  #color_vec = c("black", "gray50", "darkblue", "blue", "cyan", "aquamarine2", "darkgreen", "green", "darkolivegreen1", "burlywood1", "red", "darkkhaki", "darkorange3", "saddlebrown", "orange", "yellow", "plum1", "lightseagreen", "purple", "darkmagenta", "gray80", "hotpink", "lightblue", "gold2", "lightpink", "seagreen1", "khaki", "magenta", "olivedrab")
  #linetype_vec = rep(c("solid", "twodash"), times = length(names_vector))
  
  #melted_smoothed_residual_used = melted_smoothed_residual[which(!melted_smoothed_residual$variable %in% c("residu_AWV-Nevele_smoothed")),]
  #melted_smoothed_residual_not_used = melted_smoothed_residual[which(melted_smoothed_residual$variable %in% c("residu_AWV-Nevele_smoothed")),]
  
  # Plot the temperature residuals and smoothed temperature residuals of the weather stations
  pdf(paste0("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/temp_res_weather_station.pdf"), width = 5.5, height=7)
  plot = ggplot() + geom_line(data = melted_residual, aes(x = date, y = value, color = variable, linetype = "dashed"), size = 1) + geom_line(data = melted_smoothed_residual_used, aes(x = date, y = value, color = str_sub(variable, 1, -10), linetype = "solid"), size = 1) + geom_line(data = melted_smoothed_residual_not_used, aes(x = date, y = value, color = str_sub(variable, 1, -10), linetype = "solid"), size = 2)                                                                               
  plot2 = plot + scale_color_manual(labels = names_vector, values = color_vec) + labs(x = 'Time (in UTC)', y = 'Temperature residual (°C)') + guides(color = guide_legend(ncol = 2), linetype=FALSE) + theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), legend.position = "top", legend.text=element_text(size=11), axis.title=element_text(size=15), axis.text=element_text(size=15))
  print(plot2)
  dev.off()
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------CHART FUNCTIONS-------------------------------------------------------------------------------------#

# Function that creates a dataframe with the coordinates (in Lambert-72 system) of the grid points of the temperature chart together with a number
making_chart_dataframe_BBK = function(dataframe, step_distance, number_horizontal_steps, number_vertical_steps, chart_extension){
  
  # Initialize a dataframe
  chart_dataframe = data.frame(c(), c(), c())
  
  # Initialize the starting point used for creating the grid
  point = c(min(dataframe$X) - chart_extension - step_distance, max(dataframe$Y) + chart_extension + step_distance, 0) 
  
  # Create a grid of equally spaced points. The distance between the gridpoints is specified by the 'step_distance' variable in the general settings
  # The coordinates of the grid points are stored in a dataframe together with a number. The grid covers the area of the route with an additional spacing that is given by the 'chart_extension' variable starting from the four extremity points of the route
  # The grid is created from the upper left corner to the bottom right corner
  for (step in 1:number_vertical_steps){
    point[2] = point[2] - step_distance
    for (step in 1:number_horizontal_steps){
      point[1] = point[1] + step_distance
      point[3] = point[3] + 1
      chart_dataframe = rbind(chart_dataframe, point)
    }
    point[1] = min(dataframe$X) - chart_extension - step_distance
  }
  
  names(chart_dataframe) = c("X_lambert", "Y_lambert", "Name")

  return(chart_dataframe)
}

# Function that calculates the numbers of the grid points where the weather stations used for validation of the chart temperatures are located
location_numbers_on_chart = function(chart_dataframe, locations_list, step_distance){
  
  # Initialize a vector for the numbers of the grid points where the weather stations are located
  vector_loc_numbers = c()
  
  # Calculate for each location the number of the grid point that is closest to the location
  for (location in locations_list){
    distance = 2*step_distance
    
    # Loop over all grid points
    for (row in 1:length(chart_dataframe$X)){
      # Calculate the distance between the location and the grid point
      new_distance = abs(location[1]-chart_dataframe$X[row]) + abs(location[2]-chart_dataframe$Y[row])
      if (new_distance < distance){
        distance = new_distance
        row_number_location = row
      }
      
    }
    vector_loc_numbers = c(vector_loc_numbers, row_number_location)
  }
  
  return(vector_loc_numbers)
}

# Function that calculates the indices of the grid points where the weather stations used for validation of the chart temperatures are located
matrix_indices_locations = function(vector_location_numbers, number_horizontal_steps, number_vertical_steps){
  
  # Initialize vectors for the X and Y indices of the grid points of the locations of the weather stations
  X_indices = c()
  Y_indices = c()
  
  # Transform the numbers of the grid points of the locations to pairs of indices
  for (loc_number in vector_location_numbers){
    row = floor(loc_number/number_horizontal_steps)
    if (mod(loc_number, number_horizontal_steps) != 0){
      row = row + 1
      column = mod(loc_number, number_horizontal_steps)
    }
    if (mod(loc_number, number_horizontal_steps) == 0){
      column = number_horizontal_steps
    }
    X_indices = c(X_indices, column)
    Y_indices = c(Y_indices, row)
  }
  Y_indices = abs(Y_indices - number_vertical_steps) + 1
  
  # Put the X and Y indices of the locations in a list
  indices_list = list("X_indices" = X_indices, "Y_indices" = Y_indices)
  return(indices_list)
}

# Function that removes weather stations used for validation of the chart temperatures from which the location doesn't fit in the temperature chart
remove_locations = function(dataframe_chart, location_list, stations_validation_names, step_distance){
  
  # Create a vector with the four extremity points of the temperature chart
  chart_extremities = c(min(dataframe_chart$X), max(dataframe_chart$X), min(dataframe_chart$Y), max(dataframe_chart$Y))
  
  # Initialize a vector for the indices of locations that should be removed
  removed_indices = c()
  
  # Check for each location if it fits in the temperature chart
  for (index in 1:length(location_list)){
    location = unlist(location_list[index])
    if (location[1] < (chart_extremities[1] - step_distance/2) | location[1] > (chart_extremities[2] + step_distance/2) | location[2] < (chart_extremities[3] - step_distance/2) | location[2] > (chart_extremities[4] + step_distance/2)){
      print(paste0("Location: ", stations_validation_names[index], " doesn't fit in the temperature chart"))
      removed_indices = c(removed_indices, index)
    }
  }
  
  # Remove locations that doesn't fit in the temperature chart
  if (length(removed_indices) > 0){
    location_list = location_list[-removed_indices]
    stations_validation_names = stations_validation_names[-removed_indices]
  }
  
  returnlist = list("locations" = location_list, "names" = stations_validation_names, "removed_indices" = removed_indices)
  return(returnlist)
}

# Function that creates a dataframe with the model variables for the grid points of the temperature chart
# These model variables are used to predict the temperatures of the grid points of the chart
chart_variables_selection =  function(chart_dataframe, model_variables){
  
  # Remove rows in the dataframe with NA values
  chart_dataframe = drop_na(chart_dataframe)
  
  # Select model variables and create a new dataframe with these variables 
  chart_variables_dataframe = chart_dataframe %>% dplyr::select(model_variables)
  chart_variables_dataframe = chart_variables_dataframe[!is.infinite(rowSums(chart_variables_dataframe)),]   # Remove rows with infinite values
  
  return(chart_variables_dataframe)
}

# Function that creates a dataframe with the model variables and temperature
data_variables =  function(dataframe, model_variables){
  
  # Remove rows in the dataframe with NA values
  dataframe = drop_na(dataframe)
  
  # Select model variables and create a new dataframe with these variables and temperature
  variables_dataframe = dataframe %>% dplyr::select(model_variables)
  variables_dataframe$temperature = dataframe$temperature
  variables_dataframe = variables_dataframe[!is.infinite(rowSums(variables_dataframe)),]   # Remove rows with infinite values
  
  return(variables_dataframe)
}

# Function that creates a dataframe with the model variables for the grid points of the temperature chart
# The land cover fractions of the grid points of the temperature chart are calculated based on the BBK map 
# The indices of the grid points where the weather stations used for validation of the chart temperatures are located are calculated as well
chart_variables_and_locations = function(landcover_dataframe, num_horizontal_steps, num_vertical_steps, step_distance, extension_chart, buffer_distances_chart, loc_list, stations_validation_names, model_variables){
  
  # Create a dataframe with the coordinates (in Lambert-72 system) of the grid points of the temperature chart together with a number
  chart_dataframe = making_chart_dataframe_BBK(landcover_dataframe, step_distance, num_horizontal_steps, num_vertical_steps, extension_chart)   
    
  # Create buffers for the grid points of the temperature chart
  chart_bufferdata = creating_buffers_BBK(chart_dataframe, buffer_distances_chart)
    
  print("Building the BBK small land cover map for the temperature charts")
    
  # Create a small land cover map for the temperature chart (based on the BBK land cover map of Flanders)
  chart_small_landcover_map = creating_small_landcover_map_BBK(BBK_directory, chart_bufferdata, buffer_distances_chart)
    
  print("The small land cover map is built")
    
  print("Calculating the land cover fractions for the grid points of the temperature charts (based on the BBK land cover map of Flanders)")
    
  # Calculate land cover fractions for the grid points of the temperature chart (based on the BBK land cover map of Flanders)
  chart_landcover = calculate_landcover_BBK(chart_bufferdata, buffer_distances_chart, chart_small_landcover_map, cores)
  chart_total_landcover = calculate_combined_landcover(chart_landcover, buffer_distances_chart)
    
  print("The land cover fractions are calculated")
  
  print("Removing locations that doesn't fit in the temperature charts")
  
  # Remove weather stations used for validation of the chart temperatures from which the location doesn't fit in the temperature chart
  loc_list = remove_locations(chart_total_landcover, loc_list, stations_validation_names, step_distance)
  # Calculate the numbers of the grid points where the weather stations used for validation of the chart temperatures are located
  vector_location_numbers = location_numbers_on_chart(chart_total_landcover, loc_list$locations, step_distance)
  
  # Calculate the indices of the grid points where the weather stations used for validation of the chart temperatures are located
  list_indices_loc = matrix_indices_locations(vector_location_numbers, num_horizontal_steps, num_vertical_steps)
  
  # Create a dataframe with the model variables for the grid points of the temperature chart
  # These model variables are used to predict the temperatures of the grid points of the chart
  chart_variables = chart_variables_selection(chart_total_landcover, model_variables)
  
  chart_variables_loc_indices = list("chart_variables" = chart_variables, "indices_of_locations" = list_indices_loc, "validation_names" = loc_list$names, "removed_indices" = loc_list$removed_indices)
  
  return(chart_variables_loc_indices)
}

# Function that calculates the positions of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements
# Temperature can not accurately be predicted for these grid points
chart_selection = function(route_variables, chart_variables, max_distance){
  
  # Initialize a vector for the positions of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements
  positions = c()
  
  # Compare the land cover fractions of each grid point to the land cover fractions of the measured data points
  for (chart_row in 1:length(chart_variables[,1])){
    for (route_row in 1:length(route_variables[,1])){
      
      # Calculate a measure for the similarity between the land cover of a grid point and a measured data point
      distancevector = route_variables[route_row,] - chart_variables[chart_row,]
      distance = sqrt(sum((distancevector)^2))
      
      # If the similarity measure exceeds a particular specified limit, go to the next grid point
      if (distance < max_distance){
        break
      }
      
      # If the similarity measure doesn't exceed the particular specified limit for all data points, the position of that grid point is stored
      if (route_row == length(route_variables[,1])){
        positions = c(positions, chart_row)
      }
    }
  }
  
  return(positions)
}

# Function that predicts the temperatures of the grid points of the temperature chart based on a linear model with as variables the variables in 'model_variables'
# The linear model is built with temperatures corrected for temperature decline at a specific reference time specified by the 'reference_number' variable
# The positions of grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements are also calculated
pred_chart_temp_and_selection = function(landcover_dataframe, chart_landcover_data, model_variables, reference_number, corr_temp_timeseries, step_distance){
  
  # Extract the temperatures corrected for temperature decline at the specific reference time specified by the 'reference_number' variable
  corr_temperatures = corr_temp_timeseries[[paste0("corrected_", reference_number)]]
  landcover_dataframe$temperature = corr_temperatures
  
  # Create a dataframe with the model variables and temperature of the measurements
  data_variables_temp = data_variables(landcover_dataframe, model_variables)
  
  # Build a linear model based on the model variables and temperature of the measurements
  chart_linear_model = lm(temperature ~ ., data = data_variables_temp)
  
  # Predict the temperatures of the grid points of the temperature chart based on the linear model
  chart_temperatures = predict(chart_linear_model, newdata = chart_landcover_data)
  
  # Create a dataframe with only the model variables of the measurements
  data_variables_landcover = data_variables_temp[1:(length(data_variables_temp) - 1)]  
  if (reference_number == 1){
    # Calculate the positions of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements
    positions_no_temp = chart_selection(data_variables_landcover, chart_landcover_data, max_distance = 0.11)
  }
  else{
    positions_no_temp = 0
  }
  
  info = list("chart_temperatures" = chart_temperatures, "positions_no_temp" = positions_no_temp, "linear_model" = chart_linear_model)
  return(info)
}

# Function that creates a list of temperatures that will be used in the temperature charts 
chart_timeseries = function(landcover_dataframe, num_horizontal_steps, num_vertical_steps, model_variables, step_distance, corr_temp_timeseries, chart_variables){
  
  # Create a vector of the numbers of the reference data points that are used for different temperature decline corrections
  # The different temperature charts are made based on these different temperature decline corrected temperatures
  # The time interval between two consecutive reference data points is specified by the 'seconds_per_interval' variable
  vector_ref_numbers = seq(1, length(landcover_dataframe$fulldate), by = seconds_per_interval/10)
  
  # Initialize a list for the temperatures of the different temperature charts
  chart_list = vector(mode = "list", length = length(vector_ref_numbers))
  index_list = 1
  
  # Calculate the temperatures of the grid points of the temperature chart for each reference time
  # Calculate the positions of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements
  for (ref_number in vector_ref_numbers){
    
    # Calculate the temperatures of the grid points of the temperature chart and the positions of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements
    temperatures_positions = pred_chart_temp_and_selection(landcover_dataframe, chart_variables, model_variables, ref_number, corr_temp_timeseries, step_distance)
    
    # Add the temperatures of the grid points to the list for the temperatures of the different temperature charts
    chart_list[[index_list]] = temperatures_positions$chart_temperatures
    if (ref_number == 1){
      
      # Extract the positions of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements
      positions_no_temp = temperatures_positions$positions_no_temp
    }
    index_list = index_list + 1
  }
  
  returnlist = list("chart_temperatures_list" = chart_list, "positions_no_temp" = positions_no_temp)
  return(returnlist)
}

# Function that transforms a vector with certain values to a matrix of these values
making_chart = function(vector, number_horizontal_steps, number_vertical_steps){
  
  # Initialize an empty matrix with the specified number of rows and columns
  chart_matrix = matrix(nrow=number_vertical_steps, ncol=number_horizontal_steps)
  
  # Fill the matrix with the values that are in the vector 
  # The matrix is filled from the upper left corner to the bottom right corner
  for (rownumber in 1:number_vertical_steps){
    for (columnnumber in 1:number_horizontal_steps){
      number = (rownumber-1)*number_horizontal_steps + columnnumber
      value = vector[number]
      chart_matrix[rownumber, columnnumber] = value
    }
  }
  return(chart_matrix)
}

# Function that plots the temperature charts 
plot_chart = function(chart, dataframe, number_horizontal_steps, number_vertical_steps, step_distance, number_reference, min_temp, max_temp, indices_of_locations, stations_validation_names, correction_type){
  
  # This is needed for plotting the scalebar at the right place
  if (number_horizontal_steps > 0 & number_horizontal_steps < 11){
    index=2
  }
  if (number_horizontal_steps > 10 & number_horizontal_steps < 15){
    index=3
  }
  if (number_horizontal_steps > 14 & number_horizontal_steps < 19){
    index=4
  }
  if (number_horizontal_steps > 18 & number_horizontal_steps < 23){
    index=5
  }
  if (number_horizontal_steps > 22 & number_horizontal_steps < 30){
    index=6
  }
  if (number_horizontal_steps > 29 & number_horizontal_steps < 37){
    index=8
  }
  if (number_horizontal_steps > 36 & number_horizontal_steps < 65){
    index=10
  }
  if (number_horizontal_steps > 64 & number_horizontal_steps < 89){
    index=20
  }
  if (number_horizontal_steps > 88 & number_horizontal_steps < 109){
    index=30
  }
  if (number_horizontal_steps > 108 & number_horizontal_steps < 128){
    index=40
  }
  if (number_horizontal_steps > 127 & number_horizontal_steps < 145){
    index=50
  }
  if (number_horizontal_steps > 144){
    index=60
  }

  
  # Create a color palette used for coloring the grid points of the temperature chart
  colors_chart = c("ivory", "lemonchiffon1","#FFEDA0","#FED976","#FEB24C","#FD8D3C","orangered","firebrick")

  # Get the indices and names of the stations in Antwerp that should be plotted at the top, bottom, right and left of the marking spot on the temperature charts
  X_indices_top = indices_of_locations$X_indices[c(1, 7, 12)]
  X_indices_bottom = indices_of_locations$X_indices[15]
  X_indices_right = indices_of_locations$X_indices[c(3, 4, 5, 6, 8, 9, 13, 14)]
  X_indices_left = indices_of_locations$X_indices[c(2, 10, 11)]
  Y_indices_top = indices_of_locations$Y_indices[c(1, 7, 12)]
  Y_indices_bottom = indices_of_locations$Y_indices[15]
  Y_indices_right = indices_of_locations$Y_indices[c(3, 4, 5, 6, 8, 9, 13, 14)]
  Y_indices_left = indices_of_locations$Y_indices[c(2, 10, 11)]
  stations_top = stations_validation_names[c(1, 7, 12)]
  stations_bottom = stations_validation_names[15]
  stations_right = stations_validation_names[c(3, 4, 5, 6, 8, 9, 13, 14)]
  stations_left = stations_validation_names[c(2, 10, 11)]

  # Get the indices and names of the stations in Ghent that should be plotted at the top, bottom, right and left of the marking spot on the temperature charts
  #X_indices_top = indices_of_locations$X_indices[c(4, 5, 15, 16, 19, 21)]
  #X_indices_bottom = indices_of_locations$X_indices[c(9, 11, 22, 26)]
  #X_indices_right = indices_of_locations$X_indices[c(1, 2, 7, 8, 10, 12, 17, 20, 23, 27)]
  #X_indices_left = indices_of_locations$X_indices[c(3, 6, 13, 14, 18, 24, 25)]
  #Y_indices_top = indices_of_locations$Y_indices[c(4, 5, 15, 16, 19, 21)]
  #Y_indices_bottom = indices_of_locations$Y_indices[c(9, 11, 22, 26)]
  #Y_indices_right = indices_of_locations$Y_indices[c(1, 2, 7, 8, 10, 12, 17, 20, 23, 27)]
  #Y_indices_left = indices_of_locations$Y_indices[c(3, 6, 13, 14, 18, 24, 25)]
  #stations_top = stations_validation_names[c(4, 5, 15, 16, 19, 21)]
  #stations_bottom = stations_validation_names[c(9, 11, 22, 26)]
  #stations_right = stations_validation_names[c(1, 2, 7, 8, 10, 12, 17, 20, 23, 27)]
  #stations_left = stations_validation_names[c(3, 6, 13, 14, 18, 24, 25)]
  
  # Plot the temperature chart
  png(paste0("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/temp_chart_LINEAR_", number_reference, ".png"), width = 720, height = 720)   
  chartplot = plot(chart, axis.col = NULL, axis.row = NULL, na.col = "gray60", main = "", xlab = "", ylab = "", key = list(side = 2, cex.axis = 1.2, font = 2), spacing.key = c(2,1.5,0), col = colors_chart, breaks = seq(from=min_temp-0.1, to=max_temp+0.1, length.out=9), border = NA) 
  # Plot the locations of the weather stations used for validation of the temperatures of the temperature chart
  points(x = indices_of_locations$X_indices, y = indices_of_locations$Y_indices, col = "black", pch = 5, cex = 100/max(number_horizontal_steps, number_vertical_steps), lwd = 4-(max(number_horizontal_steps, number_vertical_steps)*0.01))
  # Add a north arrow to the plot 
  addnortharrow(pos = "topright", cols = c("black", "white"), border = "black", text.col = "black", lwd = 3)
  # Add a scalebar to the plot
  plotscalebar(x = round(number_horizontal_steps/15), y = round(number_vertical_steps/15), ht = ceil(number_vertical_steps/40), params = scalebarparams(plotunit = "m", widthhint = 0.2 + 0.1*number_horizontal_steps/70), style = "bar")
  textxy(X = c(round(number_horizontal_steps/15), round(number_horizontal_steps/15) + index), Y = c(round(number_vertical_steps/15) + number_vertical_steps/16, round(number_vertical_steps/15) + number_vertical_steps/16), labs = c("0", paste0(index*step_distance/1000, "")), font = 2, cex = 2, offset=0)
  # Add the names of the weather stations used for validation of the temperatures of the temperature chart
  textxy(X = X_indices_top, Y = Y_indices_top, labs = stations_top, cex = 2-(max(number_horizontal_steps, number_vertical_steps)*0.005), offset = 0.55, col = "black", font = 2, pos = 3)
  textxy(X = X_indices_bottom, Y = Y_indices_bottom, labs = stations_bottom, cex = 2-(max(number_horizontal_steps, number_vertical_steps)*0.005), offset = 0.55, col = "black", font = 2, pos = 1)
  textxy(X = X_indices_right, Y = Y_indices_right, labs = stations_right, cex = 2-(max(number_horizontal_steps, number_vertical_steps)*0.005), offset = 0.55, col = "black", font = 2, pos = 4)
  textxy(X = X_indices_left, Y = Y_indices_left, labs = stations_left, cex = 2-(max(number_horizontal_steps, number_vertical_steps)*0.005), offset = 0.55, col = "black", font = 2, pos = 2)
  textxy(X = round(ncol(chart)/15) + index + round(ncol(chart)/20) - 1, Y = round(nrow(chart)/15) + 2, labs = c("km"), font = 2, cex = 2, offset=0)
  print(chartplot)
  dev.off()
  
}

# Function that creates the names of the validation weather stations that are plotted on the temperature charts
# The name consists of the observed temperature for the specific station followed by a number
create_station_names = function(landcover_dataframe, seconds_per_interval, removed_indices, names_vector, val_files){
  
  # Remove stations that doesn't fit in the temperature charts
  if (length(removed_indices) != 0){
    val_files = val_files[-removed_indices]
  }
 
  # Extract the reference times and put them in a dataframe
  vector_ref_numbers = seq(1, length(landcover_dataframe$fulldate), by = seconds_per_interval/10)
  ref_times = landcover_dataframe$fulldate[vector_ref_numbers]
  obs_data = data.frame(ref_times)
  
  # Create a dataframe with the observations of the validation weather stations at the reference times
  for (index in 1:length(val_files)){
    # Read the data of the particular validation station
    dataframe_measurements = read.table(val_files[index], header = FALSE, sep = ',', skip = 1)
    dataframe_measurements$fulldate = as.POSIXct(dataframe_measurements[,1], format="%Y-%m-%d %H:%M:%S", tz=Sys.timezone(), origin)
    dataframe_measurements$fulldate = as.POSIXct(dataframe_measurements$fulldate, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    
    # Initialize a vector with the temperatures of the particular validation station at the reference times
    ref_temperatures = c()
    
    # Calculate the temperatures of the particular validation station at the reference times
    for (ref_number in 1:length(ref_times)){
      ref_time = ref_times[ref_number]
      for (number in 1:length(dataframe_measurements$fulldate)){
        time = dataframe_measurements$fulldate[number]
        if (number == 1){
          min_time = as.numeric(difftime(time, ref_time, units = "secs"))
          number_min = 1
        }
        else{
          time_diff = as.numeric(difftime(time, ref_time, units = "secs"))
          if (abs(time_diff) < abs(min_time)){
            number_min = number
            min_time = time_diff
          }
        }
      }
      ref_temperatures = c(ref_temperatures, dataframe_measurements[number_min, 2])
    }
    
    # Add the vector with the temperatures at the reference times to a dataframe
    obs_data = cbind(obs_data, ref_temperatures)
    
  }
  
  # Initialize a list for the names of the validation stations
  station_names_list = vector(mode = "list", length = length(vector_ref_numbers))
  # Fill the list
  for (i in 1:length(vector_ref_numbers)){
    station_names_list[[i]] = paste0(as.vector(format(round(obs_data[i, 2:(1+length(names_vector))], digits = 1), nsmall = 1)), "°C [", names_vector, "]")
  }
  
  return(station_names_list)
}

# Function that creates the temperature charts of the environment of the bicycle route  
# The measurements of the route are used as input for the linear model with which the charts are made
# The charts are made at specific reference times. The begin time of the measurements is always included and the time interval is specified by the 'seconds_per_interval' variable in the general settings of the main script
create_temp_chart = function(chart_list, num_horizontal_steps, num_vertical_steps, seconds_per_interval, landcover_dataframe, step_distance, minimum_temp, maximum_temp, list_indices_loc, list_names, correction_type){
  
  # Initialize a list for the temperatures of the different temperature charts
  chart_list_temp = vector(mode = "list", length = length(chart_list$chart_temperatures_list))
  
  # Fill the list for the temperatures of the different temperature charts
  for (index in 1:length(chart_list$chart_temperatures_list)){
    # Extract the temperatures of the grid points of the temperature chart
    chart_temperatures = chart_list$chart_temperatures_list[[index]]
    # Set temperatures of the grid points with land cover fractions that are not likely to the land cover fractions that are encountered in the measurements to NA value
    chart_temperatures[chart_list$positions_no_temp] = NA
    
    # Add the temperatures of the temperature chart to the list for the temperatures of the different temperature charts
    chart_list_temp[[index]] = chart_temperatures
  }
  
  # Plot each of the temperature charts
  for (i in 1:length(chart_list_temp)){
    # Extract the temperatures of the chart
    chart_temp = chart_list_temp[[i]]
    
    # Transform the vector with the chart temperatures to a matrix
    chart = making_chart(chart_temp, num_horizontal_steps, num_vertical_steps)
    
    # Calculate the reference number of the data point that is used as reference for the temperature decline correction and thus as reference for the temperature chart
    reference_num = 1 + (i-1)*(seconds_per_interval/10)
    
    # Plot the temperature chart
    plot_chart(chart, landcover_dataframe, num_horizontal_steps, num_vertical_steps, step_distance, reference_num, minimum_temp, maximum_temp, list_indices_loc, list_names[[i]], correction_type) 
  }
  
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------VALIDATION FUNCTIONS-----------------------------------------------------------------------------#

# Function that creates a dataframe of the observed and predicted temperatures of the validation stations at the reference times
validation = function(val_files, landcover_dataframe, removed_indices, names_vector, landcover_data, chart_variables, model_variables, corr_temp_timeseries, step_distance, distance_vec){
  
  # Remove stations that doesn't fit in the temperature charts
  if (length(removed_indices) != 0){
    val_files = val_files[-removed_indices]
    landcover_data = landcover_data[-removed_indices,]
  }
  
  # Extract the reference times and put them in a dataframe
  vector_ref_numbers = seq(1, length(landcover_dataframe$fulldate), by = seconds_per_interval/10)
  ref_times = landcover_dataframe$fulldate[vector_ref_numbers]
  obs_data = data.frame(ref_times)
  
  # Create a dataframe with the observations of the validation weather stations at the reference times
  for (index in 1:length(val_files)){
    # Read the data of the particular validation station
    dataframe_measurements = read.table(val_files[index], header = FALSE, sep = ',', skip = 1)
    
    dataframe_measurements$fulldate = as.POSIXct(dataframe_measurements[,1], format="%Y-%m-%d %H:%M:%S", tz=Sys.timezone(), origin)
    dataframe_measurements$fulldate = as.POSIXct(dataframe_measurements$fulldate, format="%Y-%m-%d %H:%M:%S", tz="UTC")
    
    # Initialize a vector for the temperatures of the particular validation station at the reference times
    
    # Calculate the temperatures of the particular validation station at the reference times
    ref_temperatures = c()
    for (ref_number in 1:length(ref_times)){
      ref_time = ref_times[ref_number]
      for (number in 1:length(dataframe_measurements$fulldate)){
        time = dataframe_measurements$fulldate[number]
        if (number == 1){
          min_time = as.numeric(difftime(time, ref_time, units = "secs"))
          number_min = 1
        }
        else{
          time_diff = as.numeric(difftime(time, ref_time, units = "secs"))
          if (abs(time_diff) < abs(min_time)){
            number_min = number
            min_time = time_diff
          }
        }
      }
      ref_temperatures = c(ref_temperatures, dataframe_measurements[number_min, 2])
      
    }
    
    # Add the vector with the temperatures at the reference times to a dataframe
    obs_data = cbind(obs_data, ref_temperatures)
    
  }
  
  # Create a dataframe with the predictions of the validation weather stations at the reference times
  pred_dataframe = data.frame(matrix(ncol = length(landcover_data$Name), nrow=0))
  for (ref_number in vector_ref_numbers){
    # Create a linear model based on the temperatures corrected for temperature decline at the reference time
    temp_data_and_model = pred_chart_temp_and_selection(landcover_dataframe, chart_variables, model_variables, ref_number, corr_temp_timeseries, step_distance)
    lin_model = temp_data_and_model$linear_model
    
    # Select land cover types used to predict the temperature
    col_names = c()
    for (distance in distance_vec){
      col_names = c(col_names, paste(c("water", "impervious", "green"), distance, sep=""))
    }
    final_landcover_data = select(landcover_data, col_names)
    
    # Calculate the predicted temperatures and store them in a dataframe
    pred_temp = predict(lin_model, newdata = final_landcover_data)
    pred_dataframe = rbind(pred_dataframe, pred_temp)
    
  }
  
  pred_data = cbind(obs_data[,1], pred_dataframe)
  names(obs_data) = c("fulldate", paste0(names_vector, "_obs"))
  names(pred_data) = c("fulldate", paste0(names_vector, "_pred"))
  
  melted_obs = melt(obs_data, id.vars = "fulldate", variable.name = "station")
  melted_pred = melt(pred_data, id.vars = "fulldate", variable.name = "station")
  
  return(list("observed" = melted_obs, "predicted" = melted_pred))
}

# Function that calculates the RMSE, BIAS and MAE scores based on the observed and predicted temperatures of the validation weather stations
rmse_mae_bias = function(validation_data){
  
  # Extract the observations and predictions of temperature of the validation stations
  observations = unlist(validation_data$observed[3])
  predictions = unlist(validation_data$predicted[3])
  
  # Calculate the RMSE, BIAS and MAE scores
  rmse = sqrt(mean((observations - predictions)^2))
  mae = mean(abs(observations - predictions))
  bias = mean(predictions - observations)
  
  return(list("rmse" = rmse, "bias" = bias, "mae" = mae))
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------ROUTE VISUALISATION FUNCTIONS------------------------------------------------------------------------------#

# Function that visualises the route on the BBK land cover map of Flanders
visualisation_route_BBK = function(landcover_dataframe, data_landcover_map){
  
  names(data_landcover_map) = c("long", "lat", "BBK_Flanders2015")
  # Visualise the route
  pdf("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/visualisation_route_BBK.pdf", width = 10, height = 8) 
  plot = ggplot() + geom_raster(data = data_landcover_map, aes(x = long, y = lat, fill = cut(round(BBK_Flanders2015), breaks = c(0, 1, 4, 5, 14)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  plot2 = plot + scale_fill_manual(values = c("#b36b6b", "#bfbfbf", "#99dfff", "#cfff8c"), guide=FALSE)                                                                                                    
  plot3 = plot2 + geom_point(data = landcover_dataframe, aes(x = X, y = Y), col = "black", shape = 19, size = 1) + guides(alpha = FALSE) + geom_point(data = landcover_dataframe, aes(x = X[1], y = Y[1]), col = "cyan", shape = 5, size = 3.5, stroke = 1.5) + geom_point(data = landcover_dataframe, aes(x = X[length(landcover_dataframe$X)], y = Y[length(landcover_dataframe$Y)]), col = "black", shape = 5, size = 3.5, stroke = 1.5) + theme(plot.title = element_text(size = 13), legend.title=element_text(size=13), legend.text=element_text(size=13))
  plot4 = plot3 + north(data_landcover_map, symbol=12) + scalebar(data_landcover_map, location = "bottomright", dist = 2, dist_unit = "km", transform = FALSE, st.size=7)
  print(plot4)
  dev.off()
}

# Function that visualises the measured uncorrected temperatures of the route 
# The route is plotted on the BBK land cover map of Flanders and covers a detailed section of the complete route
visualisation_route_measured_temp_BBK_detailed = function(landcover_dataframe, data_landcover_map, inertia_bool, decline_bool, corr_temp_timeseries = 0){
  
  if (inertia_bool == FALSE & decline_bool == FALSE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$measured_temp)
    max_temp = max(landcover_dataframe$measured_temp)
  }
  
  if (inertia_bool == TRUE & decline_bool == FALSE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
    max_temp = max(landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
  }
  
  if (inertia_bool == FALSE & decline_bool == TRUE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$measured_temp, corr_temp_timeseries)
    max_temp = max(landcover_dataframe$measured_temp, corr_temp_timeseries)
  }
  
  if (inertia_bool == TRUE & decline_bool == TRUE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(corr_temp_timeseries, landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
    max_temp = max(corr_temp_timeseries, landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
  }  
  # Calculate temperature breaks for legend
  breaks_temp = seq(from=min_temp, to=max_temp, length.out=11)[-c(1,11)]
    
  names(data_landcover_map) = c("long", "lat", "BBK_Flanders2015")
  # Visualise the measured temperatures 
  pdf("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/visualisation_measured_temp_BBK_detailed.pdf", width = 10, height = 8) 
  plot = ggplot() + geom_raster(data = data_landcover_map, aes(x = long, y = lat, fill = cut(round(BBK_Flanders2015), breaks = c(0, 1, 4, 5, 14)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  plot2 = plot + scale_fill_manual(values = c("#b36b6b", "#bfbfbf", "#99dfff", "#cfff8c"), guide=FALSE)                                                                                                    
  plot3 = plot2 + geom_point(data = landcover_dataframe, aes(x = X, y = Y, col = measured_temp), shape = 19, size = 3) + scale_color_stepsn(colors = c("blue4", "steelblue", "skyblue3", "paleturquoise", "beige", "khaki1", "sandybrown", "sienna2", "red", "firebrick4"), breaks = breaks_temp, labels = paste0(round(breaks_temp, digits=1), "°C"), limits = c(min_temp, max_temp), name = "Temperature") + guides(alpha = FALSE) + theme(plot.title = element_text(size = 13), legend.title=element_text(size=19), legend.text=element_text(size=19), legend.position = c(0.1, 0.73), legend.key.height=unit(3,"line"), legend.key.width=unit(1, 'cm'))
  plot4 = plot3 + north(data_landcover_map, symbol=12) + scalebar(data_landcover_map, location = "bottomright", dist = 2, dist_unit = "km", transform = FALSE, st.size=7)
  print(plot4)
  dev.off()
  
}

# Function that visualises the measured uncorrected temperatures of the route 
# The route is plotted on the BBK land cover map of Flanders
visualisation_route_measured_temp_BBK = function(landcover_dataframe, data_landcover_map, inertia_bool, decline_bool, corr_temp_timeseries = 0){
  
  if (inertia_bool == FALSE & decline_bool == FALSE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$measured_temp)
    max_temp = max(landcover_dataframe$measured_temp)
  }
  
  if (inertia_bool == TRUE & decline_bool == FALSE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
    max_temp = max(landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
  }
  
  if (inertia_bool == FALSE & decline_bool == TRUE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$measured_temp, corr_temp_timeseries)
    max_temp = max(landcover_dataframe$measured_temp, corr_temp_timeseries)
  }
  
  if (inertia_bool == TRUE & decline_bool == TRUE){
    # Calculate the minimum and maximum measured temperatures
    min_temp = min(corr_temp_timeseries, landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
    max_temp = max(corr_temp_timeseries, landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
  }
  # Calculate temperature breaks for legend
  breaks_temp = seq(from=min_temp, to=max_temp, length.out=11)[-c(1,11)]

  names(data_landcover_map) = c("long", "lat", "BBK_Flanders2015")

  # Visualise the measured temperatures 
  pdf("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/visualisation_measured_temp_BBK.pdf", width = 10, height = 8) 
  plot = ggplot() + geom_raster(data = data_landcover_map, aes(x = long, y = lat, fill = cut(round(BBK_Flanders2015), breaks = c(0, 1, 4, 5, 14)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  plot2 = plot + scale_fill_manual(values = c("#b36b6b", "#bfbfbf", "#99dfff", "#cfff8c"), guide=FALSE)                                                                                                    
  plot3 = plot2 + geom_point(data = landcover_dataframe, aes(x = X, y = Y, col = measured_temp), shape = 19, size = 2) + scale_color_stepsn(colors = c("blue4", "steelblue", "skyblue3", "paleturquoise", "beige", "khaki1", "sandybrown", "sienna2", "red", "firebrick4"), breaks = breaks_temp, labels = paste0(round(breaks_temp, digits=1), "°C"), limits = c(min_temp, max_temp), name = "Temperature") + guides(alpha = FALSE) + geom_point(data = landcover_dataframe, aes(x = X[1], y = Y[1]), col = "cyan", shape = 5, size = 5, stroke = 1.5) + geom_point(data = landcover_dataframe, aes(x = X[length(landcover_dataframe$X)], y = Y[length(landcover_dataframe$Y)]), col = "black", shape = 5, size = 5, stroke = 1.5) + theme(plot.title = element_text(size = 13), legend.title=element_text(size=15), legend.text=element_text(size=15), legend.position = c(0.16, 0.79), legend.key.height=unit(3,"line"), legend.key.width=unit(1, 'cm')) 
  plot4 = plot3 + north(data_landcover_map, symbol=12) + scalebar(data_landcover_map, location = "bottomright", dist = 2, dist_unit = "km", transform = FALSE, st.size=7)
  print(plot4)
  dev.off()
  
}

# Function that visualises the temperatures of the route that are corrected for thermal inertia
# The route is plotted on the BBK land cover map of Flanders and covers a detailed section of the complete route
visualisation_route_thermal_inertia_BBK_detailed = function(landcover_dataframe, data_landcover_map, decline_bool, corr_temp_timeseries = 0){
  
  if (decline_bool == FALSE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
    max_temp = max(landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
  }
  
  if (decline_bool == TRUE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(corr_temp_timeseries, landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
    max_temp = max(corr_temp_timeseries, landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
  }
  # Calculate temperature breaks for legend
  breaks_temp = seq(from=min_temp, to=max_temp, length.out=11)[-c(1,11)]
  
  names(data_landcover_map) = c("long", "lat", "BBK_Flanders2015")
  # Visualise the temperatures corrected for thermal inertia
  pdf("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/visualisation_thermal_inertia_BBK_detailed.pdf", width = 10, height = 8) 
  plot = ggplot() + geom_raster(data = data_landcover_map, aes(x = long, y = lat, fill = cut(round(BBK_Flanders2015), breaks = c(0, 1, 4, 5, 14)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  plot2 = plot + scale_fill_manual(values = c("#b36b6b", "#bfbfbf", "#99dfff", "#cfff8c"), guide=FALSE)                                                                                                    
  plot3 = plot2 + geom_point(data = landcover_dataframe, aes(x = X, y = Y, col = thermal_inertia), shape = 19, size = 3) + scale_color_stepsn(colors = c("blue4", "steelblue", "skyblue3", "paleturquoise", "beige", "khaki1", "sandybrown", "sienna2", "red", "firebrick4"), breaks = breaks_temp, labels = paste0(round(breaks_temp, digits=1), "°C"), limits = c(min_temp, max_temp), name = "Temperature") + guides(alpha = FALSE) + theme(plot.title = element_text(size = 13), legend.title=element_text(size=19), legend.text=element_text(size=19), legend.position = c(0.1, 0.73), legend.key.height=unit(3,"line"), legend.key.width=unit(1, 'cm'))
  plot4 = plot3 + north(data_landcover_map, symbol=12) + scalebar(data_landcover_map, location = "bottomright", dist = 2, dist_unit = "km", transform = FALSE, st.size=7)
  print(plot4)
  dev.off()
}

# Function that visualises the temperatures of the route that are corrected for thermal inertia
# The route is plotted on the BBK land cover map of Flanders
visualisation_route_thermal_inertia_BBK = function(landcover_dataframe, data_landcover_map, decline_bool, corr_temp_timeseries = 0){
  
  if (decline_bool == FALSE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
    max_temp = max(landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
  }
  
  if (decline_bool == TRUE){
    # Calculate the minimum and maximum temperatures
    min_temp = min(corr_temp_timeseries, landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
    max_temp = max(corr_temp_timeseries, landcover_dataframe$thermal_inertia, landcover_dataframe$measured_temp)
  }
  # Calculate temperature breaks for legend
  breaks_temp = seq(from=min_temp, to=max_temp, length.out=11)[-c(1,11)]
  
  names(data_landcover_map) = c("long", "lat", "BBK_Flanders2015")
  # Visualise the temperatures corrected for thermal inertia
  pdf("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/visualisation_thermal_inertia_BBK.pdf", width = 10, height = 8) 
  plot = ggplot() + geom_raster(data = data_landcover_map, aes(x = long, y = lat, fill = cut(round(BBK_Flanders2015), breaks = c(0, 1, 4, 5, 14)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  plot2 = plot + scale_fill_manual(values = c("#b36b6b", "#bfbfbf", "#99dfff", "#cfff8c"), guide=FALSE)                                                                                                    
  plot3 = plot2 + geom_point(data = landcover_dataframe, aes(x = X, y = Y, col = thermal_inertia), shape = 19, size = 3) + scale_color_stepsn(colors = c("blue4", "steelblue", "skyblue3", "paleturquoise", "beige", "khaki1", "sandybrown", "sienna2", "red", "firebrick4"), breaks = breaks_temp, labels = paste0(round(breaks_temp, digits=1), "°C"), limits = c(min_temp, max_temp), name = "Temperature") + geom_point(data = landcover_dataframe, aes(x = X[1], y = Y[1]), col = "cyan", shape = 5, size = 5, stroke = 1.5) + geom_point(data = landcover_dataframe, aes(x = X[length(landcover_dataframe$X)], y = Y[length(landcover_dataframe$Y)]), col = "black", shape = 5, size = 5, stroke = 1.5) + guides(alpha = FALSE) + theme(plot.title = element_text(size = 13), legend.title=element_text(size=19), legend.text=element_text(size=19), legend.position = c(0.1, 0.73), legend.key.height=unit(3,"line"), legend.key.width=unit(1, 'cm'))
  plot4 = plot3 + north(data_landcover_map, symbol=12) + scalebar(data_landcover_map, location = "bottomright", dist = 2, dist_unit = "km", transform = FALSE, st.size=7)
  print(plot4)
  dev.off()
}

# Function that visualises the temperatures of the route for a specific reference time
# The temperatures are corrected for temperature decline (based on the BBK land cover map of Flanders) and the route is plotted on the BBK land cover map of Europe
visualisation_route_temp_decline_BBK = function(landcover_dataframe, data_landcover_map, ref_number, min_temp, max_temp, correction_type){
  
  # Calculate temperature breaks for legend
  breaks_temp = seq(from=min_temp, to=max_temp, length.out=11)[-c(1,11)]
  
  names(data_landcover_map) = c("long", "lat", "BBK_Flanders2015")
  # Visualise the temperatures corrected for temperature decline
  pdf(paste0("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/visualisation_temp_decline_BBK_",correction_type, "_", ref_number, ".pdf"), width = 10, height = 8)
  plot = ggplot() + geom_raster(data = data_landcover_map, aes(x = long, y = lat, fill = cut(round(BBK_Flanders2015), breaks = c(0, 1, 4, 5, 14)))) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), plot.title = element_text(hjust = 0), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
  plot2 = plot + scale_fill_manual(values = c("#b36b6b", "#bfbfbf", "#99dfff", "#cfff8c"), guide=FALSE)                                                                                                    
  plot3 = plot2 + geom_point(data = landcover_dataframe, aes(x = X, y = Y, col = temperature), shape = 19, size = 3) + scale_color_stepsn(colors = c("blue4", "steelblue", "skyblue3", "paleturquoise", "beige", "khaki1", "sandybrown", "sienna2", "red", "firebrick4"), breaks = breaks_temp, labels = paste0(round(breaks_temp, digits=1), "°C"), limits = c(min_temp, max_temp), name = "Temperature") + geom_point(data = landcover_dataframe, aes(x = X[1], y = Y[1]), col = "cyan", shape = 5, size = 5, stroke = 1.5) + geom_point(data = landcover_dataframe, aes(x = X[length(landcover_dataframe$X)], y = Y[length(landcover_dataframe$Y)]), col = "black", shape = 5, size = 5, stroke = 1.5) + guides(alpha = FALSE) + theme(plot.title = element_text(size = 13), legend.title=element_text(size=19), legend.text=element_text(size=19), legend.position = c(0.1, 0.73), legend.key.height=unit(3,"line"), legend.key.width=unit(1, 'cm'))
  plot4 = plot3 + north(data_landcover_map, symbol=12) + scalebar(data_landcover_map, location = "bottomright", dist = 2, dist_unit = "km", transform = FALSE, st.size=7)
  print(plot4)
  dev.off()
}

# Function that visualises the temperatures of the route corrected for temperature decline (based on the BBK land cover map of Flanders)
# Different temperature decline corrections are used based on different reference data points
# The result is a series of plots of the temperatures of the route at different reference times and the route is plotted on the BBK land cover map of Flanders
# The begin time of the measurements is always included and the time interval is specified by the 'seconds_per_interval' variable in the general settings of the main script
visualisation_route_timeseries_BBK = function(landcover_dataframe, dataframe_landcover_map, seconds_per_interval, corr_temp_timeseries, correction_type){
  
  # Create a vector of the numbers of the reference data points that are used for different temperature decline corrections
  # The time interval between two consecutive reference data points is specified by the 'seconds_per_interval' variable
  vector_ref_numbers = seq(1, length(landcover_dataframe$fulldate), by = seconds_per_interval/10)
  
  # Calculate the minimum and maximum temperatures 
  minimum_temp = min(corr_temp_timeseries, landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
  maximum_temp = max(corr_temp_timeseries, landcover_dataframe$measured_temp, landcover_dataframe$thermal_inertia)
  
  # Plot the temperatures of the route for each temperature decline correction (reference time)
  for (reference_number in vector_ref_numbers){
    
    # Extract the temperature decline correction corresponding to the specified reference number
    corr_temperatures = corr_temp_timeseries[[paste0("corrected_", reference_number)]]
    landcover_dataframe$temperature = corr_temperatures
    
    # Visualise the temperatures of the route for the specific reference time
    visualisation_route_temp_decline_BBK(landcover_dataframe, dataframe_landcover_map, reference_number, minimum_temp, maximum_temp, correction_type)
  }
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
