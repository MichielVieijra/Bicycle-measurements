# THIS IS A SCRIPT THAT CAN BE USED TO ANALYSE MOBILE BICYCLE MEASUREMENTS AND TO CREATE TEMPERATURE CHARTS BASED ON THESE MEASUREMENTS

# Load packages that are required to run the program
library(numbers)
library(plyr)
library(ggsn)
library(devtools)
library(raster)
library(rgeos)
library(data.table)
library(ggplot2)
#library(ggtern)  # UNCOMMENT THIS FOR CREATING THE TRIANGLE PLOTS
library(tibble)
library(GISTools)
library(maps)
library(base)
library(prettymapr)
library(postGIStools)
library(RColorBrewer)
library(proj4)
library(rgdal)
library(sf)
library(foreach)
library(Metrics)
library(doParallel)
library(tcltk2)
library(spex)
library(fasterize)
library(Hmisc)
library(PerformanceAnalytics)
library(plotly)
library(dplyr)
library(geosphere)
library(stringr)
library(lubridate)
library(pspline)
library(reshape2)
library(rasterVis)
library(doMC)
library(viridis)
library(tmap)
library(leaflet)
library(rpart)
library(rpart.plot)
library(plotrix)
library(fields)
library(tidyr)
library(nlstools)
library(plot.matrix)
library(easynls)
library(calibrate)

# Make a connection to the script with all the necessary functions
source("/home/mivieijra/Documents/VLINDER/VLINDER/Programs/functions_script.R")

#--------------------------------------------------------------------------------GENERAL SETTINGS---------------------------------------------------------------------------------#

print("General settings are loaded")
thermal_inertia = TRUE    # If TRUE the temperature is corrected for thermal inertia
temp_decline = TRUE    # If TRUE the temperature is corrected for temperature decline
chart = TRUE   # If TRUE temperature charts are created

# Specify the number of cores that will be used for calculation of the land cover fractions
cores = 3
registerDoMC(cores)

# Specify the distance(s) for which circular buffers have to be created (in meter)
buffer_distances = c(250)

# Specify the variables that will be used to predict the temperatures of the grid points of the temperature chart
# These variables are given as input to the linear model from which the temperature chart is created
model_variables = c("impervious250", "green250", "water250")

# Create a vector with all distances that are appearing in the land cover input variables for the linear model (model variables)
buffer_distances_chart = c(250)

step_distance = 250   # The distance between two grid points (in meter)

seconds_per_interval = 1800   # The time interval between reference times for which temperature charts are created (in seconds)
if (mod(seconds_per_interval, 10) != 0){
  stop("The number of seconds of the time interval between two temperature charts should be divisible by 10", call. = FALSE)
}

# The distance that is added to the four extremity points of the route (north, east, south and west) to set the scale of the temperature chart
extension_of_chart = 11000

datafile_measurements = "/home/mivieijra/Documents/VLINDER/VLINDER/Bicycle_measurements/Antwerp_measurements.csv"   # Directory of the datafile with the mobile measurements
BBK_directory = "/home/mivieijra/Documents/VLINDER/VLINDER/Landcover_maps/BBK_Flanders2015.tif"   # Directory of the BBK land cover map of Flanders

# Directory of the folder with the datafiles used for the thermal inertia correction
# MAKE SURE THE INERTIA FILES IN THIS FOLDER ARE IN TXT-FORMAT
inertiafiles_directory = "/home/mivieijra/Documents/VLINDER/VLINDER/Inertia_files/"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# PUT HERE THE DIRECTORIES OF THE DATAFILES WITH THE TEMPERATURE MEASUREMENTS OF THE WEATHER STATIONS 
# THAT SHOULD BE TAKEN INTO ACCOUNT FOR THE TEMPERATURE DECLINE CORRECTION
# MAKE SURE THAT THE FILES HAVE THE FOLLOWING FORMAT: FIRST COLUMN THE DATE (IN FORMAT YYYY-MM-DD HH:MM:SS), SECOND COLUMN THE MEASURED TEMPERATURE
# THE TWO COLUMNS SHOULD BE SEPARATED BY A "," SIGN WITH NO SPACINGS

# Stations for temperature decline of Antwerp measurements
vlinder_kontich = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Kontich.csv"
vlinder_beveren = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Beveren.csv"
vlinder_ITG = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_ITG.csv"
vlinder_linkeroever = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Linkeroever.csv"
vlinder_eilandje = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Eilandje.csv"
vlinder_zoo = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Zoo.csv"
infr_berchem = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Infrabel_Berchem.csv"
infr_ekeren = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Infrabel_Ekeren.csv"
infr_lier = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Infrabel_Lier.csv"
AWV_boechout = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Boechout.csv"
AWV_kontich = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Kontich.csv"
AWV_ranst = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Ranst.csv"
AWV_stabroek = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Stabroek.csv"
VMM_lier = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/VMM_Lier.csv"
VMM_melsele = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/VMM_Melsele.csv"
VMM_vremde = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/VMM_Vremde.csv"
VMM_stabroek = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/VMM_Stabroek.csv"
VMM_luchtbal = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/VMM_Luchtbal.csv"
WOW_kruibeke = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Kruibeke.csv"
WOW_borgerhout = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Borgerhout.csv"

# Stations for temperature decline of Ghent measurements
#infr_st_pieters = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Infrabel_St_Pieters.csv"
#AWV_destelbergen = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Destelbergen.csv"
#AWV_evergem = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Evergem.csv"
#AWV_oosterzele = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Oosterzele.csv"
#AWV_terdonk = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Terdonk.csv"
#AWV_nevele = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Nevele.csv"
#AWV_deinze = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/AWV_Deinze.csv"
#WOW_belfort = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Belfort.csv"
#WOW_zaffelare = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Zaffelare.csv"
#WOW_nazareth = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Nazareth.csv"
#WOW_de_pinte = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_De_Pinte.csv"
#WOW_wetteren = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Wetteren.csv"
#WOW_lovendegem = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Lovendegem.csv"
#WOW_merendree = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Merendree.csv"
#WOW_zevergem = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Zevergem.csv"
#WOW_zwijnaarde = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Zwijnaarde.csv"
#WOW_destelbergen = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/WOW_Destelbergen.csv"
#vlinder_melle = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Melle.csv"
#vlinder_oostakker = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Oostakker.csv"
#vlinder_gentbrugge = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Gentbrugge.csv"
#vlinder_ugent = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_UGent.csv"
#vlinder_ottogracht = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Ottogracht.csv"
#vlinder_observatory = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/Vlinder_Observatory.csv"
#mocca_melle = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/MOCCA_Melle.csv"
#mocca_st_bavo = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/MOCCA_St_Bavo.csv"
#mocca_wondelgem = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/MOCCA_Wondelgem.csv"
#mocca_harbour = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/MOCCA_harbour.csv"
#mocca_botan = "/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/MOCCA_botan.csv"

# CREATE A VECTOR WITH THE DIRECTORIES OF THE WEATHER STATIONS DATA

# Stations for temperature decline of Antwerp measurements
station_path_list = c(AWV_boechout, AWV_kontich, AWV_ranst, AWV_stabroek, infr_berchem, infr_ekeren, infr_lier, vlinder_beveren, vlinder_eilandje, vlinder_ITG, vlinder_kontich, vlinder_linkeroever, VMM_luchtbal, VMM_melsele, WOW_borgerhout)

# Stations for temperature decline of Ghent measurements
#station_path_list = c(AWV_deinze, AWV_destelbergen, AWV_evergem, AWV_oosterzele, AWV_terdonk, infr_st_pieters, mocca_botan, mocca_harbour, mocca_melle, mocca_st_bavo, mocca_wondelgem, vlinder_gentbrugge, vlinder_melle, vlinder_observatory, vlinder_oostakker, vlinder_ottogracht, vlinder_ugent, WOW_belfort, WOW_de_pinte, WOW_destelbergen, WOW_lovendegem, WOW_merendree, WOW_nazareth, WOW_wetteren, WOW_zaffelare, WOW_zevergem, WOW_zwijnaarde)

# CREATE A VECTOR WITH THE NAMES OF THE WEATHER STATIONS
# MAKE SURE THAT THE ORDER OF THE NAMES CORRESPONDS TO THE ORDER OF THE DIRECTORIES AND THE NAMES ARE ORDERED ALPHABETICALLY

# Stations for temperature decline of Antwerp measurements
station_name_list = c("AWV-Boechout", "AWV-Kontich", "AWV-Ranst", "AWV-Stabroek", "Infrabel-Berchem", "Infrabel-Ekeren", "Infrabel-Lier", "VLINDER-Beveren", "VLINDER-Eilandje", "VLINDER-ITG", "VLINDER-Kontich", "VLINDER-Linkeroever", "VMM-Luchtbal", "VMM-Melsele", "WOW-Borgerhout")

# Stations for temperature decline of Ghent measurements
#station_name_list = c("AWV-Deinze", "AWV-Destelbergen", "AWV-Evergem", "AWV-Oosterzele", "AWV-Terdonk", "Infrabel-St-Pieters", "MOCCA-Botanic garden", "MOCCA-Harbour", "MOCCA-Melle", "MOCCA-St-Bavo", "MOCCA-Wondelgem", "VLINDER-Gentbrugge", "VLINDER-Melle", "VLINDER-Observatory", "VLINDER-Oostakker", "VLINDER-Ottogracht", "VLINDER-UGent", "WOW-Belfort", "WOW-De Pinte", "WOW-Destelbergen", "WOW-Lovendegem", "WOW-Merendree", "WOW-Nazareth", "WOW-Wetteren", "WOW-Zaffelare", "WOW-Zevergem", "WOW-Zwijnaarde")

# CREATE A DATAFRAME WITH THE LOCATIONS OF THE WEATHER STATIONS THAT SHOULD BE USED FOR THE TEMPERATURE DECLINE CORRECTION
# THE FORMAT OF THE DATAFRAME SHOULD BE: (NAMES OF STATIONS, X WGS COORDINATES, Y WGS COORDINATES)
# MAKE SURE THAT THE NAMES ARE THE SAME AS IN THE 'STATION_NAME_LIST' VECTOR

# Stations for temperature decline of Antwerp measurements
locations_of_stations = data.frame(c("VLINDER-Kontich", "VLINDER-Beveren", "VLINDER-ITG", "VLINDER-Linkeroever", "VLINDER-Eilandje", "VLINDER-Zoo", "Infrabel-Berchem", "Infrabel-Ekeren", "Infrabel-Lier", "AWV-Boechout", "AWV-Kontich", "AWV-Ranst", "AWV-Stabroek", "VMM-Lier", "VMM-Melsele", "VMM-Vremde", "VMM-Stabroek", "VMM-Luchtbal", "WOW-Kruibeke", "WOW-Borgerhout"), c(4.449670, 4.293436, 4.398711, 4.381624, 4.417377, 4.423266, 4.441457, 4.421149, 4.555863, 4.492301, 4.42919, 4.54813, 4.35373, 4.572627, 4.266921, 4.54872, 4.363919, 4.424406, 4.308052, 4.439055), c(51.13478, 51.26685, 51.21209, 51.22249, 51.23884, 51.21730, 51.18855, 51.28671, 51.13534, 51.168832, 51.128294, 51.209277, 51.316133, 51.131951, 51.243935, 51.16601, 51.324885, 51.260992, 51.171938, 51.212704))
colnames(locations_of_stations) = c("Name", "X_WGS_GPS", "Y_WGS_GPS")
# Order the stations alphabetically
locations_of_stations = arrange(locations_of_stations, locations_of_stations$Name)
# Remove invalid stations
locations_of_stations = locations_of_stations[which(!locations_of_stations$Name %in% c("VLINDER-Zoo", "VMM-Lier", "VMM-Stabroek", "VMM-Vremde", "WOW-Kruibeke")),]
# Check if the names of the set of weather stations that are taken into account are valid
names_check(locations_of_stations, station_name_list)

# Stations for temperature decline of Ghent measurements
#locations_of_stations = data.frame(c("AWV-Deinze", "AWV-Destelbergen", "AWV-Evergem", "AWV-Nevele", "AWV-Oosterzele", "AWV-Terdonk", "Infrabel-St-Pieters", "MOCCA-Botanic garden", "MOCCA-Harbour", "MOCCA-Melle", "MOCCA-St-Bavo", "MOCCA-Wondelgem", "VLINDER-Gentbrugge", "VLINDER-Melle", "VLINDER-Observatory", "VLINDER-Oostakker", "VLINDER-Ottogracht", "VLINDER-UGent", "WOW-Belfort", "WOW-De Pinte", "WOW-Destelbergen", "WOW-Lovendegem", "WOW-Merendree", "WOW-Nazareth", "WOW-Wetteren", "WOW-Zaffelare", "WOW-Zevergem", "WOW-Zwijnaarde"), c(3.539556, 3.819851, 3.65903, 3.534848, 3.811545, 3.775575, 3.704171, 3.722470, 3.749000, 3.815744, 3.732000, 3.702875, 3.791667, 3.815763, 3.724490, 3.789441, 3.728000, 3.709695, 3.7248, 3.64, 3.7827, 3.6398, 3.5749, 3.59059, 3.8807, 3.8811, 3.6815, 3.7139), c(50.97128, 51.049651, 51.087875, 51.057791, 50.950897, 51.142188, 51.03698, 51.03570, 51.10900, 50.98043, 51.05200, 51.08400, 51.04350, 50.98044, 51.04573, 51.08864, 51.05808, 51.02238, 51.05314, 50.99, 51.0558, 51.1036, 51.0741, 50.95331, 50.9923, 51.1286, 50.9866, 51.0045))
#colnames(locations_of_stations) = c("Name", "X_WGS_GPS", "Y_WGS_GPS")
# Order the stations alphabetically
#locations_of_stations = arrange(locations_of_stations, locations_of_stations$Name)
# Remove invalid stations
#locations_of_stations = locations_of_stations[which(!locations_of_stations$Name %in% c("AWV-Nevele")),]
# Check if the names of the set of weather stations that are taken into account are valid
#names_check(locations_of_stations, station_name_list)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# CREATE A VECTOR WITH THE DIRECTORIES OF THE DATAFILES OF THE WEATHER STATIONS THAT ARE USED FOR VALIDATION OF THE CHART TEMPERATURES 
# THE ORDER OF THE STATIONS SHOULD BE THE SAME AS IN THE 'STATION_PATH_LIST'
validation_files = station_path_list

# CREATE A DATAFRAME WITH WGS84 COORDINATES OF THE LOCATIONS OF THE WEATHER STATIONS THAT WILL BE USED FOR VALIDATING THE CHART TEMPERATURES 
# THE ORDER OF THE STATIONS SHOULD BE THE SAME AS THE CORRESPONDING VALIDATION FILES
validation_loc = locations_of_stations

# Convert the dataframe of WGS84 coordinates to a list of Lambert-72 coordinates 
validation_loc_lambert = WGS_to_lambert(validation_loc)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------READ OUT SECTION-----------------------------------------------------------------------------------#

print("Read out section")

# Reading out the data
dataframe_measurements = reading_data(datafile_measurements)
# Data used for the detailed plot in 'Rivierenhof'
dataframe_measurements_detailed = dataframe_measurements[1:105,]

print("Data is correctly read out")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------LAND COVER SECTION----------------------------------------------------------------------------------#

print("Land cover section")

# Create circular buffers for the distance(s) specified in the general settings (based on the BBK land cover map of Flanders)
buffers_dataframe_BBK = creating_buffers_BBK(dataframe_measurements, buffer_distances)
buffers_dataframe_BBK_detailed = creating_buffers_BBK(dataframe_measurements_detailed, buffer_distances)

# Build the BBK small land cover map
print("Building the BBK small land cover map of the environment of the route")
  
small_landcover_map_BBK = creating_small_landcover_map_BBK(BBK_directory, buffers_dataframe_BBK, buffer_distances)
small_landcover_map_BBK_detailed = creating_small_landcover_map_BBK_detailed(BBK_directory, buffers_dataframe_BBK_detailed, buffer_distances)
  
print("The small land cover map is built")

# Calculate the land cover fractions for all measured data points (based on the BBK land cover map)
print("Calculating the land cover fractions of the measured data points (based on the BBK land cover map of Flanders)")
  
landcover_dataframe_BBK = calculate_landcover_BBK(buffers_dataframe_BBK, buffer_distances, small_landcover_map_BBK, cores)
  
# Calculate the green and impervious land cover fractions for all measured data points
landcover_dataframe_BBK = calculate_combined_landcover(landcover_dataframe_BBK, buffer_distances)

print("The land cover fractions of the measured data points are calculated")

print("Calculating the land cover fractions of the set of weather stations (based on the BBK land cover map of Flanders)")

# Calculate the land cover fractions of the set of weather stations based on the BBK land cover map of Flanders
station_landcover_BBK = landcover_weather_stations_BBK(buffer_distances, locations_of_stations, BBK_directory, cores)
station_landcover_BBK = calculate_combined_landcover(station_landcover_BBK, buffer_distances)

print("The land cover fractions of the weather stations are calculated")

# Create triangle plot for the Antwerp case
#landcover_table = read.table("/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/landcover_of_stations_Antwerp.csv", header = TRUE, sep = ",", dec = ".")
#landcover_table = landcover_table[which(!landcover_table$Name %in% c("VLINDER-Zoo", "VMM-Lier", "VMM-Stabroek", "VMM-Vremde", "WOW-Kruibeke")),]
#colors = c("black", "gray50", "darkblue", "blue", "cyan", "aquamarine2", "darkgreen", "green", "darkolivegreen1", "burlywood1", "red", "darkkhaki", "orange", "yellow", "purple")
#names_vector = c("AWV-Boechout (1)", "AWV-Kontich (2)", "AWV-Ranst (3)", "AWV-Stabroek (4)", "Infrabel-Berchem (5)", "Infrabel-Ekeren (6)", "Infrabel-Lier (7)", "VLINDER-Beveren (8)", "VLINDER-Eilandje (9)", "VLINDER-ITG (10)", "VLINDER-Kontich (11)", "VLINDER-Linkeroever (12)", "VMM-Luchtbal (13)", "VMM-Melsele (14)", "WOW-Borgerhout (15)")
#landcover_table$Name = names_vector

# Create triangle plot for the Ghent case
#landcover_table = read.table("/home/mivieijra/Documents/VLINDER/VLINDER/Station_data/landcover_of_stations_Ghent.csv", header = TRUE, sep = ",", dec = ".")
#landcover_table = landcover_table[which(!landcover_table$Name %in% c("AWV-Nevele")),]
#colors = c("black", "gray50", "darkblue", "cyan", "aquamarine2", "darkgreen", "green", "darkolivegreen1", "burlywood1", "red", "darkkhaki", "darkorange3", "saddlebrown", "orange", "yellow", "plum1", "lightseagreen", "purple", "darkmagenta", "gray80", "hotpink", "lightblue", "gold2", "lightpink", "seagreen1", "khaki", "magenta")
#names_vector = c("AWV-Deinze (1)", "AWV-Destelbergen (2)", "AWV-Evergem (3)", "AWV-Oosterzele (4)", "AWV-Terdonk (5)", "Infrabel-St-Pieters (6)", "MOCCA-Botanic garden (7)", "MOCCA-Harbour (8)", "MOCCA-Melle (9)", "MOCCA-St-Bavo (10)", "MOCCA-Wondelgem (11)", "VLINDER-Gentbrugge (12)", "VLINDER-Melle (13)", "VLINDER-Observatory (14)", "VLINDER-Oostakker (15)", "VLINDER-Ottogracht (16)", "VLINDER-UGent (17)", "WOW-Belfort (18)", "WOW-De Pinte (19)", "WOW-Destelbergen (20)", "WOW-Lovendegem (21)", "WOW-Merendree (22)", "WOW-Nazareth (23)", "WOW-Wetteren (24)", "WOW-Zaffelare (25)", "WOW-Zevergem (26)", "WOW-Zwijnaarde (27)")
#landcover_table$Name = names_vector

#pdf("/home/mivieijra/Documents/VLINDER/VLINDER/Produced_data/triangle_plot.pdf", width=6, height=12)
#plot = ggtern() + geom_mask() + geom_point(data = landcover_table, aes(x=water250, y=impervious250, z=green250, color=landcover_table$Name), pch=17, cex=5) + scale_color_manual(labels = landcover_table$Name, values=colors) + geom_point(data=landcover_dataframe_BBK, aes(x=water250, y=impervious250, z=green250, fill=landcover_dataframe_BBK$Name), cex=1.2, pch=21) + scale_fill_gradientn(colours = brewer.pal(11, "Spectral"), guide="colourbar", breaks=c(1, 571, 1140), labels=c("18:14 UTC", "19:49 UTC", "21:24 UTC")) 
#plot2 = plot + labs(x = "", xarrow = "Water 250m", y = "",  yarrow = "Impervious 250m", z = "", zarrow = "Green 250m")  + theme_rgbw() + theme_arrowlength(start=c(0.25,0.25,0.25), finish=c(0.75,0.75,0.75)) + theme(tern.axis.title.L = element_text(vjust = 2)) + theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), legend.position="top", legend.box = "vertical", legend.text=element_text(size=11.5), axis.title=element_text(size=15), axis.text=element_text(size=15)) + guides(color = guide_legend(ncol = 2), fill = guide_colourbar(barwidth=10, direction = "horizontal"))
#print(plot2)
#dev.off()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------VISUALISATION OF ROUTE----------------------------------------------------------------------------------#

print("Visualising the route on the BBK land cover map of Flanders")
  
# Create a duplicate column with the measured temperatures
landcover_dataframe_BBK$measured_temp = landcover_dataframe_BBK$temperature

# Lower the resolution of the small land cover map
small_landcover_map_low_resolution_BBK = aggregate_blocks(small_landcover_map_BBK, res_factor=10)
dataframe_landcover_map_BBK = as.data.frame(small_landcover_map_low_resolution_BBK, xy = TRUE)
  
small_landcover_map_low_resolution_BBK_detailed = aggregate_blocks(small_landcover_map_BBK_detailed, res_factor=10)
dataframe_landcover_map_BBK_detailed = as.data.frame(small_landcover_map_low_resolution_BBK_detailed, xy = TRUE)

# Visualise the route on the BBK land cover map
visualisation_route_measured_temp_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK, FALSE, FALSE)
visualisation_route_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK)

print("The route is visualised")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------TEMPERATURE CORRECTION SECTION--------------------------------------------------------------------------------#

# Apply a correction for thermal inertia 
if (thermal_inertia){
  
  print("Thermal inertia section")
  
  # Specify the degrees of freedom used for smoothing of the temperature derivative
  degrees_of_freedom = round(NROW(landcover_dataframe_BBK)/3)
    
  # Correct the temperatures for thermal inertia
  print("Correcting the temperature for thermal inertia")
  thermal_inertia_corrected_temp = thermal_inertia_correction(landcover_dataframe_BBK, inertiafiles_directory, degrees_of_freedom)
  landcover_dataframe_BBK$temperature = thermal_inertia_corrected_temp
  landcover_dataframe_BBK$thermal_inertia = thermal_inertia_corrected_temp
  
  # Data used for the detailed plot in 'Rivierenhof'
  landcover_dataframe_BBK_detailed = landcover_dataframe_BBK[1:105,]
    
  # Visualise the corrected and uncorrected temperatures on the BBK land cover map
  visualisation_route_measured_temp_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK, TRUE, FALSE)
  visualisation_route_thermal_inertia_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK, FALSE)
    
  visualisation_route_measured_temp_BBK_detailed(landcover_dataframe_BBK_detailed, dataframe_landcover_map_BBK_detailed, TRUE, FALSE)
  visualisation_route_thermal_inertia_BBK_detailed(landcover_dataframe_BBK_detailed, dataframe_landcover_map_BBK_detailed, FALSE)
  
  print("The temperature is corrected for thermal inertia")
}

# Apply a correction for temperature decline 
if(temp_decline){
  
  print("Temperature decline section")
  
  print("Correcting the temperatures for temperature decline based on the set of weather stations (based on the BBK land cover map of Flanders)")
  
  # Create a vector with uniformly increasing temperature residuals from begin to end of route (used in traditional decline corrections)
  temp_diff = landcover_dataframe_BBK$measured_temp[1] - landcover_dataframe_BBK$measured_temp[length(landcover_dataframe_BBK$temperature)]
  vector_residuals = seq(0, temp_diff, length.out = length(landcover_dataframe_BBK$temperature))
  
  # Calculate the temperature residuals and smoothed temperature residuals for the set of weather stations
  station_smoothed_residual_dataframe = station_smoothed_residuals(landcover_dataframe_BBK, station_path_list, station_name_list, 1)
  # Plot the temperature residuals and smoothed temperature residuals of the set of weather stations 
  residu_plot_stations(station_smoothed_residual_dataframe, station_name_list, vector_residuals)
  
  # Apply the temperature decline correction, temperatures are referred to times with a time interval given by 'seconds_per_interval' in the general settings with the begin time as starting point
  corr_temp_timeseries_BBK = temp_decline_timeseries(landcover_dataframe_BBK, seconds_per_interval, station_landcover_BBK, station_path_list, station_name_list, buffer_distances)
  
  # Create a plot of the temperature evolution with and without thermal inertia and temperature decline correction applied
  decline_corrected_temperatures_BBK = temp_decline_corr(landcover_dataframe_BBK, station_landcover_BBK, station_path_list, station_name_list, buffer_distances, 1, TRUE)
  decline_corr_temp_no_smoothing_BBK = temp_decline_corr(landcover_dataframe_BBK, station_landcover_BBK, station_path_list, station_name_list, buffer_distances, 1, FALSE)
  dataframe_temperatures = data.frame("date" = landcover_dataframe_BBK$fulldate, "decline_corrected" = decline_corrected_temperatures_BBK, "decline_corrected_no_smoothing" = decline_corr_temp_no_smoothing_BBK, "inertia_corrected" = landcover_dataframe_BBK$thermal_inertia, "measured" = landcover_dataframe_BBK$measured_temp)
  dataframe_temperatures = plot_temp_evolutions(dataframe_temperatures)
  
  print("The temperatures are corrected for temperature decline based on the set of weather stations")
  
  if (chart){
    
    # Specify the number of grid points in the horizontal and vertical directions
    number_horizontal_steps = ceil((max(landcover_dataframe_BBK$X) - min(landcover_dataframe_BBK$X) + (2*extension_of_chart))/step_distance)
    number_vertical_steps = ceil((max(landcover_dataframe_BBK$Y) - min(landcover_dataframe_BBK$Y) + (2*extension_of_chart))/step_distance)
      
    # Calculate the model variables for the grid points of the temperature chart
    # The land cover fractions of the grid points of the temperature chart are calculated based on the BBK land cover map of Flanders
    # The indices of the grid points where the validation weather stations are located are calculated as well
    chart_variables_loc_indices_BBK = chart_variables_and_locations(landcover_dataframe_BBK, number_horizontal_steps, number_vertical_steps, step_distance, extension_of_chart, buffer_distances_chart, validation_loc_lambert, validation_loc$Name, model_variables)
    chart_variables_BBK = chart_variables_loc_indices_BBK$chart_variables
    indices_of_locations = chart_variables_loc_indices_BBK$indices_of_locations
    indices_removed_loc = chart_variables_loc_indices_BBK$removed_indices
    
    print("Creating temperature charts")
    
    # Predict the temperatures of the grid points
    # The temperature measurements of the mobile bicycle campaign are used as input for the linear model together with the 'model_variables' selected in the general settings 
    # The grid point temperatures are calculated at all reference times for which decline corrected temperatures are available 
    charts_temp_data_BBK = chart_timeseries(landcover_dataframe_BBK, number_horizontal_steps, number_vertical_steps, model_variables, step_distance, corr_temp_timeseries_BBK, chart_variables_BBK)
    
    # Extract the names of the validation stations
    stations_validation_names = chart_variables_loc_indices_BBK$validation_names
    # Number the validation stations
    stations_validation_numbers = as.numeric(seq(1,length(stations_validation_names),by=1))
    # Create names of the stations that will be plotted on the temperature charts
    station_names_list = create_station_names(landcover_dataframe_BBK, seconds_per_interval, indices_removed_loc, stations_validation_numbers, validation_files)
    
    # Calculate the minimum and maximum grid point temperatures
    min_temp = min(unlist(charts_temp_data_BBK$chart_temperatures_list))
    max_temp = max(unlist(charts_temp_data_BBK$chart_temperatures_list))

    # Create the temperature charts
    create_temp_chart(charts_temp_data_BBK, number_horizontal_steps, number_vertical_steps, seconds_per_interval, landcover_dataframe_BBK, step_distance, min_temp, max_temp, indices_of_locations, station_names_list, "weather_stations")
    
    print("The temperature charts are created")
    
    print("Validating the weather station temperatures predicted based on the linear models used for the temperature charts")
    
    # Extract the observed temperatures of the validation weather stations at the reference times
    # Extract the predicted temperatures of the validation weather stations at the reference times based on the linear models used to create the temperature charts
    val_data = validation(validation_files, landcover_dataframe_BBK, indices_removed_loc, stations_validation_names, station_landcover_BBK, chart_variables_BBK, model_variables, corr_temp_timeseries_BBK, step_distance, buffer_distances)  
    
    # Calculate the RMSE, BIAS and MAE scores based on the observed and predicted temperatures of the validation weather stations
    val_scores = rmse_mae_bias(val_data)
    
    print(paste0("RMSE: ", val_scores$rmse, " BIAS: ", val_scores$bias, " MAE: ", val_scores$mae))
    print("The validation is done")
  }
  
  print("Visualising the temperatures of the route corrected for temperature decline at different times")
        
  # Visualise the uncorrected measured temperatures and the temperatures only corrected for thermal inertia on the BBK land cover map 
  visualisation_route_measured_temp_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK, TRUE, TRUE, corr_temp_timeseries_BBK)
  visualisation_route_thermal_inertia_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK, TRUE, corr_temp_timeseries_BBK)
        
  # Visualise the temperatures of the route corrected for temperature decline 
  # There is a visualisation at each time for which decline corrected temperatures are calculated
  visualisation_route_timeseries_BBK(landcover_dataframe_BBK, dataframe_landcover_map_BBK, seconds_per_interval, corr_temp_timeseries_BBK, "weather_stations")
        
  print("The visualisation is done")
}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#  

print("Program is done")


