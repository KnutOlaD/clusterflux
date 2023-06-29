Script for doing flare clustering of single beam echosounder 
flare data

Contact info: 
Knut Ola Dølven
knut.o.dolven@uit.no

######################
##### WHAT IS IT #####
######################

Clusters flare observations according to method described in 
Veloso et al. (2015), doi: 10.1002/lom3.10024 , with some 
added features. Sorts all the flares according to clusters. 
Returns an excel sheet with the flowrate, area, and depth of 
the clusters as well as a .txt file listing the name of flare 
observations which were included in clusters and which flare 
observations were non-clustered (lonely flareas). 


########################################
##### CONTENTS AND HOW TO INITIATE #####
########################################

Includes:

-----------------------------------------------
clustering.py:
-----------------------------------------------
Author: Knut Ola Dølven, knut.o.dolven@uit.no
License: MIT

The main file which contains a FUNCTIONS section and 
an INITIATION section for the clustering program. 
The script starts with the FUNCTIONS section and the 
INITIATION section is at the end of the script. The
script is set up to work with output excel files from ESP3
with the following headers: 

Field_Name,Average_Lat_C_Foot,Average_Lon_C_Foot,Average_X_C_Foot
Average_Y_C_Foot,UTM_zone,Average_Radius_Foot,Flow_Rate_realBRS,

but this can easily be changed within the code and 
support for various header and filetypes is easy to implement
and will be implemented soon. 


FUNCTIONS:
load_flare_data: Loads the FlareHunter/ESP3 excel 
dataset. The excel file contains data from the
FlareHunter/ESP3 experiment. The data is stored 
in a pandas dataframe. The function takes a filepath 
and a column identified input for the excel file
which defines the column names in the excel file. 
It is also possible to define the UTM zone explicitly 
instead of loading it into the dataframe from the excel file.

get_distance: Calculates the distance between the points 
in a vector with coordinates xlocs and ylocs.
The function returns a matrix with the distances
between all the points in the vector.

get_shared_area: Calculates the shared area between two
 circles with radius R1 and R2 and mto different origins 
located at x_loc and y_loc. The function returns the shared area.
Reference: Weisstein, Eric W. "Circle-Circle Intersection." 
From MathWorld--A Wolfram Web Resource. 
http://mathworld.wolfram.com/Circle-CircleIntersection.html

get_area_of_footprint: Calculates the area of the footprint
 of the acoustic beam of the SB echosounder at a specified 
depth. The function returns the area of the footprint.

get_close_flares: Finds all the flares that are within 
a certain distance of each other or shares a certain
amount of area. The function use UTM x an y or lan/lot
 coordinates and footprint radius to return a list of pairs 
with indices of flares that are defined as close to each 
other according to the closeness_param and threshold 
parameters. Can calculate closeness based on distance between 
footprint center or the fractional (%) shared area between 
the footprints.

cluster_flares_mario: Clusters flares that have overlapping 
areas/are within a certain distance of each other (this is 
done through get_close_flares) and calculates the total area 
and flowrate of each cluster. The function returns a dictionary 
with the cluster info. Based on method described in
Veloso et al., 2015, doi: 10.1002/lom3.10024

save_clustered_data: Creates a dataframe with clustered 
flares and their locations and flowrates together with
the flares that were not clustered as individual clusters. 
Saves the dataframe to an excel file with path/filename filepath.

write_cluster_textfile: Writes a textfile listing all the 
flare observations included in each cluster as well as all 
the flare observations that were not clustered.


INITIATION: 
There are two initation options which is set by the 
beginning of the INITIATION section in the script by 
choosing the "runGUI" parameter to runGUI == False or 
runGUI == True . The former makes the script run dynamically
as a normal python script while the latter returns a simple
gui for giving the input parameters. 

Inputs: 
filepath: string
	Path and filename of the excel file to be stored. 
closeness_parameter: string
	The preferred method for determining wether two flares
	should be clustered or not. Options are 'distance' or 
	'area'. The former uses the distance between the center of 
	the flare observations as clustering criterion while the 
	latter uses the amount of shared area between two flare 
	observations
threshold: The threshold tha determines if two flare observations are to 
	be clustered. If closeness_param = 'area' the threshold is the 
	fractional (between 0 and 1) overlap between the flares. If 
	closeness_param = 'distance' the threshold is the distance between '
	the flares number of flare footprint radii.

Ouputs: 
DFclustered: pandas dataframe
	Dataframe containing the clustered flare data.

Creates an .xslx file with clustered flare data named "filename + clustered.xlsx"
and stores it in the same folder as filename. Creates a plot if plot=True. 
Creates a .txt file with the names of the flares in each cluster and the cluster 
name and the names of the flares that were not clustered.

IMPORTANT NOTES: 

1. In the GUI the strings should not have quotation marks 
2. The path and filename of the excel file needs to have double 
backslashe's if you're running windows, i.e. the path C:\Users\
should be written C:\\Users\\.
--------------------------------------------------------------


--------------------------------------------------------------
conversion.py: 
--------------------------------------------------------------
Copied from the utm package by Tobias Bieniek
Bidirectional UTM-WGS84 converter for python
Author: Tobias Bieniek
Version: 0.7.0
License: MIT
Source: https://pypi.org/project/utm/

Copied from the python utm package to avoid any
dependency issues.


 
