flare_clustering
Script to do flare clustering of single beam echosounder flare data

Conversion.py is copied from the utm package by Tobias Bieniek
Bidirectional UTM-WGS84 converter for python
Author: Tobias Bieniek
Version: 0.7.0
License: MIT
Source: https://pypi.org/project/utm/

Everything else (except the MIT-license text) is written by me:
Knut Ola Dølven, knut.o.dolven@uit.no

######################
##### WHAT IS IT #####
######################

Clusters flare observations according to method described in Veloso et al. (2015), 
doi: 10.1002/lom3.10024 , with some added features. Sorts all the flares
according to clusters. Returns an excel sheet with the flowrate, area, and 
depth of the clusters as well as a .txt file listing the name of flare observations 
which were included in clusters and which flare observations were non-clustered 
(lonely flareas). 


####################
##### CONTENTS #####
####################

Includes:

clustering.py - the main script 

 
