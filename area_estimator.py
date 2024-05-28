'''
Script that estimates seepage area based on single beam echosounder data and the output
area from FlareHunter/ESP-3, being the accumulated acoustic footprint covered during 
a seep observation.

See bottom of the script for initiation. It works by loading a flarehunter/ESP3 output 
excel file and estimates the seep area based on the area of the acoustic footprint from
FlareHunter/ESP-3. It then outputs a new dataframe to an .xlsx file with the estimated
seep area instead of the accumulated acoustic footprint area (legacy area).

Uses functions from the clustering.py script.

Author: Knut Ola Dølven

FUNCTIONS:

load_flare_data(filepath,columns): 
Loads the FlareHunter/ESP3 excel dataset. The excel file contains output
data of flare observations from FlareHunter/ESP3. The data is stored in a pandas dataframe. 
The function takes a filepath and a column identified input for the excel file
which defines the column names in the excel file. It is also possible to define
the UTM zone explicitly instead of loading it into the dataframe from the excel file.

'''
#########################
### Importing modules ###
#########################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#add folder to aux_functions and clustering.py
import sys
sys.path.append(r'C:\Users\kdo000\Dropbox\post_doc\Marie_project\src\area_estimator')
sys.path.append(r'C:\Users\kdo000\Dropbox\post_doc\Marie_project\src\flare_clustering')

import aux_functions as aux


##############################    
######### FUNCTIONS ##########
##############################

def get_area_from_avg_radius(radius_array):
    '''
    Calculates the area of a circle based on the average radius of the circle.
    
    Input:
    radius_array: numpy array.

    Output:
    area_array: numpy array.
    '''

    area_array = np.pi*radius_array**2

    return area_array


########################################################################

def get_L_from_area(area_array,depth_array,opening_angle = 7):
    '''
    Calculates the length between the centers of the circles based on the area of the 
    footprint from FlareHunter/ESP-3

    Input:
    area_array: numpy array.

    depth_array: numpy array.

    opening_angle: float.

    Output:
    L_array: numpy array.
    
    '''

    #Calculate the radius of the acoustic footprint using the depth vector
    radius_afp = depth_array*np.tan(np.deg2rad(opening_angle/2))

    #Define the L_array
    L_array = np.zeros(len(area_array))

    #Calculate the length between the centers of the circles
    for i in range(len(area_array)):
        L_array[i] = area_array[i]/(2*radius_afp[i])-np.pi*radius_afp[i]/2
        
    #    (area_array[i]-np.pi*(depth_array[i]*np.tan(np.deg2rad(opening_angle/2))))\
    #    /(2*depth_array[i]*np.tan(np.deg2rad(opening_angle/2)))

    return L_array
    
########################################################################

def get_area_from_L(L,R,acc_fp_area,depth,opening_angle = 7,delta=None):
    '''
    Calculates the radius of the seep area based on  the length between the centers of the footprint
    circles. 

    Input:
    L: float.
        The length between the centers of the footprint circles.
    R: float.
        The radius of the footprint of the instrument
    acc_fp_area: float.
        The area of the accumulated seep specific footprint.
    depth_array: float.
        The depth of the seep area.
    opening_angle: float.
        The opening angle of the footprint circles.
    delta: float.
        The distance between the shiptrack and center of seep area center.
        Default: 0.5*R

    Output:
    A_seep: float.
        The area of the seep area.
    '''
    
    #Determine delta
    if (4+np.pi)*R**2 >= acc_fp_area:
        Delta = 0.5*R
    else:
        Delta = 0.5*(np.sqrt(R**2-0.25*L**2+R))
        
    A_seep = np.pi*(np.sqrt(0.25*L**2+Delta)-R)**2

    return A_seep

########################################################################

def get_area_from_r(r):
    '''
    Calculates the area of the seep area based on the radius of the seep area.

    Input:
    r: float.
        The radius of the seep area.

    Output:
    area: float.
        The area of the seep area.
    '''
    area = np.pi*r**2

    return area

########################################################################

def master_function_area_est(data,opening_angle = 7,delta=None):
    '''
    Master function for estimating seep area based on the output from FlareHunter/ESP-3.

    Input:
    data: pandas dataframe.
        The dataframe containing the data from FlareHunter/ESP-3.
    opening_angle: float.
        The opening angle of the footprint circles.
    delta: float.
        The distance between the shiptrack and center of seep area center.
        Default: 0.5*R

    Output:
    r: float.
        The radius of the seep area.
    '''

    #Get the radius of each flare observation and put them into an array 

    radius_array = data['Average_Radius_Foot'].values
    depth_array = data['Depth'].values

    area_array = get_area_from_avg_radius(radius_array)

    #Get the length between the centers of the circles

    L_array = get_L_from_area(area_array,depth_array,opening_angle)

    #Calculate the radius of the seep areas

    seep_area_array = np.zeros(len(L_array))
    seep_area_array = np.zeros(len(L_array))

    for i in range(len(L_array)):
        seep_area_array[i] = get_area_from_L(L_array[i],radius_array[i],area_array[i],depth_array[i],opening_angle,delta)
        #seep_area_array[i] = get_area_from_r(r_array[i])

    #Create a new column in the dataframe containing the seep radius and save 
    #the dataframe to a new excel file

    data['Seep_Radius_Est'] = np.sqrt(seep_area_array)/np.pi

    return data

########################################################################

########################
### Executing script ###
########################

if __name__ == '__main__':
    #load data
    filepath = r'C:\Users\kdo000\Dropbox\post_doc\Marie_project\src\flare_clustering\test_data\test_data.xlsx'

    #Load the data
    DFdata,varstrings = aux.load_flare_data(filepath)

    #use the masterfunction to estimate the seep radius
    DFdata = master_function_area_est(DFdata)

    #Replace the values in the Average_Radius_Foot column with the r_array values
    DFdata['Average_Radius_Foot'] = DFdata['Seep_Radius_Est']

    #output to output file. 
    DFdata.to_excel(r'C:\Users\kdo000\Dropbox\post_doc\Marie_project\src\flare_clustering\test_data\test_data_newest.xlsx')

    #display the dataframe
    
    



########################################################################

