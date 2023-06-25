'''
Functions and initation for clustering techniques of flare data

contains:

load_flare_data(filenam): loads the flowrates and areas from the file filename into a 
dataframe.




'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


################################################
################# FUNCTIONS ####################
################################################

def load_flare_data(filepath):
    '''
    Function that loads the FlareHunter/ESP3 excel dataset. The excel file contains data from the 
    FlareHunter/ESP3 experiment. The data is stored in a pandas dataframe. This is a 8 column
    data file with rows equal to the number of datasets. The columns contain the following 
    information:
    [flare_identifier,Lat,Lon,UTM_X,UTM_Y,UTM_Zone,Radius,Flowrate]

    Input:
    dirpath: string.
        The path to the directory where the excel file is stored.

    Output: 
    DFdata: pandas dataframe
        Containes data for all the 8 columns in row "inputrow"
    '''

    #Load the data
    DFdata = pd.read_excel(filepath,header = 0, index_col = 0, usecols = [0,1,2,3,4,5,6,7])

    return DFdata

###############################################################################################

def shared_area(x_loc1,xloc2,y_loc1,yloc2,R1,R2):
    '''
    Function that calculates the shared area between two circles with radius R1 and R2 and
    centers at x_loc and y_loc. The function returns the shared area in m^2.

    Input:
    x_loc1: float
        The x location of the center of circle 1.
    x_loc2: float
        The x location of the center of circle 2.
    y_loc1: float
        The y location of the center of circle1.
    y_loc2: float
        The y location of the center of circle2.
    R1: float
        The radius of the first circle.
    R2: float
        The radius of the second circle.

    Output:
    shared_area: float
        The shared area between the two circles in m^2.
    '''

    #Calculate the distance between the two circles
    dist = np.sqrt((x_loc1-x_loc2)**2+(y_loc1-y_loc2)**2)

    #Calculate the shared area between the two circles
    if dist < R1+R2:
        #Calculate the shared area
        

        return shared_area




################################################
################# INITIATION ###################
################################################

if __name__ == '__main__':
    
    #Load the data
    filepath = 'C:\\Users\\kdo000\\Dropbox\\post_doc\\Marie_project\\data\\CAGE_18_02_FlareHunt-D20180523-T080503_4FR.xlsx'

    DFdata = load_flare_data(filepath)

    #Load the UTM_X and UTM_Y data into a numpy array

    UTM_X = DFdata['Average_X_C_Foot'].values
    UTM_Y = DFdata['Average_Y_C_Foot'].values

    #Calculate the distance between all the flares and subtract the radius of the flare area
    #and divide by the combined flare area of the two flares.

    #Create a matrix with the distances between all the flares in the dataset
    #indices follows the field_name order in the excel file

    flare_dists = np.zeros((len(UTM_X),len(UTM_X)))

    for i in range(len(UTM_X)):
        for j in range(len(UTM_X)):
            flare_dists[i,j] = np.sqrt((UTM_X[i]-UTM_X[j])**2+(UTM_Y[i]-UTM_Y[j])**2)  

    #Create a matrix with the combined flare areas of all the flares in the dataset

    flare_areas = np.zeros((len(UTM_X),len(UTM_X)))


