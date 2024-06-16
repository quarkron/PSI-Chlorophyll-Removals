import numpy as np
import chimera
import os
import csv
from hamiltonianv2 import extract
from StructBio.Scenographics.solids import Solids
from StructBio.Utilities.miscFunctions import colorByName

csv_file_path = os.getcwd() + '\data\PSI_chls.csv'

def center(chl):
    mg = chl.findAtom('MG')
    return np.array(mg.coord())

def second_smallest_distinct(numbers):
    # Check if the list has at least two distinct elements
    if len(set(numbers)) < 2:
        return "List doesn't have a second distinct smallest number"

    # Initialize the smallest and second smallest variables
    smallest = second = float('inf')  # Set to positive infinity

    # Iterate through the list
    for num in numbers:
        if num < smallest:
            second = smallest
            smallest = num
        elif num < second and num != smallest:
            second = num

    return second

# runCommand("set bg_color white")

def interpolate_color(value, min_val, max_val):
    # Define the colors: blue, white, and red
    blue = (0, 0, 1)    # RGB for blue
    white = (1, 1, 1)   # RGB for white
    red = (1, 0, 0)     # RGB for red

    # Define the ranges where the gradient transitions occur
    blue_max = (min_val + max_val) / 2
    red_min = (min_val + max_val) / 2

    # Calculate the ratio of the value within the appropriate range
    if value <= blue_max:
        ratio = value / blue_max
        start_color = blue
        end_color = white
    else:
        ratio = (value - red_min) / (max_val - red_min)
        start_color = white
        end_color = red

    # Interpolate between start and end color based on the ratio
    r = (1 - ratio) * start_color[0] + ratio * end_color[0]
    g = (1 - ratio) * start_color[1] + ratio * end_color[1]
    b = (1 - ratio) * start_color[2] + ratio * end_color[2]

    return r, g, b

pdbID = raw_input("Type in four character PDB Identifier: \n")

p = chimera.openModels.open(pdbID, type = "PDB")[0]

cla_array = extract(p)

cla_dict = {}

with open(csv_file_path, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    first_row = True
    for row in csvreader:
        if first_row:
            first_row = False
            continue
        cla_dict[row[0]]=row[1::]

for chl in cla_array:
    chl_ID = chl.type + str(chl.id.position) + '.' + chl.id.chainId
    cla_dict[chl_ID].append(center(chl))


def construct_solids():
    BC_solids = Solids()
    inStrength_solids = Solids()
    outStrength_solids = Solids()
    farness_solids = Solids()
    transferP700_solids = Solids()
    euclid_distances_solids = Solids()
    solids_list = [inStrength_solids, outStrength_solids, farness_solids, BC_solids, transferP700_solids, euclid_distances_solids]

    special_pair = ['CLA1011.A', 'CLA1021.B']
    RC = ['CLA1012.B', 'CLA1013.A', 'CLA1022.A', 'CLA1023.B']
    

    for index_pos in [0,1,2,3,4,5]:
        # print(index_pos)

        #first calculate the range of values for a particular characteristic to determine proper scaling
        vals=[]
        for i in cla_dict.keys():
            if i not in special_pair:
                vals.append((float(cla_dict[i][index_pos])))


        min_xval = float(min(vals))
        max_xval = float(max(vals))
        # print((min_xval, max_xval))
                        
        if index_pos == 4: #log scale the transfer rate data
            max_xval = np.log(max_xval)+12
            min_xval = np.log(min_xval)+12
            print(min_xval, max_xval)
            
        for i in cla_dict:
            cen_i = cla_dict[i][-1]
            if i not in special_pair:
                val = float(cla_dict[i][index_pos])
                if index_pos == 4:
                    val = np.log(val)+12 
                color = interpolate_color(val, min_xval, max_xval)
                solids_list[index_pos].addSphere(cen_i, 1.25, color)
            if i in special_pair:
                solids_list[index_pos].addSphere(cen_i, 2, colorByName('black'))
        solids_list[index_pos].display()
        print('finished')

construct_solids()

#Index key: 1st = inStrength, 2nd = outStrength, 3rd= Farness, 4th = BC, 5th = P700transfer
        







    
    





