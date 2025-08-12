# Description
#------------
# Simple script to deform structures. Very basic, one has to implement
# ones own deformations etc manually. 
from GOATpy import dh as dh 
from GOATpy import pl as pl
import numpy as np
from GOATpy import gt as gt
from python_source import StructureDeformer as sd 
import copy
import os 

# Load structure
root = '/mnt/c/Users/u0110555/Desktop/code_werk/goatf'
structuredir = root + '/goatf/Runs/JT60SA/Giulio/structure.dat'
writedir  = root + '/goatf/Runs/JT60SA/10degreesBot'
topomeshdir = root + '/goatf/Runs/JT60SA/Giulio/output/topomesh.dat'
structure = dh.ReadStructureFile(structuredir)
newstructure = copy.deepcopy(structure)
topomesh = dh.ReadTopomeshFile(topomeshdir)

pl.PlotStructure(structure, 0)

# Set operations (sequential!):
# 1: rotate point
# 2: translate point 
# 3: rotate entire vessel
# 4: translate entire vessel
# 5: add point before or after point  
# 6: set point coordinates equal to intersection of two lines between two structure edges
# 7: translate along edge vector
# 8: set point coordinates equal to other point coordinates

# Operations to apply translations and rotations
ApplyDeformations = 3
if ApplyDeformations == 1:
    # 6cm shift of entire vessel, then rotation
    operations = [4, 1, 1]
    structIDs = [0, 2, 2]
    points = [0, 3, 4]

    xdata = [-0.5, structure[1].x[4], structure[1].x[4]]
    ydata = [-0.3, structure[1].y[4], structure[1].y[4]]
    zdata = [0.06, 10, 10]

elif ApplyDeformations == 2:
    # Extension along limiter direction
    operations = [2, 1, 1]
    structIDs = [2, 2, 2]
    points = [2, 3, 4]

    dx = structure[1].x[1] - structure[1].x[0]
    dy = structure[1].y[1] - structure[1].y[0]
    dn = np.sqrt(dx**2 + dy**2)
    dx = dx/dn
    dy = dy/dn 
    xdata = [dx, structure[1].x[4], structure[1].x[4]]
    ydata = [dy, structure[1].y[4], structure[1].y[4]]
    #zdata = [0.12, 10, 10]
    zdata = [0.05, 5, 5]

elif ApplyDeformations == 3:
    # Rotate at bottom point, find new point on baffle as intersection, fix other point at intersection coordinates, translate point 2 to be below point 1 again
    operations = [1, 1, 1, 6, 8, 8, 7, 7, 8]
    structIDs = [10, 10, 10, 10, 21, 21, 21, 10, 14]
    points = [3, 4, 5, 5, 1, 2, 1, 5, 1]

    auxstructID1 = [0, 0, 0, 10, 10, 10, 10, 10, 21]
    auxstructID2 = [0, 0, 0, 14, 0, 0, 0, 0, 0]
    auxpoint1 = [0, 0, 0, 2, 5, 5, 4, 4, 2]
    auxpoint2 = [0, 0, 0, 1, 0, 0, 0, 0, 0]

    xdata = [structure[9].x[1], structure[9].x[1], structure[9].x[1], 0, 
        0, 0, 0, 0, 0]
    ydata = [structure[9].y[1], structure[9].y[1], structure[9].y[1], 0, 
        0, 0, 0, 0, 0]
    tempdist = np.sqrt((structure[20].x[0] - structure[20].x[1])**2 + (structure[20].y[0] - structure[20].y[1])**2)
    zdata = [10, 10, 10, 0, 0, 0, -tempdist, -tempdist, 0]



# Operations to adjust original structure
elif ApplyDeformations == 0:
    operations = [2, 2, 2, 2, 2, 2, 2]
    structIDs = [7, 7, 7, 7, 7, 7, 7]
    points = [25, 26, 27, 28, 29, 30, 31]

    xdata = [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 1.0]
    ydata = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0]
    zdata = [0.07, 0.07, 0.07, 0.07, 0.12, 0.03, 0.05]
else:
    # Dry run, do nothing
    operations = []




for i in np.arange(0, len(operations), 1):
    ts = structIDs[i]-1
    tp = points[i]-1
    tx = xdata[i]
    ty = ydata[i]
    tz = zdata[i]
    ts1 = auxstructID1[i]-1 
    ts2 = auxstructID2[i]-1 
    pID1 = auxpoint1[i]-1 
    pID2 = auxpoint2[i]-1
    if operations[i] == 1:
        newx, newy = sd.RotatePoint(newstructure[ts].x[tp], newstructure[ts].y[tp], 
            tx, ty, tz)
        newstructure[ts].x[tp] = newx 
        newstructure[ts].y[tp] = newy
    elif operations[i] == 2:
        newx, newy = sd.TranslatePoint(newstructure[ts].x[tp], newstructure[ts].y[tp], 
            tx, ty, tz)
        newstructure[ts].x[tp] = newx 
        newstructure[ts].y[tp] = newy   
    elif operations[i] == 3:
        for i in newstructure:
            newx, newy = sd.RotatePoint(i.x, i.y, tx, ty, tz)
            i.x = newx 
            i.y = newy
    elif operations[i] == 4:
        for i in newstructure:
            newx, newy = sd.TranslatePoint(i.x, i.y, tx, ty, tz)
            i.x = newx 
            i.y = newy
    elif operations[i] == 5: 
        newx, newy = sd.AddPoint(newstructure[ts].x, newstructure[ts].y, 
            tx, ty, tz, tp)
        newstructure[ts].x = newx 
        newstructure[ts].y = newy
        if newstructure[ts].n > 0:
            newstructure[ts].n = newstructure[ts].n + 1
        else:
            newstructure[ts].n = newstructure[ts].n - 1
    elif operations[i] == 6:
        # Intersection between two edges forms new point
        x11 = newstructure[ts1].x[pID1]
        x12 = newstructure[ts1].x[pID1+1]
        y11 = newstructure[ts1].y[pID1]
        y12 = newstructure[ts1].y[pID1+1]
        x21 = newstructure[ts2].x[pID2]
        x22 = newstructure[ts2].x[pID2+1]
        y21 = newstructure[ts2].y[pID2]
        y22 = newstructure[ts2].y[pID2+1]
        newx, newy = sd.LineIntersections(x11, y11, x12, y12, x21, y21, x22, y22)
        newstructure[ts].x[tp] = newx 
        newstructure[ts].y[tp] = newy
    elif operations[i] == 7:
        # Translate along vector
        tx = newstructure[ts1].x[pID1+1] -  newstructure[ts1].x[pID1]
        ty = newstructure[ts1].y[pID1+1] -  newstructure[ts1].y[pID1]
        newx, newy = sd.TranslatePoint(newstructure[ts].x[tp], newstructure[ts].y[tp], 
            tx, ty, tz)
        newstructure[ts].x[tp] = newx 
        newstructure[ts].y[tp] = newy  
    elif operations[i] == 8:
        # Set vertex coordinates equal
        newstructure[ts].x[tp] = newstructure[ts1].x[pID1]
        newstructure[ts].y[tp] = newstructure[ts1].y[pID1]


    


# Visualize structure, with vertex/structure IDs 
pl.PlotStructure(structure, 1, color='r')
pl.PlotStructure(newstructure, 1, color='g')
#pl.PlotTopologicalMesh(topomesh, 1)

# Write output
try:
    os.mkdir(writedir)
except:
    print("could not make directory")
dh.WriteStructureFile(writedir, newstructure)

pl.ShowFigures()
