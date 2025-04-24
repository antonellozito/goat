# Description
#------------
# Simple script to deform structures. Very basic, one has to implement
# ones own deformations etc manually. 
from . import Datahandler as dh 
from . import Plotter as pl
import numpy as np
from . import goat_types as gt
from . import StructureDeformer as sd 
import copy
import os 

# Load structure
structuredir = './goatf/Examples/JT60SA/Baseline/structure.dat'
writedir  = './goatf/Examples/JT60SA/5degreesRP'
topomeshdir = './goatf/Examples/JT60SA/Baseline/output/topomesh.dat'
structure = dh.ReadStructureFile(structuredir)
newstructure = copy.deepcopy(structure)
topomesh = dh.ReadTopomeshFile(topomeshdir)

pl.PlotStructure(structure, 0)

# Set operations (sequential!):
# 1: rotate point
# 2: translate point 
# 3: rotate entire vessel
# 4: translate entire vessel

# Operations to apply translations and rotations
ApplyDeformations = 2
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

# Operations to adjust original structure
else:
    operations = [2, 2, 2, 2, 2, 2, 2]
    structIDs = [7, 7, 7, 7, 7, 7, 7]
    points = [25, 26, 27, 28, 29, 30, 31]

    xdata = [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 1.0]
    ydata = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0]
    zdata = [0.07, 0.07, 0.07, 0.07, 0.12, 0.03, 0.05]




for i in np.arange(0, len(operations), 1):
    ts = structIDs[i]-1
    tp = points[i]-1
    tx = xdata[i]
    ty = ydata[i]
    tz = zdata[i]
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
    


# Visualize structure, with vertex/structure IDs 
pl.PlotStructure(structure, 1, color='r')
pl.PlotStructure(newstructure, 1, color='g')
pl.PlotTopologicalMesh(topomesh, 1)

# Write output
try:
    os.mkdir(writedir)
except:
    print("could not make directory")
dh.WriteStructureFile(writedir, newstructure)

pl.ShowFigures()
