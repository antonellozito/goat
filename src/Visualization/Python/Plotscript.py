# Description
#------------
# Simple plot script that uses the other python functionality
import Plotter as pl 

# Set paths
#----------
# Simulation directory
simdir = './goatf/Examples/ASDEX_Sander'

# Plot stuff
#-----------
mylevels = [-0.1, -0.01, -0.001, -0.0001, 0, 0.0001, 0.001, 0.01, 0.1]
pl.Plot2DSurfaceDataContourf(simdir + '/closedexactplf.dat', 1, levels=mylevels)
pl.ShowFigures()
