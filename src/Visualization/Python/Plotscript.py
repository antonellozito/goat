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
mylevels = [-1, -0.1, -0.01, -0.001, -0.0001, 0, 0.0001, 0.001, 0.01, 0.1, 1]
mylevelsbias = [0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05]
#pl.Plot2DSurfaceDataContourf(simdir + '/generalplf.dat', 0, levels=mylevels)
#pl.Plot2DSurfaceDataContourf(simdir + '/closedexactplf.dat', 1, levels=mylevels)
#pl.Plot2DSurfaceDataContourf(simdir + '/closedapproximationplf.dat', 2, levels=mylevels)
# pl.PlotPolygonData(simdir + '/testpolyg.dat', 2)
# pl.Plot2DSurfaceDataContourf(simdir + '/costfunctionLR_vesselcontours.dat', 2)
pl.Plot2DSurfaceDataContourf(simdir + '/costfunctionFA_weights.dat', 2)

pl.ShowFigures()
