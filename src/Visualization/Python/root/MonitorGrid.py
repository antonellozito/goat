from src import Plotter as plotter
from src import Datahandler as dh

pausetime = 1
maxtime = 1000

# Assumed to be ran from run directory
datadir = dh.GetDataDirectory()

plotter.MonitorGrid(datadir, 1, pausetime, maxtime)
