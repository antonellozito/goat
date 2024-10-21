import Plotter as plotter
import Datahandler as dh

pausetime = 1
maxtime = 1000

# Assumed to be ran from run directory
datadir = dh.GetDataDirectory()

plotter.MonitorGrid(datadir, 1, pausetime, maxtime)
