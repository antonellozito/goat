from python_source import Plotter as plotter
from python_source import Datahandler as dh


pausetime = 1
maxtime = 1000

# Assumed to be ran from run directory
datadir = dh.GetDataDirectory()

plotter.MonitorGridAndVessel(datadir, 1, pausetime, maxtime)
