from . import Plotter as plotter
from . import Datahandler as dh

pausetime = 1
maxtime = 1000

# Assumed to be ran from run directory
datadir = dh.GetDataDirectory()

plotter.MonitorGridAndVessel(datadir, 1, pausetime, maxtime)
