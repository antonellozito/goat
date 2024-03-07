import Plotter as plotter

pausetime = 1
maxtime = 1000

datadir = './goatf/Examples/ASDEX_Sander'

plotter.MonitorGrid(datadir, 1, pausetime, maxtime)
