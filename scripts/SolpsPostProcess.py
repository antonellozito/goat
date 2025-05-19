# Description
#------------
# This script visualizes the most important SOLPS 
# simulation quantities based on state files provided by the 
# simulation. No attempt at recomputation of derived quantities is made
# as this will depend on the version of SOLPS and the available 
# physics therein. 

from python_source import solps_types as st 
from python_source import Datahandler as dh
from python_source import Plotter as pl 
from matplotlib import pyplot as plt
import numpy as np 
import sys 
import copy

# Initialize
#-----------
# Set initial paths
statefile = 'b2fstati' # state
gridfile = 'traduit.out.b2us_gg' # grid
plasmafile = 'b2fplasmf'
#gridfile = 'traduit.out.b2us_norefxp2_gd2_correctedb'

# Read command line arguments
narg = len(sys.argv)
for i in range(1, narg):
    # First argument is state file
    if (i == 1):
        statefile = sys.argv[i]
    elif (i == 2):
        gridfile = sys.argv[i]
    elif (i == 3):
        plasmafile = sys.argv[i]

simdir = '/mnt/c/Users/u0110555/Desktop/code_werk/SOLPS/runs/DEMO_2025/fully_extended_DEMO_Donly/afn_adjustedvessel/'
griddir = '/mnt/c/Users/u0110555/Desktop/code_werk/SOLPS/runs/DEMO_2025/fully_extended_DEMO_Donly/afn_adjustedvessel/'
#simdir = '/mnt/c/Users/u0110555/Desktop/code_werk/SOLPS/runs/DEMO_2025/extended2_DEMO_Donly/afn/'
#griddir = '/mnt/c/Users/u0110555/Desktop/code_werk/SOLPS/runs/DEMO_2025/extended2_DEMO_Donly/baserun/'
simdir = '/mnt/c/Users/u0110555/Desktop/code_werk/SOLPS/runs/DEMO_2025/fully_extended_DEMO_Donly/afn/'
griddir = '/mnt/c/Users/u0110555/Desktop/code_werk/SOLPS/runs/DEMO_2025/fully_extended_DEMO_Donly/baserun/'
statedir = simdir + statefile
plasmadir = simdir + plasmafile
griddir = griddir + gridfile


# Print out paths 
print('SolpsPostProcess: reading state file from: ' + statedir)
print('SolpsPostProcess: reading grid file from: ' + griddir)


# Initialize objects
state = st.PlasmaState()
const = st.PhysicalConstants()

# Load state and geometry
state.ReadB2fstatefile(statedir)
grid = dh.ReadTraduitOutB2us(griddir)

# Try loading the residuals
plotresiduals = True 
try:
    state.ReadB2fplasmfFile(plasmadir)
except:
    plotresiduals = False
    print("SolpsPostProcess: could not read b2fplasmf, not plotting residuals")

#state.WriteB2fstatefile(simdir + 'teststate')

# Figure counter
fignum = 0

# Visualize grid
#---------------
# Face labels
pl.PlotGridFaceLabels(grid, fignum)
pl.PlotGridFaces(grid, fignum)
# Set title and other descriptors
thisaxes = plt.gca()
thisaxes.set_title('Grid face labels')
thisaxes.set_xlabel('x [m]')
thisaxes.set_ylabel('y [m]')
thisaxes.legend(loc='upper right')
fignum = fignum + 1

# Visualize plasma state
#-----------------------
# Density
for i in range(0, state.ns):
    pl.PlotCellBasedQuantity2D(grid, np.log10(state.na[0:grid.cell.ntot, i]), fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('log10(density) [# m^-3], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

# Velocity
for i in range(0, state.ns):
    pl.PlotCellBasedQuantity2D(grid, state.ua[0:grid.cell.ntot, i], fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('Parallel velocity [m s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1 

# Temperature 
pl.PlotCellBasedQuantity2D(grid, state.ti[0:grid.cell.ntot]*const.eVperJoule, fignum)
thisaxes = plt.gca()
thisaxes.set_title('Ion temperature [eV], species: ' + str(i))
thisaxes.set_xlabel('x [m]')
thisaxes.set_ylabel('y [m]')
thisaxes.legend(loc='upper right')
fignum = fignum + 1
pl.PlotCellBasedQuantity2D(grid, state.te[0:grid.cell.ntot]*const.eVperJoule, fignum)
thisaxes = plt.gca()
thisaxes.set_title('Electron temperature [eV], species: ' + str(i))
thisaxes.set_xlabel('x [m]')
thisaxes.set_ylabel('y [m]')
thisaxes.legend(loc='upper right')
fignum = fignum + 1
pl.PlotCellBasedQuantity2D(grid, state.tn[0:grid.cell.ntot]*const.eVperJoule, fignum)
thisaxes = plt.gca()
thisaxes.set_title('Neutral temperature [eV], species: ' + str(i))
thisaxes.set_xlabel('x [m]')
thisaxes.set_ylabel('y [m]')
thisaxes.legend(loc='upper right')
fignum = fignum + 1

# Residuals
#----------
if plotresiduals:
    # Density
    for i in range(0, state.ns):
        pl.PlotCellBasedQuantity2D(grid, (state.resco[0:grid.cell.ntot, i]), fignum)
        thisaxes = plt.gca()
        thisaxes.set_title('Density residual [s^-1], species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 
    
    # Momentum
    for i in range(0, state.ns):
        pl.PlotCellBasedQuantity2D(grid, (state.resmo[0:grid.cell.ntot, i]), fignum)
        thisaxes = plt.gca()
        thisaxes.set_title('Momentum residual [s^-1], species: ' + str(i))
        thisaxes.set_xlabel('x [m]')
        thisaxes.set_ylabel('y [m]')
        thisaxes.legend(loc='upper right')
        fignum = fignum + 1 

    # Temperatures
    pl.PlotCellBasedQuantity2D(grid, state.reshi[0:grid.cell.ntot]*const.eVperJoule, fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('Ion temperature residual [s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1
    pl.PlotCellBasedQuantity2D(grid, state.reshe[0:grid.cell.ntot]*const.eVperJoule, fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('Electron temperature residual [s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1
    pl.PlotCellBasedQuantity2D(grid, state.reshn[0:grid.cell.ntot]*const.eVperJoule, fignum)
    thisaxes = plt.gca()
    thisaxes.set_title('Neutral temperature residual [s^-1], species: ' + str(i))
    thisaxes.set_xlabel('x [m]')
    thisaxes.set_ylabel('y [m]')
    thisaxes.legend(loc='upper right')
    fignum = fignum + 1

        
    
pl.ShowFigures()