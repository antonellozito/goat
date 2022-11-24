# Set plotting window
set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "PlotGridCells"
set nokey
set grid
set xlabel "x"
set ylabel "y"

# First, plot the grid, then points
m = "./src/Visualization/plotgridcells.dat"
n = "./src/Visualization/plotgridwithpoints.dat"
plot m using 1:2 with lines, n using 1:2 with points pt 7

