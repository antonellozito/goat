set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "PlotGridCells"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="./src/Visualization/plotgridcells.dat"
plot m using 1:2 with lines