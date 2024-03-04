set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "PlotGrid"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="./src/Visualization/plotgrid.dat"
plot m using 1:2 with points