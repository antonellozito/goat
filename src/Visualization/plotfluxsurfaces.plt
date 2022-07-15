set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "PlotFluxSurfaces"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="./src/Visualization/plotfluxsurfaces.dat"
plot m using 1:2 with points
plot m using 1:2:3 with labels 