set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "Plot2DPolygon"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="./src/Visualization/quiverplot2d.dat"
plot m using 1:2:3:4 with vectors