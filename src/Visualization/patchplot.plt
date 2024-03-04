set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "Patchplot"
set nokey
set grid
set xlabel "x"
set ylabel "y"
set size ratio -1
m="./src/Visualization/patchplot.dat"
plot m using 1:2:3 lc pal 
