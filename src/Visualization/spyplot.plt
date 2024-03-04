set term x11 font "-*-helvetica-medium-r-*-*-14-*-*-*-*-*-*-*" 
set title "Spyplot"
set nokey
set grid
set xlabel "x"
set ylabel "y"
m="./src/Visualization/spyplot.dat"
plot m using 1:2 with points