set term pngcairo transparent
set output "trjpc.png"

#set title "Free energy landscape"
set nokey

#set xlabel "PC1"
#set ylabel "PC2"  
#set ylabel "PC2"  offset -70,0
#set zlabel "PMF (kcal/mol)"

set xrange[-15:15]
set yrange[-10:10]
set view map

unset xtics
unset ytics

set label 3 point pt 4 ps 3  at 3.37282,-4.21373 
set label 1 point pt 2 ps 2  at 4.39508081616,   -2.52739873541
set label 2 point pt 1 ps 3  at 4.32654538931,   -2.69707799014

plot "out.dat" every 1000 w l  lc rgb "salmon"

