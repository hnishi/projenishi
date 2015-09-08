set term pngcairo transparent
set output "trjpc.png"

#set title "Free energy landscape"
#set nokey

#set xlabel "PC1"
#set ylabel "PC2"  
#set ylabel "PC2"  offset -70,0
#set zlabel "PMF (kcal/mol)"

set xrange[-15:15]
set yrange[-10:10]
set view map

set label 1 point pt 2 ps 3  at 3.37282,-4.21373 

plot "out.dat" every 100 w l

