set terminal png size 500, 350 font 'Verdana, 10'

set xzeroaxis lt -1
set yzeroaxis lt -1

set grid xtics lc rgb  '#555555' lw 1 lt 0
set grid ytics lc rgb  '#555555' lw 1 lt 0

set xtics axis
set ytics axis

set xrange[0:1]
set yrange[0:3]

set output 'Splines.png'
plot 'Splines.txt' using 1:2 with linespoints lw 1 lt rgb 'blue'

set output 'Runge_Kutt.png'
plot 'Runge_Kutt.txt' using 1:2 lt rgb 'red'

set output 'Diffirent.png'
plot 'Splines.txt' using 1:2 lt rgb 'blue', \
    'Runge_Kutt.txt' using 1:2 lt rgb 'red'