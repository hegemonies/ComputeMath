set terminal png size 500, 350 font 'Verdana, 10'
set output 't.png'
plot 't.txt' using 1:2 with linespoints lw 1 lt rgb 'purple'