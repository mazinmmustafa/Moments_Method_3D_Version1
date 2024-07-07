reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'figure1.tex'

set key top right

# set xrange [0:180]
# set yrange [-20:+30]

# set xtics 30
# set ytics 10

# set mxtics 20
# set mytics 10

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

# set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc 'black' lw 2

set xlabel '$z$ [m]'
set ylabel '$|E|$ [V/m]'
# set title 'Bessel Functions'

plot 'figure1.txt' using 1:2 with lines lw 2 dt 1 lc 16 title '$E_{x}$',\
     'figure1.txt' using 1:3 with lines lw 2 dt 2 lc 16 title '$E_{y}$',\
     'figure1.txt' using 1:4 with lines lw 2 dt 4 lc 16 title '$E_{z}$' 