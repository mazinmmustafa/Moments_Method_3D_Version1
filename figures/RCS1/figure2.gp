reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'figure2.tex'

set key top right

set xrange [0:180]
set yrange [-20:+30]

set xtics 30
set ytics 10

# set mxtics 20
# set mytics 10

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

# set arrow from 0.0, 0.0 to 20.0, 0.0 nohead linestyle 1 lc 'black' lw 2

set xlabel '$\theta$ [deg]'
set ylabel '$\sigma_{\varphi\varphi}$ [dB]'
# set title 'Bessel Functions'

plot 'figure2.txt' using 1:3 with lines lw 2 dt 1 lc 16 title ''