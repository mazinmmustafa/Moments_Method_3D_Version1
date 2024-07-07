reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'figure3.tex'

# set palette rgb 7,5,15; set title "traditional pm3d\n(black-blue-red-yellow)"
# set palette rgb 3,11,6; set title "green-red-violet"
# set palette rgb 23,28,3; set title "ocean (green-blue-white)\ntry also other permutations"
# set palette rgb 21,22,23; set title "hot (black-red-yellow-white)"
# set palette rgb 30,31,32; set title "color printable on gray\n(black-blue-violet-yellow-white)"
# set palette rgb 33,13,10; set title "rainbow (blue-green-yellow-red)"
# set palette rgb 34,35,36; set title "AFM hot (black-red-yellow-white)"
set palette defined ( 0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000') ; set title "MATLAB (jet)"

set size ratio 1

set xrange [-1:1]
set yrange [-1:1]

set cbrange [0:+4]

# set xtics 2
# set ytics 2
# set cbtics 0.5
set tics front

# set arrow front from -4, 0 to +4, 0 nohead linestyle 1 lc 'black' lw 2
# set arrow front from 0, -4 to 0, +4 nohead linestyle 1 lc 'black' lw 2

set xlabel '$x$ [m]'
set ylabel '$z$ [m]'
set title '$|\mathcal{R}eH|$ [mA/m]'

plot 'figure3.txt' with image notitle
