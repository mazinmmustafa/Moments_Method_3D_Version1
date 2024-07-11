reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'figure1.tex'

set key top right outside

set polar
set angle radian
set size square

set tmargin 3
set bmargin 3

set grid polar ls 1 dt 1 lc 'black'

unset border
unset xtics
unset ytics
                                    
tickstep = 10                                   

r_min = -30
r_max = 0
set trange [-pi+pi/2:pi+pi/2]
set rrange [r_min:r_max]
set rtics tickstep

r_add = (r_max-r_min)/5
set size ratio -1

set label '$+90^{\circ}$' at (r_add+r_max-r_min)*cos(0*pi/180), (r_add+r_max-r_min)*sin(0*pi/180) center
set label '$+60^{\circ}$' at (r_add+r_max-r_min)*cos(30*pi/180), (r_add+r_max-r_min)*sin(30*pi/180) center 
set label '$+30^{\circ}$' at (r_add+r_max-r_min)*cos(60*pi/180), (r_add+r_max-r_min)*sin(60*pi/180) center
set label '$0^{\circ}$' at (r_add+r_max-r_min)*cos(90*pi/180), (r_add+r_max-r_min)*sin(90*pi/180) center
set label '$-30^{\circ}$' at (r_add+r_max-r_min)*cos(120*pi/180), (r_add+r_max-r_min)*sin(120*pi/180) center 
set label '$-60^{\circ}$' at (r_add+r_max-r_min)*cos(150*pi/180), (r_add+r_max-r_min)*sin(150*pi/180) center
set label '$-90^{\circ}$' at (r_add+r_max-r_min)*cos(180*pi/180), (r_add+r_max-r_min)*sin(180*pi/180) center
set label '$+120^{\circ}$' at (r_add+r_max-r_min)*cos(-30*pi/180), (r_add+r_max-r_min)*sin(-30*pi/180) center 
set label '$+150^{\circ}$' at (r_add+r_max-r_min)*cos(-60*pi/180), (r_add+r_max-r_min)*sin(-60*pi/180) center
set label '$\pm 180^{\circ}$' at (r_add+r_max-r_min)*cos(-90*pi/180), (r_add+r_max-r_min)*sin(-90*pi/180) center
set label '$-150^{\circ}$' at (r_add+r_max-r_min)*cos(-120*pi/180), (r_add+r_max-r_min)*sin(-120*pi/180) center 
set label '$-120^{\circ}$' at (r_add+r_max-r_min)*cos(-150*pi/180), (r_add+r_max-r_min)*sin(-150*pi/180) center

plot 'figure1.txt' using (pi/2-$1):($2 > r_min ? $2 : r_min) with lines lw 2 lc 16 dt 1 title '$20\log_{10}|E_{\theta}|$ [dB]'
