#!/bin/bash

tools_contrib_path=$( cd "$( dirname "$0" )" && pwd )

if [ $# -lt 1 ] ; then
  echo "USAGE : bash $tools_contrib_path/$0 arg=path_to_Temperature.dat"
  exit
fi

gnuplot -persist << EOF
fontname = 'Helvetica,12'
set terminal pngcairo dashed enhanced color font fontname size 1200,800

set title "Temperature Map"

# set grid to help reading :
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"

# Set colors rainbow :
set pm3d map
set pm3d interpolate 0,0
set palette rgbformulae 35,13,10

# Set transparent fill for colors :
set style fill transparent solid 0.70 noborder

# Labels :
set xlabel 'x [m]'
set ylabel 'y [m]'
set cblabel "Temperature"

# Ranges
set xrange []
set yrange []

set output 'Temperature.png'

# splot is the command for drawing 3D plots (well, actually projections on a 2D surface)
splot '$1' using 1:2:3 notitle

unset output
unset terminal
EOF

display Temperature.png
