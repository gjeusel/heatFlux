#!/bin/bash

tools_contrib_path=$( cd "$( dirname "$0" )" && pwd )

if [ $# -lt 1 ] ; then
  echo "USAGE : bash $tools_contrib_path/$0 arg=path_to_Temperature.dat"
  exit
fi

input_file_name=$1

# ${string/%substring/replacement}
# If $substring matches back end of $string, substitute $replacement for $substring.
output_file_name=${input_file_name/%csv/png}

gnuplot -persist << EOF
fontname = 'Helvetica,12'
set terminal pngcairo dashed enhanced color font fontname size 1200,800

set title "Temperature Map"

# set grid to help reading :
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"

# Set colors rainbow :
set pm3d map
#set pm3d interpolate 0,0
set palette rgbformulae 35,13,10

# Input file contains comma-separated values fields
set datafile separator comma

# Set transparent fill for colors :
set style fill transparent solid 0.70 noborder

# Labels :
set xlabel 'x [m]'
set ylabel 'y [m]'
set cblabel "Temperature"

# Ranges
set xrange []
set yrange []

set output '$output_file_name'

# splot is the command for drawing 3D plots (well, actually projections on a 2D surface)
splot '$1' using 1:2:3 notitle

unset output
unset terminal
EOF

display $output_file_name
