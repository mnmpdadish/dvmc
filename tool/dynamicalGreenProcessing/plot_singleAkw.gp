#!/usr/bin/gnuplot 


w = 8.5
h = 12
set terminal pdf enhanced color size w cm,h cm
set output 'Akw.pdf'
#set terminal epslatex size w cm,h cm color
#set output 'plot.tex'

set pm3d map interpolate 0,0;
set macros
set multiplot

set xrange [0:32]
set yrange [0:16]
set cbrange [0.00:0.7]

top = 0.95
left = 0.10
right = 0.85
bottom = 0.05

set palette defined (0 'white',0.02 'white',0.06 '#ffea4a',0.1 '#fdb416',0.18 '#fa3100',0.27 '#c31b00',0.6 'black')

set colorbox vertical user origin right+0.02, top-0.4 size .04,.4
offset1 = 5.7
set ytics ('-4' offset1-4, '0' offset1, '4' offset1+4.0, '8' offset1+8.0, '12' offset1+12.0)
set tics front
set arrow from 0,offset1 to 32,offset1 nohead lw 2 dt 2 lc 'black' front
#set xtics ('$0$' 0, '$\pi$' 32 )
set xtics ('0' 0, '{/Symbol p}' 32)
set lmargin at screen left; set rmargin at screen right
set tmargin at screen top; set bmargin at screen bottom
plot 'output/Akw.dat' u 1:($2*16.0/1500.0):3 matrix notitle w image 

