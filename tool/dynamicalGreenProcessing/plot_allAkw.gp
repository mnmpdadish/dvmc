#!/usr/bin/gnuplot 


w = 2*8.5
h = 8
offset1 = 0.0


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

# A B C
top = 0.95
left = 0.09
right = 0.90
bottom = 0.14

margin1 = 0.05
width1 = (right-left-2.0*margin1)/3.0

rightA = left+width1
leftB = rightA+margin1
rightB = leftB+width1
leftC = rightB+margin1


set palette defined (0 'white',0.02 'white',0.06 '#ffea4a',0.1 '#fdb416',0.18 '#fa3100',0.27 '#c31b00',0.6 'black')

set colorbox vertical user origin right+0.02, (top+bottom)/2 size .02,.4
set ytics ('-12' offset1-12,'-8' offset1-8,'-4' offset1-4, '0' offset1, '4' offset1+4.0, '8' offset1+8.0, '12' offset1+12.0)
set tics front
set arrow from 0,offset1 to 32,offset1 nohead lw 2 dt 2 lc 'black' front
#set xtics ('$0$' 0, '$\pi$' 32 )
set xtics ('0' 0, '{/Symbol p}' 32)

set xlabel 'k' 
set ylabel '{/Symbol w}' 
set label 1 'A(k,{/Symbol w})' at screen left+0.4*width1, screen top+0.03 front
set label 2 'electron' at screen leftB+0.4*width1, screen top+0.03 front
set label 3 'hole' at screen leftC+0.4*width1, screen top+0.03 front

#set label 2 'salut' at screen 0.5, screen 0.5 front

#show label

set lmargin at screen left; set rmargin at screen rightA
set tmargin at screen top; set bmargin at screen bottom
plot 'output/Akw.dat' u 1:($2*16.0/1500.0):3 matrix notitle w image 

unset ylabel

set ytics ('' offset1-4, '' offset1, '' offset1+4.0, '' offset1+8.0, '' offset1+12.0)
set lmargin at screen leftB; set rmargin at screen rightB
set tmargin at screen top; set bmargin at screen bottom
plot 'output/Akw_e.dat' u 1:($2*16.0/1500.0):3 matrix notitle w image 

set lmargin at screen leftC; set rmargin at screen right
set tmargin at screen top; set bmargin at screen bottom
plot 'output/Akw_h.dat' u 1:($2*16.0/1500.0):3 matrix notitle w image 






