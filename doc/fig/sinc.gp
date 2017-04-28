#!/usr/bin/env gnuplot



# window kernel function and updating matrix
# for the bandpass updating schemes


set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3 make dashed line longer
set terminal postscript eps enhanced dl 3 size 3.7, 3.7 font "Times, 24"
set output "sinc.eps"
set multiplot


dx = 0.01
dy = 0.05

normd(x, sig) = exp(-x*x/2/sig/sig)/sqrt(2*pi)/sig;

# bottom panel
set size 1, 1
set origin 0, 0

set tmargin 0.5
set lmargin 7

set xtics offset 0, 0.2
set mxtics 10
set xlabel "{/Times-Italic i}" offset 0, 0.5

set ytics 0.05 offset 0.2, 0
set mytics 5
set ylabel "{/Symbol-Oblique m}_{/Times-Italic i}" offset 2, 0

set key right bottom Left reverse invert spacing 1.2 width -8 font "Times, 18"

plot [:100][-0.045:0.12] \
    0 lt 1 lw 0.1 lc rgb "#cccccc" notitle, \
    normd(x, 10) lt 1 lw 5 lc rgb "#aaaaaa" t "Gaussian, {/Symbol-Oblique s} = 10", \
    "../../data/sinc/sinc_nonpbc_win.dat" u 1:($2) w lp pt 4 ps 0.8 lt 1 lw 1.0 t "Bandpass, non{/*0.7 -}periodic, {/Times-Italic K} = 10", \
    "../../data/sinc/sinc_pbc_win.dat"    u 1:($2) w lp pt 5 ps 0.8 lt 4 lw 1.0 t "Bandpass, periodic, {/Times-Italic K} = 5"

unset xlabel
unset ylabel

# inset
set origin 0.37, 0.44
set size 0.6, 0.5
set size ratio 1

set title "{/Times-Italic w_{ij}}" offset 0, -0.7 font "Times, 20"
#unset colorbox

set palette defined (0 1 1 1, 1 0.3 0.3 0.3, 2 0 0 0)
#set palette gray negative

set cbrange [-0.046:0.21]
set cbtics 0.10 offset -0.8, 0 font "Times, 14"
set mcbtics 5

set xtics 20,20,100 offset 0, 0.6 font "Times, 14"
#set mxtics 2
#set format x ""
#set xlabel "{/Times-Italic*0.7 j}" offset 0, 1008

set ytics 20,20,100 offset 0.8, 0 font "Times, 14"
#set mytics 2
#set format y ""
#set ylabel "{/Times-Italic*0.7 i}" offset 3.2, 0

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

#set pm3d map
plot [:] "../../data/sinc/sinc_nonpbc_winmat.dat" matrix u ($1+1):($2+1):($3) with image notitle

#unset pm3d

unset multiplot
unset output
set terminal pop
reset
