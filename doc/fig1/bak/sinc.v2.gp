#!/usr/bin/env gnuplot



# optimized schedule and the bandpass schemes
# for Gaussian updating scheme with sigma = 5



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4 make dashed line longer
set terminal postscript eps enhanced dl 4 size 3.7, 6 font "Times, 24"
set output "sinc.eps"
set multiplot


reset

htop = 0.45
hbot = 1 - htop

dx = 0.01
dy = 0.05

set label "(a)" at screen dx, 1 - dy    font "Times, 32"
set label "(b)" at screen dx, hbot - dy font "Times, 32"

normd(x, sig) = exp(-x*x/2/sig/sig)/sqrt(2*pi)/sig;

# bottom panel
set size 1, hbot
set origin 0, 0

set tmargin 0.5
set lmargin 7

set xtics offset 0, 0.2
set mxtics 10
set xlabel "{/Times-Italic i}" offset 0, 0.5

set ytics 0.05 offset 0.2, 0
set mytics 5
set ylabel "{/Symbol-Oblique m}_{/Times-Italic i}" offset 2, -1

set key right bottom Left reverse invert spacing 1.2 width -8

plot [:100][-0.06:0.1] \
    0 lt 4 lw 0.1 lc rgb "#cccccc" notitle, \
    normd(x, 6) lt 1 lw 5 lc rgb "#aaaaaa" t "Gaussian, {/Symbol-Oblique s} = 6", \
    "../../data/sinc/sinc_nonpbc_win.dat" u 1:($2) w lp pt 1 ps 1.0 lt 1 lw 1.0 t "Bandpass, non{/*0.7 -}periodic, {/Times-Italic K} = 6", \
    "../../data/sinc/sinc_pbc_win.dat"    u 1:($2) w lp pt 6        lt 1 lw 1.0 t "Bandpass, periodic, {/Times-Italic K} = 3"


# inset
set origin 0.38, hbot*0.51
set size 0.48, hbot*0.57


set title "{/Times-Italic w_{ij}}" offset 0, -0.7 font "Times, 20"
#unset colorbox
set cbrange [-0.05:0.2]
set cbtics 0.1 offset -1, 0 font "Times, 14"
set mcbtics 2

set xtics 20 offset 0, 0.8 font "Times, 14"
set mxtics 2
#set format x ""
set xlabel "{/Times-Italic*0.7 j}" offset 0, 1.8

set ytics 20 offset 1, 0 font "Times, 14"
set mytics 2
#set format y ""
set ylabel "{/Times-Italic*0.7 i}" offset 3.2, 0

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set pm3d map
splot [:] "../../data/sinc/sinc_nonpbc_winmat.dat" u 2:1:3 notitle

unset pm3d
reset



# top panel
set origin 0, hbot
set size 1, htop

set lmargin 7

unset logscale x
set xtics 0.2
set mxtics 2
#set format x "10^{%T}"
#set xrange [1e4:1e8]
set xlabel "Simulation time, {/Times-Italic t} ({/Symbol \264} 10^{/*0.7 8})" offset 0, 0

set logscale y
set format y "10^{%T}"
set mytics 10
#set yrange [1e4:1e9]
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})" offset 0, 0

set key right top Left reverse noinvert spacing 1.2 width -4

T = 1e8
plot [][7e-9:1e-4] \
    1/(x*T + 20000) lw 2 t "1 / ({/Times-Italic t} + 2/{/Symbol-Oblique a}_{/*0.7 0})", \
    "../../data/sinc/sinc_pbc_alpha.dat"    u ($1/T):($2) every 300 w p pt 6 ps 1.5 lw 2 t "Periodic, {/Times-Italic K} = 3", \
    "../../data/sinc/sinc_nonpbc_alpha.dat" u ($1/T):($2) every 300 w p pt 1 ps 1.5 lw 3 t "Non{/*0.7 -}periodic, {/Times-Italic K} = 6", \
    -1 notitle



unset multiplot
unset output
set terminal pop
reset
