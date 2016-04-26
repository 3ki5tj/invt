#!/usr/bin/env gnuplot



# Error components for the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 3.0 font "Times, 24"
set output "xerr.eps"
set multiplot


wleft = 0.54
wright = 1 - wleft

# left panel
set size wleft, 1
set origin 0, 0

set rmargin 0.5
set bmargin 5

set xtics 5 offset 0, 0.3
set mxtics 5
#set xlabel "{/Symbol-Oblique s} (Gaussian) or {/Times-Italic n}/{/Times-Italic K}/{/Symbol \326}(2{/Symbol-Italic p}) (Bandpass)" offset 0, 0.5
set xlabel "Eigenmode" offset 0, 0.7
set xrange [0:17]

set logscale y
set format y "10^{/*0.8 %T}"
set ytics 100 offset 0.5, 0
#set format y "%.0t{/Symbol \264}10^{%T}"
#set mytics 2
set ylabel "Error component" offset 1, 0

a0 = 0.0001
fac = 2*100/sqrt(2*pi)

set title "Global" offset 0, -0.5

set key at screen 0.19, 0.1 Left reverse samplen 1.5 width -4 maxrows 1

plot [:][1e-11:1e-4] \
    "../../data/xerr/xerr_sig10_t1e8_g.dat"   u 1:($2) w lp lt 1 lw 1 pt 12                  t "Initial", \
    "../../data/xerr/xerr_sig10_t1e8_g.dat"   u 1:($3) w lp lt 1 lw 1 pt 13 lc rgb "#808080" notitle, \
    "../../data/xerr/xerr_sig10_t1e10_g.dat"  u 1:($3) w lp lt 1 lw 1 pt 13                  notitle, \
    "../../data/xerr/xerr_sinc_t1e8_g.dat"    u 1:($2) w lp lt 4 lw 1 pt  4 ps 0.8                  t "Initial", \
    "../../data/xerr/xerr_sinc_t1e8_g.dat"    u 1:($3) w lp lt 4 lw 1 pt  5 ps 0.8 lc rgb "#808080" notitle, \
    "../../data/xerr/xerr_sinc_t1e10_g.dat"   u 1:($3) w lp lt 4 lw 1 pt  5 ps 0.8                  notitle, \
    -1 notitle



# right panel
set origin wleft, 0
set size wright, 1

set lmargin 4

unset ylabel
#set format y ""

set title "Local"

set key at screen 0.98, 0.1 Left reverse samplen 1.5 width -2 maxrows 1

plot [:][1e-9:1] \
    "../../data/xerr/xerr_sig10_t1e8_l.dat"   u 1:($2) w lp lt 1 lw 1 pt 12                  notitle, \
    "../../data/xerr/xerr_sig10_t1e8_l.dat"   u 1:($3) w lp lt 1 lw 1 pt 13 lc rgb "#808080" t "{/Times-Italic T} = 10^{/*0.8 8}", \
    "../../data/xerr/xerr_sig10_t1e10_l.dat"  u 1:($3) w lp lt 1 lw 1 pt 13                  t "{/Times-Italic T} = 10^{/*0.8 10}", \
    "../../data/xerr/xerr_sinc_t1e8_l.dat"    u 1:($2) w lp lt 4 lw 1 pt  4 ps 0.8                  notitle, \
    "../../data/xerr/xerr_sinc_t1e8_l.dat"    u 1:($3) w lp lt 4 lw 1 pt  5 ps 0.8 lc rgb "#808080" notitle, \
    "../../data/xerr/xerr_sinc_t1e10_l.dat"   u 1:($3) w lp lt 4 lw 1 pt  5 ps 0.8                  notitle, \
    -1 notitle


unset multiplot
unset output
set terminal pop
reset
