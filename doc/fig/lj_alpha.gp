#!/usr/bin/env gnuplot



# optimized schedule for the Lennard-Jones system



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced font "Helvetica, 32"
set output "lj_alpha.eps"


reset


htop = 0.5
hbot = 1 - htop

dx = 0.01
dy = 0.05

set lmargin 7

set logscale x
#set xtics 0.2 offset 0, 0.2
#set mxtics 2
set format x "10^{/*0.8 %T}"
set mxtics 10
set xrange [1e4:1e8]
set xlabel "Shifted simulation time, {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0}"

set logscale y
set format y "10^{/*0.8 %T}"
set mytics 10
set yrange [6e-9:1e-4]
set ylabel "{/Symbol-Oblique a}({/Times-Italic t})"


set key left bottom Left reverse spacing 1.5

alpha0 = 1e-4
t0 = 2 / alpha0

plot [:][:] \
    "../../data/lj/rho0.1/alpha_sig0.2_t1e8.dat"      u ($1+t0):($2) w l lt 1 lw 5 t "Optimal", \
    1/x lt 2 lw 2 t "1 / ({/Times-Italic t} + {/Times-Italic t}_{/*0.8 0})", \
    -1 notitle

unset output
set terminal pop
reset
