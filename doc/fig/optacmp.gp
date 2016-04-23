#!/usr/bin/env gnuplot



# optimized schedule versus the optimal inverse-time schedule
# for Gaussian updating scheme with sigma = 10



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 2.7 font "Times, 24"
set output "optacmp.eps"


reset


htop = 0.5
hbot = 1 - htop

dx = 0.01
dy = 0.05

set lmargin 7

set logscale x
#set xtics 0.2 offset 0, 0.2
#set mxtics 2
set format x "10^{/*0.7 %T}"
set mxtics 10
set xrange [1e4:1e8]
set xlabel "Shifted simulation time, {/Times-Italic t} + 2/{/Symbol-Oblique a}_{/*0.7 0}"

set logscale y
set format y "10^{/*0.7 %T}"
set mytics 10
set yrange [3e-9:1e-4]
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})"


set key at 4e7, 2e-7 Left reverse font "Times, 22"

alpha0 = 1e-4
t0 = 2 / alpha0

plot [:][:] \
    "../../data/opta/sig10_g_alpha_q.dat"      u ($1+t0):($2) w l lt 1 lw 5 lc rgb "#000000"    t "Gaussian, global", \
    "../../data/opta/sig10_l_alpha_q.dat"      u ($1+t0):($2) w l lt 1 lw 2 lc rgb "#000000"    t "Gaussian, local", \
    "../../data/opta/sinc_g_alpha_q.dat"       u ($1+t0):($2) w l lt 2 lw 7 lc rgb "#000000"    t "Bandpass, global", \
    "../../data/opta/sinc_l_alpha_q.dat"       u ($1+t0):($2) w l lt 4 lw 2 lc rgb "#888888"    t "Bandpass, local", \
    -1 notitle




unset output
set terminal pop
reset
