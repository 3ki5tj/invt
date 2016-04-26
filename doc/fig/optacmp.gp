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
set format x "10^{/*0.8 %T}"
set mxtics 10
set xrange [1e4:1e8]
set xlabel "Shifted simulation time, {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0}"

set logscale y
set format y "10^{/*0.8 %T}"
set mytics 10
set yrange [3e-9:1e-4]
set ylabel "{/Symbol-Oblique a} ({/Times-Italic t})"


set key at 4.5e7, 14e-7 Left reverse font "Times, 22" spacing 1.2

alpha0 = 1e-4
t0 = 2 / alpha0

plot [:][:] \
    "../../data/opta/sig10_g_alpha_q.dat"      u ($1+t0):($2) w l lt 1 lw 5 lc rgb "#000000"    t "Global", \
    "../../data/opta/sig10_l_alpha_q.dat"      u ($1+t0):($2) w l lt 1 lw 2 lc rgb "#000000"    t "Local", \
    "../../data/opta/sig10k4_g_alpha_q.dat"    u ($1+t0):($2) w l lt 2 lw 5 lc rgb "#000000"    t "Global, modified", \
    "../../data/opta/sig10k4_l_alpha_q.dat"    u ($1+t0):($2) w l lt 2 lw 2 lc rgb "#000000"    t "Local, modified", \
    1/x lc rgb "#606060" t "1 / ( {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0 })", \
    -1 notitle




unset output
set terminal pop
reset
