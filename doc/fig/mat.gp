#!/usr/bin/env gnuplot



# updating matrix function



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced size 5, 4 font "Times, 28"
set output "mat.eps"


reset

htop = 0.4
hbot = 1 - htop
dx = 0.01
dy = 0.02

set label "(a)" at screen dx, 1 - dy
set label "(b)" at screen dx, hbot - dy

set multiplot
set size 1, htop
set origin 0, hbot

set tmargin 0.5
set bmargin 0

set xtics 100 offset 0, 0.2
set mxtics 2
set format x ""
unset xlabel

#set format y "10^{/*0.8 %T}"
set ytics 0.02
set mytics 2
set ylabel "{/Times-Italic w_{ij}}" offset 1.5, 0


set key right top Left reverse font "Times, 22" spacing 1.0

plot [:][-0.:0.04] \
    "../../data/lj/rho0.1/win_sig0.2.dat"    u ($2):($1==  1?$3:1/0) w l lt 1 lw 2 t "{/Times-Italic j} = 1", \
    "../../data/lj/rho0.1/win_sig0.2.dat"    u ($2):($1==150?$3:1/0) w l lt 2 lw 2 t "{/Times-Italic j} = 150", \
    "../../data/lj/rho0.1/win_sig0.2.dat"    u ($2):($1==300?$3:1/0) w l lt 3 lw 2 t "{/Times-Italic j} = 300", \
    "../../data/lj/rho0.1/win_sig0.2.dat"    u ($2):($1==450?$3:1/0) w l lt 4 lw 2 t "{/Times-Italic j} = 450", \
    -1 notitle

set size 1, hbot
set origin 0, 0

set bmargin 3
set format x "%g"
set xlabel "{/Times-Italic i}" offset 0, 0.5

plot [:][-0.025:0.085] \
    "../../data/lj/rho0.1/win_kc20.dat"      u ($2):($1==  1?$3:1/0) w l lt 1 lw 2 t "{/Times-Italic j} = 1", \
    "../../data/lj/rho0.1/win_kc20.dat"      u ($2):($1==150?$3:1/0) w l lt 2 lw 2 t "{/Times-Italic j} = 150", \
    "../../data/lj/rho0.1/win_kc20.dat"      u ($2):($1==300?$3:1/0) w l lt 3 lw 2 t "{/Times-Italic j} = 300", \
    "../../data/lj/rho0.1/win_kc20.dat"      u ($2):($1==450?$3:1/0) w l lt 4 lw 2 t "{/Times-Italic j} = 450", \
    -1 notitle


unset multiplot
unset output
set terminal pop
reset
