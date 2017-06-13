#!/usr/bin/env gnuplot



# Error components for the Gaussian updating scheme
# On the Lennard-Jones system



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 font "Times, 28"
set output "lj_xerr.eps"


set size 1, 1
set origin 0, 0

set xtics 10 offset 0, 0.3
set mxtics 10
set xlabel "Eigenmode, {/Times-Italic k}" offset 0, 0.7
set xrange [0:35]

set logscale y
set format y "10^{/*0.8 %T}"
set ytics 100 offset 0.5, 0
set ylabel "Error component" offset 1, 0

set key Left reverse samplen 3.0 width -8

plot [:][1e-8:3e-1] \
    "../../data/lj/rho0.1/xerr_sig0.2_opt.dat"       u 1:($4+$5**2) w lp lt 1 lw 2 pt  2 ps 1.4 t "Gaussian, initial", \
    "../../data/lj/rho0.1/xerr_sig0.2_invt.dat"      u 1:($2+$3**2) w lp lt 1 lw 1 pt  6 ps 1.4 t "Gaussian, 1/{/Times-Italic t}, final", \
    "../../data/lj/rho0.1/xerr_sig0.2_opt.dat"       u 1:($2+$3**2) w lp lt 1 lw 1 pt  7 ps 1.4 t "Gaussian, opt., final", \
    "../../data/lj/rho0.1/xerr_sig0.2_invt.dat"      u 1:($6+$7**2) w lp lt 1 lw 1 pt 12 ps 1.4 t "Gaussian, 1/{/Times-Italic t}, corrected", \
    "../../data/lj/rho0.1/xerr_sbin.dat"             u 1:($4+$5**2) w lp lt 1 lw 2 pt  1 ps 1.8 t "Single-bin, initial", \
    "../../data/lj/rho0.1/xerr_sbin.dat"             u 1:($2+$3**2) w lp lt 1 lw 1 pt  4 ps 1.4 t "Single-bin, final", \
    -1 notitle


unset output
set terminal pop
reset
