#!/usr/bin/env gnuplot



# Error components for the Gaussian updating scheme
# On the Lennard-Jones system



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced font "Helvetica, 32"
set output "lj_xerr.eps"


set size 1, 1
set origin 0, 0

set rmargin 0.5

set xtics 5 offset 0, 0.3
set mxtics 5
set xlabel "Eigenmode, {/Times-Italic k}" offset 0, 0.5
set xrange [0:31]

set logscale y
set format y "10^{/*0.8 %T}"
set ytics 100 offset 0.5, 0
set ylabel "Error component" offset 1, 0

set key Left reverse samplen 3.0 width -6 font "Helvetica, 26"

plot [:][1e-9:1e-2] \
    "../../data/lj/rho0.1/xerr_sig0.28_t1e8_opt.dat"       u 1:($4+$5**2) w lp lt 1 lw 1 pt 12 ps 1.5 t "Gaussian, initial", \
    "../../data/lj/rho0.1/xerr_sig0.28_t1e8_opt.dat"       u 1:($2+$3**2) w lp lt 2 lw 1 pt  7 ps 1.0 t "Gaussian, optimal, final", \
    "../../data/lj/rho0.1/xerr_sig0.28_t1e8_invt.dat"      u 1:($2+$3**2) w lp lt 3 lw 1 pt  6 ps 1.5 t "Gaussian, {/Times 1/}{/Times-Italic t}, final", \
    "../../data/lj/rho0.1/xerr_sig0.28_t1e8_invt.dat"      u 1:($6+$7**2) w lp lt 4 lw 2 pt  2 ps 1.5 t "Gaussian, {/Times 1/}{/Times-Italic t}, corrected", \
    "../../data/lj/rho0.1/xerr_sbin_t1e8.dat"              u 1:($2+$3**2) w lp lt 5 lw 1 pt  4 ps 1.5 t "Single-bin, final", \
    -1 notitle


#    "../../data/lj/rho0.1/xerr_kc20_t1e8.dat"       u 1:($4+$5**2) w lp t "Bandpass, 1/{/Times-Italic t}, initial",
#    "../../data/lj/rho0.1/xerr_kc20_t1e8.dat"       u 1:($2+$3**2) w lp t "Bandpass, 1/{/Times-Italic t}, final",


unset output
set terminal pop
reset
