#!/usr/bin/env gnuplot



# error vs the proportionality constant c
# of the inverse-time schedule of the updating magnitude
#   alpha(t) = c / (t + t0)
# for the single-bin (Wang-Landau) scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 32"
set output "errsbin.eps"


set xtics 1 offset 0, 0.3
set mxtics 2
set xlabel "{/Times 1 / }{/Symbol-Oblique l}" offset 0, 0.5

set logscale y
set format y "10^{%T}"
set mytics 10
set ylabel "Error, {/Times-Italic E}"

# `width` to reduce the text length
set key Left reverse width -5 spacing 1.0

# error bars are too small, no need to have them
# to add them back use
#   "u 1:2:4 w e w lp 1"
# instead of
#   "u 1:2 w p"

plot [0:5][5e-7:0.006] \
    "../../data/singlebin/singlebin_g_err.dat" u 1:2 every 2 w p pt 5 ps 2.5  t "Global, simulation", \
    "../../data/singlebin/singlebin_g_prd.dat" u 1:($2**2)   w l lt 1 lw 2    t "Global, theory", \
    "../../data/singlebin/singlebin_l_err.dat" u 1:2 every 2 w p pt 7 ps 2.5  t "Local, simulation", \
    "../../data/singlebin/singlebin_l_prd.dat" u 1:($2**2)   w l lt 2 lw 2    t "Local, theory"

unset output
set terminal pop
reset
