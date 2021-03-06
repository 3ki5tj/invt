#!/usr/bin/env gnuplot



# error vs the proportionality constant c
# of the inverse-time schedule of the updating magnitude
#   alpha(t) = c / (t + t0)
# for the single-bin (Wang-Landau) updating scheme and
# the triple-bin updating scheme


set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3.0 make dashed line longer
set terminal postscript eps enhanced dl 3.0 size 5, 3.5 font "Times, 32"
set output "errinvt.eps"
set multiplot

set logscale x
#set xtics 1 offset 0, 0.3
#set mxtics 2
set xrange [0.05:10]
set xlabel "{/Symbol-Oblique l}" offset 0, 0.5

set logscale y
set format y "10^{/*0.8 %T}"
set mytics 10
set yrange [5e-7:0.1]
set ylabel "Error, {/Times-Italic E}"

# `width` to reduce the text length
set key left Left reverse spacing 1.2

# error bars are too small, no need to have them
# to add them back use
#   "u 1:2:4 w e w lp 1"
# instead of
#   "u 1:2 w p"

set arrow 1 from 1, 5e-7 to 1, 5e-4 lt 2 lc rgb "#808080" nohead back

plot [:][:] \
    "../../data/invt/wl_g_err.dat" u (1/$1):2         w p  pt 7  ps 2.0  notitle, \
    "../../data/invt/wl_g_prd.dat" u (1/$1):($2**2)   w l  lt 1  lw 2    t "Single{/*0.7 -}bin, global", \
    "../../data/invt/wl_l_err.dat" u (1/$1):2         w p  pt 6  ps 2.0  notitle, \
    "../../data/invt/wl_l_prd.dat" u (1/$1):($2**2)   w l  lt 1  lw 1    t "Single{/*0.7 -}bin, local"


# phantom plot for the legend
set key at 0.035, 0.083 Left reverse
plot [:][:] \
  -1 w p pt 7 ps 2.0 t " ", \
  -1 w p pt 6 ps 2.0 t " ", \
  -1 notitle

unset multiplot
unset output
set terminal pop
reset
