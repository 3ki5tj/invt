#!/usr/bin/env gnuplot



# Normalized error vs the width of the Gaussian updating scheme



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 5 font "Times, 26"
set output "errsigscan.eps"
set multiplot


hbot = 0.56
htop = 1 - hbot

# bottom panel
set size 1, hbot
set origin 0, 0

set xtics 2 offset 0, 0.3
set mxtics 2
#set xlabel "{/Symbol-Oblique s} (Gaussian) or {/Times-Italic n}/{/Times-Italic K}/{/Symbol \326}(2{/Symbol-Italic p}) (Bandpass)" offset 0, 0.5
set xlabel "Width, {/Symbol-Oblique s}" offset 0, 0.5
set xrange [0:10]

set logscale y
set format y "10^{/*0.8 %T}"
#set ytics 1e3
#set format y "%.0t{/Symbol \264}10^{/*0.8 %T}"
#set mytics 2
set ylabel "Normalized error, ({/Times-Italic T} + {/Times-Italic t}_{0}) {/Times-Italic E}"

a0 = 0.0001
fac = 2*100/sqrt(2*pi)

set title "Local" offset 0, -0.5

plot [:][5e3:1.4e4] \
    "../../data/scan/sigprd_t1e8_l.dat"  u 1:($6**2) w l lt 1 lw 4 notitle, \
    "../../data/scan/sigrun_t1e8_l.dat"  u 1:($4)    w p pt 13 ps 1.4      notitle, \
    "../../data/scan/sigprd_t1e10_l.dat" u 1:($6**2) w l lt 1 lw 2 notitle, \
    "../../data/scan/okrun_t1e8_l.dat"   u (fac/(2*$1)):($4)    w p  pt 5                             notitle, \
    "../../data/scan/okprd_t1e8_l.dat"   u (fac/(2*$1)):($6**2) w l       lt 4 lw 4  lc rgb "#000000" notitle, \
    "../../data/scan/okprd_t1e10_l.dat"  u (fac/(2*$1)):($6**2) w l       lt 2 lw 2  lc rgb "#a0a0a0" notitle, \
    -1 notitle



# top panel
set origin 0, hbot
set size 1, hbot

set tmargin 5
set bmargin 0

# `width` to reduce the text length
set key Left reverse width -7.5 spacing 1.0 font "Times, 20"

unset xlabel
set format x ""

set yrange [3:6e3]

set title "Global"

plot [:][:] \
    "../../data/scan/sigprd_t1e8_g.dat"  u 1:($6**2)            w l        lt 1 lw 4                   t "Gaussian, {/Times-Italic T} = 10^{/*0.8 8}", \
    "../../data/scan/sigrun_t1e8_g.dat"  u 1:($4)               w p  pt 13 ps 1.4                      notitle, \
    "../../data/scan/sigprd_t1e10_g.dat" u 1:($6**2)            w l        lt 1 lw 2                   t "Gaussian, {/Times-Italic T} = 10^{/*0.8 10}", \
    "../../data/scan/okrun_t1e8_g.dat"   u (fac/(2*$1)):($4)    w p  pt 5                             notitle, \
    "../../data/scan/okprd_t1e8_g.dat"   u (fac/(2*$1)):($6**2) w l        lt 4 lw 4  lc rgb "#000000" t "Gaussian, {/Times-Italic T} = 10^{/*0.8 8}", \
    "../../data/scan/okprd_t1e10_g.dat"  u (fac/(2*$1)):($6**2) w l        lt 2 lw 2  lc rgb "#a0a0a0" t "Bandpass, {/Times-Italic T} = 10^{/*0.8 10}", \
    -1 notitle

# phantom plot for the legend
set key at 1.4, 5e3
plot [:][:] \
  -1 w p pt 13 ps 1.4 t " ", \
  -1 w p pt 13 ps 1.4 lc rgb "#ffffff" t " ", \
  -1 w p pt  5 ps 1.0 t " "


unset multiplot
unset output
set terminal pop
reset
