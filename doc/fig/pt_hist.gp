#!/usr/bin/env gnuplot



set encoding cp1250 # make the minus sign longer
set terminal push
set terminal postscript eps enhanced size 5, 5.7 font 28
set output "pt_hist.eps"
set multiplot


reset

ht1 = 0.23
ht2 = 0.23
ht3 = 0.22
ht4 = 1 - ht1 - ht2 - ht3

set label 1 at screen 0.00, 0.98         "(a)"
set label 2 at screen 0.00, ht2+ht3+ht4-0.02 "(b)"
set label 3 at screen 0.00, ht3+ht4-0.03     "(c)"
set label 4 at screen 0.00, ht4-0.03     "(d)"


set style line 1 lt 1
set style line 2 lt 2 lw 3
set style line 3 lt 3 lw 1

set origin 0, ht2 + ht3 + ht4
set size 1, ht1

set lmargin 8
set bmargin 0
set tmargin 0.2
set rmargin 0.5

set format x ""
set xtics 1000
set mxtics 5
set xrange [-7600:-3100]
#set xlabel "Shifted simulation time, {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0}"
#
#set format y "10^{/*0.8 %T}"
set ytics 0.002 offset 0.5, 0
set mytics 2
set yrange [0:0.007]
set ylabel "Histogram" offset 2.5, -0.5 font "Helvetica, 28"


#set key left bottom Left reverse font "Times, 22" spacing 1.5

plot [:][:] \
    "../../data/pt/pt2gaus_L64b.his" u 1:($5 >= 0 ? $2 : 1/0) w l ls 1 t "", \
    -1 notitle

# "../../data/pt/pt2gaus_L64b.his" u 1:($5 < 0  ? ($2*20) : 1/0) w l ls 2 t "", \


set origin 0, ht3+ht4
set size 1, ht2

set ytics 5
set mytics 5
set yrange [-19:-6]
set ylabel '{/Times-Italic ~c{0.2\^}}@_{/Times 0}^{/Times ({/Times-Italic s})} - {/Symbol-Oblique b}_{/Times-Italic c }{/Times-Italic E@_{c}^{/Times  ({/Times-Italic s})}}' offset 0.0, -0.5

bc = 1.42529047
shiftlnz = 7599.78

plot [:][:] \
    "../../data/pt/pt2gaus_L64b.dat" u 2:($6-$5/sqrt(2)-bc*$2-shiftlnz) w p pt 6 lc "black" notitle, \
    "../../data/pt/lng_L64b.dat"     u 1:($2-bc*$1)      w l ls 2 t ""


#    "../../data/pt/lng_L64b.dat"     u 1:($1 < -3200 ? $3 : 1/0)                 w l ls 2 lw 1 t "", \


set origin 0, ht4
set size 1, ht3

set ytics 0.01
set mytics 2
set yrange [1.405:1.445]

set ylabel "{/Times-Italic c}@_{/Times 1}^{/Times ({/Times-Italic s})}/{/Symbol-Oblique s}_{/Times-Italic E}" offset 2.0, 0


plot [:][:] \
    "../../data/pt/pt2gaus_L64b.dat" u 2:($4/64) w p pt 6 lc "black" notitle, \
    bc w l ls 2 t ""


#    "../../data/pt/lng_L64b.dat"     u 1:($1 < -3200 ? $3 : 1/0)                 w l ls 2 lw 1 t "", \


set origin 0, 0
set size 1, ht4

set bmargin 2.8

set xtics 1000 offset 0, 0.3
set mxtics 5
set format x "%g"
set xlabel "Energy, {/Times-Italic E}" offset 0, 0.3

set ytics 0.5
set mytics 5
set yrange [0.3:1.2]

set ylabel "{/Symbol \326}{/Times ~2{.7-}} {/Times-Italic c}@_{/Times 2}^{/Times ({/Times-Italic s})}" offset 1.9, 0


plot [:][:] \
    1 w l ls 2 notitle, \
    "../../data/pt/pt2gaus_L64b.dat" u 2:($5*sqrt(2)) w p pt 6 lc "black" notitle



unset multiplot

unset output
set terminal pop
reset
