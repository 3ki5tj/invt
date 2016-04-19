#!/usr/bin/env gnuplot



# optimized schedule versus the optimal inverse-time schedule
# for Gaussian updating scheme with sigma = 10



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed line longer
set terminal postscript eps enhanced dl 2.0 size 3.5, 5.0 font "Times, 24"
set output "optacmp.eps"
set multiplot


reset


htop = 0.5
hbot = 1 - htop

dx = 0.01
dy = 0.05

set lmargin 7


set label "(a)" at screen dx, 1 - dy    font "Times, 32"
set label "(b)" at screen dx, hbot - dy font "Times, 32"

# bottom panel
set size 1, hbot
set origin 0, 0

set xtics 5 offset 0, 0.2
set mxtics 5
set xlabel "{/Times-Italic q}({/Times-Italic T}) - {/Times-Italic q}" offset 0, 0.3

set logscale y
set format y "10^{%T}" 
set ylabel "{/Times-Italic m}({/Times-Italic q}) = 1/({/Times-Italic T}{/Symbol-Oblique a})"

set key right top Left reverse spacing 1.5 width -4

plot [0:25][1e-3:] \
    "../../data/opta/sig10_l_alpha.dat"         u (162.972-$3):(1/(1e8*$2)) w l lt 1 lw 2  t "Gaussian, {/Times-Italic T} = 10^8", \
    "../../data/opta/sig10_l_t5e7_alpha_q.dat"  u (162.972-$3):(1/(5e7*$2)) every 30 w p pt 6 lw 2  t "Gaussian, {/Times-Italic T} = 5{/Symbol \264}10^7", \
    "../../data/opta/wl_l_alpha.dat"            u (8.51733-$3):(1/(1e8*$2)) w l lt 2 lw 2  t "Single{/*0.7 -}bin, {/Times-Italic T} = 10^8", \
    0 notitle




# top panel
set origin 0, hbot
set size 1, htop

#set logscale x
set xtics 0.2 offset 0, 0.2
set mxtics 2
#set format x "10^{%T}"
#set xrange [1e-4:1]
set xlabel "Simulation time, {/Times-Italic t} / {/Times-Italic T}" offset 0, 0.0

set logscale y
set format y "10^{%T}"
set mytics 10
#set yrange [1e4:1e9]
set ylabel "{/Times-Italic T} {/Symbol-Oblique a} ({/Times-Italic t})"


#set key left bottom Left reverse spacing 1.5
set key Left reverse spacing 1.5 width -5

plot [:][1:1e4] \
    "../../data/opta/sig10_l_alpha.dat"      u ($1/1e8):(1e8*$2) w l lt 1 lw 2 lc rgb "#000000"    t "Optimal, {/Times-Italic T} = 10^8", \
    "../../data/opta/sig10_l_t5e7_alpha.dat" u ($1/5e7):(5e7*$2) every 200 w p pt 6 lw 2 lc rgb "#000000"    t "Optimal, {/Times-Italic T} = 5{/Symbol \264}10^7", \
    "../../data/opta/sig10_l_alpha.dat"      u ($1/1e8):(1e8*$5) w l lt 2 lw 2 lc rgb "#000000"    t "Inverse{/*0.7 -}time, {/Times-Italic T} = 10^8", \
    -1 notitle




unset multiplot
unset output
set terminal pop
reset
