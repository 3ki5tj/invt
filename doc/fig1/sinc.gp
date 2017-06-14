#!/usr/bin/env gnuplot



# window kernel function and updating matrix
# for the bandpass updating schemes


set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3 make dashed line longer
set terminal postscript eps enhanced font "Times, 24"
set output "sinc.eps"
set multiplot


dx = 0.01
dy = 0.05

ssr(l, K, N) = (l == 0 ? (2*K+1.)/N : sin((2*K+1)*l*pi/N)/N/sin(l*pi/N));

wij_np(i, j, K, n) = ssr(i - j, K, 2*n) + ssr(i + j - 1, K, 2*n)
wij_pr(i, j, K, n) = ssr(i - j, K, n)

# bottom panel
set size 1, 1
set origin 0, 0

set tmargin 0.5
set lmargin 7

set xtics offset 0, 0.2
set mxtics 10
set xlabel "{/Times-Italic i}" offset 0, 0.5

set ytics 0.05 offset 0.2, 0
set mytics 5
set ylabel "{/Times-Italic w}_{/Times-Italic ij}" offset 2, 0

set key right top Left reverse invert spacing 1.2 width -8 font "Times, 24"

n = 100
plot [1:n][-0.03:0.17] \
    0 lt 1 lw 0.1 lc rgb "#cccccc" notitle, \
    wij_np(x, 4,   10, n) lt 1 lw 1.0 t "Bandpass, non{/*0.7 -}periodic, {/Times-Italic K} = 10, {/Times-Italic j} = 4", \
    wij_pr(x, 4,   5,  n) lt 2 lw 1.0 t "Bandpass, periodic, {/Times-Italic K} = 5, {/Times-Italic j} = 4", \
    wij_np(x, n/2, 10, n) lt 1 lw 2.0 t "Bandpass, non{/*0.7 -}periodic, {/Times-Italic K} = 10, {/Times-Italic j} = 50", \
    wij_pr(x, n/2, 5,  n) lt 2 lw 2.0 t "Bandpass, periodic, {/Times-Italic K} = 5, {/Times-Italic j} = 50"

unset xlabel
unset ylabel

unset multiplot
unset output
set terminal pop
reset
