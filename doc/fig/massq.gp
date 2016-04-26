#!/usr/bin/env gnuplot



# component of the error function



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4.0 make dashed line longer
set terminal postscript eps enhanced dl 4.0 size 3.5, 2.7 font "Times, 32"
set output "massq.eps"

xmax = 2.5
qT = 2.2

set xtics ("0" 0, "{/Times-Italic q}({/Times-Italic T})" qT) offset 0, 0.2
set xlabel "{/Times-Italic Q}" offset 0, 0.6
#set xlabel "{/Times-Italic q}({/Times-Italic T}) - {/Times-Italic q}" offset 0, 0.6

unset ytics
set ylabel "{/Times-Italic m}({/Times-Italic Q})"

# `width` to reduce the text length
set key Left reverse width -2 spacing 1.8

x1 = 0.7
f(x) = exp(-0.5*x*x)
x2 = x1 - 0.25
set arrow 1 from x2, 0.4*f(x1) to x2, 0     head lw 2
set arrow 2 from x2, 0.6*f(x1) to x2, f(x1) head lw 2

x3 = x2*2 - x1

set arrow 10 from x3, f(x1) to x1, f(x1) nohead lw 2

set label 20 "1/({/Times-Italic T}{/Symbol-Oblique a})" at x2 - 0.4, 0.5*f(x1) front


x4 = x1 + 0.4

set label 30 "{/Times-Italic t} / {/Times-Italic T}" at x4, 0.4*f(x4) front

set samples 2000

plot [0:xmax][0:1.1] \
    (x >= x1 && x <= qT ? f(x) : 1/0) w filledcu x1 lc rgb "#cccccc" notitle, \
    f(x)                   lt 4 lw 4 lc rgb "#888888" notitle, \
    (x <= qT ? f(x) : 1/0) lt 1 lw 4 notitle, \
    -1 notitle

# make the xtics on top of the plot
set grid
unset grid

unset output
set terminal pop
reset
