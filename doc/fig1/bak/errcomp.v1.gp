#!/usr/bin/env gnuplot



# component of the error function



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4.0 make dashed line longer
set terminal postscript eps enhanced dl 4.0 size 5, 3.5 font "Times, 32"
set output "errcomp.eps"


set xtics 0.5 offset 0, 0.3
set mxtics 5
set xlabel "{/Times {/Symbol-Oblique n}_{/Times-Italic k} = {/Symbolic-Oblique l}_{/Times-Italic k} / {/Symbolic-Oblique l} }" offset 0, 0.5

#set logscale y
#set format y "10^{%T}"
set ytics 2
set mytics 2
set ylabel "{/Times {/Times-Italic G}({/Symbol-Oblique n}_{/Times-Italic k}, {/Times-Italic r})}" offset 1, 0

# `width` to reduce the text length
set key Left reverse width -2 spacing 1.8

set samples 1000
G(x, y) = x * ( 1 + (1 - y**(2*x-1))/(2*x-1) )

plot [0:2][0.05:7] \
    G(x, 1e-1)  lt 1 lw 5 t "{/Times {/Times-Italic r} = 0.1}", \
    G(x, 1e-2)  lt 2 lw 5 t "{/Times {/Times-Italic r} = 0.01}", \
    G(x, 1e-3)  lt 4 lw 5 t "{/Times {/Times-Italic r} = 0.001}", \
    -1 notitle


unset output
set terminal pop
reset
