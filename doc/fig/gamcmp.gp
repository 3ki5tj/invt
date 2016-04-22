#!/usr/bin/env gnuplot



# comparison of the Gamma values



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 4.0 make dashed line longer
set terminal postscript eps enhanced dl 4.0 size 5, 3.5 font "Times, 32"
set output "gamcmp.eps"

set logscale x
set format x "10^{%T}"
set xtics 100 offset 0, 0.2
set xlabel "{/Symbol G}^{&{i}one{/*0.7 -}step}"

set logscale y
set format y "10^{%T}"
set ylabel "{/Symbol G}^{&{i}MD}"

#set key Left reverse width -2 spacing 1.8

fngamma = "../../data/gamcmp/gamma.dat"

#f(x) = a * x + b
#fit f(x) fngamma u 2:4 via a, b

plot [1e-4:1e4][:] fngamma u 2:4 w p pt 7 ps 0.8 notitle


unset output
set terminal pop
reset
