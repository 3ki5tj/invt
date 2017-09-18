#!/usr/bin/env gnuplot



# optimized schedule for the Lennard-Jones system



set encoding cp1250 # make the minus sign longer
set terminal push
# dl 2.0 make dashed lines longer
set terminal postscript eps enhanced font "Helvetica, 32"
set output "lj_alpha.eps"


reset


htop = 0.5
hbot = 1 - htop

dx = 0.01
dy = 0.05

set lmargin 8

set logscale x
#set xtics 0.2 offset 0, 0.2
#set mxtics 2
set format x "10^{/*0.8 %T}"
set mxtics 10
set xrange [1e4:1e8]
set xlabel "Shifted simulation time, {/Times-Italic t} + {/Times-Italic t}_{/*0.8 0}"

set logscale y
set format y "10^{/*0.8 %T}"
#set ytics add ("{/Times {/Symbol a}_0/2}" alpha0*0.5)
set mytics 10
set yrange [6e-9:1e-4]
set ylabel "{/Symbol-Oblique a}({/Times-Italic t})"


set key left bottom Left reverse spacing 1.5

alpha0 = 1e-4
t0 = 2 / alpha0

set arrow from 6.0e3, alpha0*0.5 to 2e4, alpha0*0.5 lt 1 lw 1 filled size screen 0.02,10,35
set label "{/Times {/Times-Italic a}_0&{/*0.5 i}/2}" at 2.0e3, alpha0*0.5

fn = "../../data/lj/rho0.1/alpha_sig0.29_t1e8.dat"
k = -2.0
b = -log(alpha0*0.5)
f(x) = k*1e-7*x - b
fit [0:5e6] f(x) fn u 1:(log($2)) via k, b
a = exp(-b)*1e5
print a, k, b

plot [t0:][:] \
    fn  u ($1+t0):($2) w l lt 1 lw 2 t "Optimal", \
    1/x lt 3 lw 2 t "1 / ({/Times-Italic t} + {/Times-Italic t}_{/*0.8 0})", \
    exp(f(x - t0)) lt 2 lw 5 t sprintf("{/Times %.1f {/Symbol \264} 10^{/*0.7-5}{/Times-Italic e}^{%.1f {/Symbol \264} 10^{/*0.7-7}{/Times-Italic t}}}", a, k), \
    -1 notitle

unset output
set terminal pop
reset
