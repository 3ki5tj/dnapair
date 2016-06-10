#!/usr/bin/env gnuplot

# plot the two-dimensional PMF

set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3 make dashed line longer
set terminal postscript eps enhanced dl 3 size 7, 6 font "Times, 24"
set output "pmf.eps"
set multiplot

set rmargin 0

set origin 0, 0
set size 0.65, 1

# draw tics along x every 1.0 A
set xtics 1.0 offset -1, -0.5
set mxtics 2

set ytics ("0" 0, \
           "" pi*1/3, \
           "" pi*2/3, \
           "{/Symbol p}" pi, \
           "" pi*4/3, \
           "" pi*5/3, \
           "2{/Symbol p}" 2*pi)

set ztics 5
set mztics 5

# set view point
set view 70, 75

set xyplane 0.2

splot [:26.5][-0.01:2*pi] "mfrt.out" w lp pt 7 ps 1.5 notitle



set origin 0.5, 0
set size 0.5, 1

# set view point
set view 80, 25

replot

unset multiplot
unset output
set terminal pop
reset

