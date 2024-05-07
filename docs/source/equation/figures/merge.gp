reset

xmin = 0.0
xmax = 9.0
ymin = 0.0
ymax = 5.25

lx = xmax - xmin
ly = ymax - ymin

set terminal epslatex standalone color size lx,ly font ',20'
set output 'result.tex'

unset border

set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

unset xlabel
unset ylabel

set xrange [xmin:xmax]
set yrange [ymin:ymax]

unset xtics
unset ytics

set size ratio -1

ref = 4.5

set label 'Physical      coordinate' center at first 0.5 * ref, first 1. * ref + 0.35 textcolor rgb '#000000'
set label 'Computational coordinate' center at first 1.5 * ref, first 1. * ref + 0.35 textcolor rgb '#000000'
set label '$\left( r, \theta \right)$'         center at first 0.5 * ref, first 1. * ref - 0.05 textcolor rgb '#000000'
set label '$\left( \xi^r, \xi^\theta \right)$' center at first 1.5 * ref, first 1. * ref - 0.05 textcolor rgb '#000000'

set label '\includegraphics[width=4.500in, height=4.500in]{phys.pdf}' center at first 0.5 * ref, first 0.5 * ref
set label '\includegraphics[width=4.500in, height=4.500in]{comp.pdf}' center at first 1.5 * ref, first 0.5 * ref

plot \
  NaN notitle
