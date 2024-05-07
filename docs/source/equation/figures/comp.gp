reset

xmin = - 0.5
xmax = + 4.0
ymin = - 0.5
ymax = + 4.0

lx = xmax-xmin
ly = ymax-ymin

set terminal epslatex standalone color size lx,ly font ',12'
set output 'comp.tex'

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
set samples 10000

set style line 1 lc rgb '#000000' lw 10

num = 8
xmin = 0.
xmax = 3.5
ymin = 0.
ymax = 3.5
deltax = (xmax - xmin) / num
deltay = (ymax - ymin) / num

do for [i = 0 : num : 1] {
  x = xmin + i * deltax
  set arrow from first x, first ymin to first x, first ymax nohead ls 1
  y = ymin + i * deltay
  set arrow from first xmin, first y to first xmax, first y nohead ls 1
}

plot \
  NaN notitle

