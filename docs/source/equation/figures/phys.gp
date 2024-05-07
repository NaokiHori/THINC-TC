reset

xmin = + 0.0
xmax = + 4.5
ymin = - 0.5
ymax = + 4.0

lx = xmax-xmin
ly = ymax-ymin

set terminal epslatex standalone color size lx,ly font ',12'
set output 'phys.tex'

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
rmin = 1.
rmax = rmin + 3.
tmin = 0.
tmax = pi / 3.
deltar = (rmax - rmin) / num
deltat = (tmax - tmin) / num

do for [i = 0 : num : 1] {
  if (0 == i) {
    r = rmin
  } else if (num == i) {
    r = rmax
  } else {
    r = i - num / 2
    r = rmin + (rmax - rmin) * 1. / (1. + exp(- 1. * r))
  }
  set object circle at first 0., first 0. size r arc [180. / pi * tmin : 180. / pi * tmax] fc rgb '#000000' lw 10 nowedge
}
do for [i = 0 : num : 1] {
  t = tmin + i * deltat
  set arrow from first rmin * cos(t), first rmin * sin(t) to first rmax * cos(t), first rmax * sin(t) nohead ls 1
}

plot \
  NaN notitle

