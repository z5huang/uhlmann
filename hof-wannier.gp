unset multiplot
reset
eval setdef('tex','0')

output='hof-wannier'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 2.5in,2.5in  #fontscale 0.6
  set output output.'.tex'
}

fn='wannier'
dat=fn.'.dat'
#edge='erg_edge.dat'
#par=fn.'.par'
#load par

###
# do stuff
q=7
set xlabel '$k_x/\pi$' offset first 0,first 0.1
set xtics 0,2
set xrange [0:2]
set mxtics q
set ylabel '$\gamma/\pi$'
set yrange [-1:1]
set ytics -1,1

p for [i=1:3] dat every 4 u ($2/pi):(column(i*2 + 5)/pi) w p lt 1 pt 6 ps 0.5 notit, \
'' u ($2/pi):($4/pi - 2*$5) w l lc 3 lt 1 lw 1 notit, \
  '' u ($2/pi):($4/pi - 2*$5 - 2) w l lc 3 lt 1 lw 1 notit

#
###

if (tex == 1){
  unset output
  set term qt
  build = buildtex(output)
  print '"eval build" to build and preview'
} else {
  #print "press enter to continue"
  #pause -1
}
