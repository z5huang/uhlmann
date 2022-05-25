unset multiplot
reset
eval setdef('tex','0')

output='wannier-bhz'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 3in,3in  #fontscale 0.6
  set output output.'.tex'
}

fn='wannier'
dat=fn.'.dat'
#par=fn.'.par'
#load par

###
# do stuff

set xlabel '$k_x/\pi$'
set ylabel '$\gamma/\pi$'
set xrange [0:2]
set yrange [-1:1]
set xtics 0,1,2
set ytics -1,1,1
p dat u ($2/pi):($7/pi) w p pt 6 ps 0.5 lc 1 notit, \
  dat u ($2/pi):($9/pi) w p pt 6 ps 0.5 lc 1 notit

#
###

if (tex == 1){
  unset output
  set term wxt
  build = buildtex(output)
  print '"eval build" to build and preview'
} else {
  #print "press enter to continue"
  #pause -1
}
