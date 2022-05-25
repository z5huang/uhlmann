unset multiplot
reset
eval setdef('tex','0')

output='hof-erg'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 2.5in,2.5in  #fontscale 0.6
  set output output.'.tex'
}

fn='erg'
dat=fn.'.dat'
edge='erg_edge.dat'
par=fn.'.par'
load par

###
# do stuff
set xlabel '$k_x/\pi$' offset first 0,first 0.3
set xtics 0,2
set xrange [0:2]
set mxtics q
set ylabel '$E$'

#set ytics -4,2,4

p for [i=0:q-1] dat every q::i u ($1/pi):(-$4) w d lc -1 notit, \
  for [i=0:q-2] edge every q-1::i u ($1/pi):(-$4/$5) w l lt 1 notit, \
  for [i=0:q-2] edge every q-1::i u ($1/pi):(-$6/$7) w l lt 3 notit

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
