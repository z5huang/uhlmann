unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-bhz-ReTr'
if (tex == 1) {
  #set lmargin at screen 0
  #set rmargin at screen 1
  #set bmargin at screen 0
  #set tmargin at screen 1
  set term lua tikz standalone createstyle size 4in,3in  #fontscale 0.6
  set output output.'.tex'
}

udat = 'uhlmann.dat'
wdat = 'wannier.dat'

###
# do stuff
set xlabel '$k_x/\pi$'
set ylabel '$\textsf{Tr}(\rho_0U_{0N})$' offset first 0.1
set xrange [0:2]
set yrange [-1:1]

set key opaque
set key bottom right spacing 2
p udat u ($2/pi):($7*cos($8)) w l lc 1 notit, \
  0 w l lt 0 notit

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
