unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-haldane'
if (tex == 1) {
  #set lmargin at screen 0
  #set rmargin at screen 1
  #set bmargin at screen 0
  #set tmargin at screen 1
  set term lua tikz standalone createstyle size 5in,3in  #fontscale 0.6
  set output output.'.tex'
}

udat = 'uhlmann.dat'
wdat = 'wannier.dat'

###
# do stuff
set multiplot
set origin 0,0
set size 1,1

### plot the phases
set xlabel '$k_x/\pi$'
set ylabel '$\gamma/\pi$'
set xrange [0:2]
set yrange [-1:1]

set key opaque
set key bottom right spacing 2
p wdat u ($2/pi):($4/pi - 2*$5) w l lt -1 notit, \
  '' u ($2/pi):($4/pi + 2 - 2*$5) w l lt -1 tit '$\gamma_b$', \
  udat u ($2/pi):($6/pi) tit '$\gamma_u^{>}$' w p lc 1 pt 4, \
  udat u ($2/pi):($4/pi) tit '$\gamma_u^{<}$' w p lc 3 pt 9, 
  

### plot the magnitudes
unset xlabel
unset ylabel
set origin 0.05, 0.7
set size 0.4, 0.4

clear # white out the region

set key center left samplen 1
set xtics 0,1
set mxtics 5
unset ytics
set y2tics 0,1 mirror
set mytics 5
set y2range [-0.1:1.1]

#p udat u ($2/pi):5 w p lc 1 pt 4 tit '$w_u^{>}$', \
#'' u ($2/pi):3 w p lc 3 pt 9 tit '$w_u^{<}$'

p udat u ($2/pi):5 axes x1y2 w l lc 1 tit '$w_u^{>}$', \
'' u ($2/pi):3 axes x1y2 w l lc 3 tit '$w_u^{<}$' 
unset multiplot
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
