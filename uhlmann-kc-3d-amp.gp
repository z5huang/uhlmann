unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-haldane-kc-3d-amp'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle solid size 3in,2.5in  #fontscale 0.6
  set output output.'.tex'
}

udat = 'uhlmann_kc.dat'
wdat = 'wannier.dat'

###
# do stuff

set xlabel '$T$' 
set ylabel '$\kappa/\pi$'
#set zlabel '$w:$' offset first 0.1, first 0, first 0.6
set xtics 0,0.5,2 offset 0,-0.2,0
set ytics 0,0.5,2 offset 0.5,0,0
set ztics 0,0.5,1
set ztics add ('$w$' 0.5)
set xrange [0:1]
set yrange [0:2]

set hidden3d
set xyplane 0.1
sp udat every :4 u 3:($2/pi):($4-0.0001) w l lt 4 notit, \
   udat every :4 u 3:($2/pi):6 w l lt 2 notit
   

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
