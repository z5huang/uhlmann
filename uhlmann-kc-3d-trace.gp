unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-haldane-kc-3d-trace'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle solid size 3in,2.5in  #fontscale 0.6
  set output output.'.tex'
}

udat = 'uhlmann_kc_mu.dat'
wdat = 'wannier.dat'

###
# do stuff

set xlabel '$T$' 
set ylabel '$\mu$'
#set zlabel '$w:$' offset first 0.1, first 0, first 0.6
set xtics 0,0.5,2 offset 0,-0.2,0
set ytics -0.5,0.5,1.5 offset first -0.05,0,0
set mytics 5
set ztics -1,0.5,1
set ztics add ('$\textsf{Tr}(\rho_0U_{0N})$' 0)
set xrange [0:1]
set yrange [-0.2:1.1]
set zrange [-1:1]

set hidden3d
set xyplane 0.1
set view 60,330
sp udat u 2:1:($7*cos($8)) w l lt 2 notit, 0 notit w l lt 0
   

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
