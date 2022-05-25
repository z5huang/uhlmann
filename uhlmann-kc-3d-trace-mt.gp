unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-bhz-kc-3d-trace'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle solid size 3in,2.5in  #fontscale 0.6
  set output output.'.tex'
}

udat = 'uhlmann_kc_m.dat'
wdat = 'wannier.dat'

###
# do stuff

set ylabel '$T$' 
set xlabel '$m$'
#set zlabel '$w:$' offset first 0.1, first 0, first 0.6
set ytics 0,0.1,1 offset 0,-0.2,0
set xtics -1,2,5 offset first -0.05,0,0
set mxtics 2
set ztics -1,0.5,1
set ztics add ('$\textsf{Tr}(\rho_0U_{0N})$' 0)
set yrange [0:0.5]
set xrange [-1:5]
set zrange [-1:1]

set hidden3d
set xyplane 0.1
set view 42,21
sp udat u 1:2:($7*cos($8)) w l lt 3 notit, 0 notit w l lt 0
   

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
