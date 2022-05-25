unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-hof-amp'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 2.5in,2.5in  #fontscale 0.6
  set output output.'.tex'
}

fn='uhlmann'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff
#set multiplot
#set origin 0,0
#set size 1,1

set xlabel '$k_x/\pi$'
set ylabel '$|m|$'
set ytics

set xrange [0:2]
set yrange [-0.05:]
set ytics 0,0.1

lev_from=1
lev_to=7
lev_split=4

### plot the rescaled magnitude
set key horizontal opaque center right spacing 1.5 maxrows 2 samplen 3
p for [i=lev_from:lev_to] dat u ($2/pi):((column(2*i+3))) w l lc i lw 2 tit (i == 1 ? sprintf('$a: %d$',i) : sprintf('%d',i))

#set logscale y
#set ytics 0, 100
#set mytics 10
#p for [i=lev_from:lev_to] dat u ($2/pi):((column(2*i+3))) w l notit

#unset xlabel
#unset ylabel
#set xtics 1
##set ytics 0.1
###set yrange [-0.1:]
#set origin 0.05,0.2
#set size 0.4,0.4
#p for [i=lev_from:lev_split] dat u ($2/pi):((column(2*i+3))) w l lc i notit
#
#set origin 0.5,0.2
#set size 0.4,0.4
#p for [i=lev_split+1:lev_to] dat u ($2/pi):((column(2*i+3))) w l lc i notit

#p for [i=lev_from:lev_to] dat u ($2/pi):(column(2*i+3)) w l notit
#
#unset multiplot

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
