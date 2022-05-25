unset multiplot
reset
eval setdef('tex','0')

output='uhlmann-hof-arg'
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
set ylabel '$\textsf{arg}/\pi$' offset first 0.1
set ytics 0.5

set xrange [0:2]
set yrange [-1:1]
#set ytics 0,0.1

lev_from=1
lev_to=7
lev_split=4

tot(n) = 18+n

### plot the rescaled magnitude
set key opaque at first 2, first 0.8 spacing 1.5
p dat u ($2/pi):(column(tot(3))/pi) w l lt 1 lw 2 tit '$\gamma_{1..3}$', \
  '' u ($2/pi):((column(tot(6)) - column(tot(3)))/pi) w l lt 3 lw 2 tit '$\gamma_{4..6}$', \
  '' u ($2/pi):((column(tot(7)) - column(tot(6)))/pi) w l lt 5 lw 2 tit '$\gamma_{7}$'

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
