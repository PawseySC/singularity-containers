# set terminal pngcairo  background "#ffffff" enhanced font "arial,8" fontscale 1.0 size 540, 384 
# set output 'hidden2.3.png'
unset border
set dummy u, v
set format cb "%.1f" 
unset key
set style increment default
set parametric
set view 60, 30, 1.5, 0.9
set isosamples 50, 20
set size ratio 0 0.55,0.9
set origin 0.2,0
unset xtics
unset ytics
unset ztics
set title "PM3D surface\ndepth sorting" 
set urange [ -3.14159 : 3.14159 ] noreverse nowriteback
set vrange [ -3.14159 : 3.14159 ] noreverse nowriteback
set xrange [ * : * ] noreverse writeback
set x2range [ * : * ] noreverse writeback
set yrange [ * : * ] noreverse writeback
set y2range [ * : * ] noreverse writeback
set zrange [ * : * ] noreverse writeback
set cbrange [ * : * ] noreverse writeback
set rrange [ * : * ] noreverse writeback
set pm3d depthorder 
set colorbox user
set colorbox vertical origin screen 0.9, 0.15 size screen 0.02, 0.5 front  noinvert bdefault
f(x,y) = sin(-sqrt((x+5)**2+(y-7)**2)*0.5)
## Last datafile plotted: "blutux.rgb"
splot cos(u)+.5*cos(u)*cos(v),sin(u)+.5*sin(u)*cos(v),.5*sin(v) with pm3d, 1+cos(u)+.5*cos(u)*cos(v),.5*sin(v),sin(u)+.5*sin(u)*cos(v) with pm3d
