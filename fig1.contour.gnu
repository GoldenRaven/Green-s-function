set terminal pdf font "Times,12" ps 0.5
set output 'Is.pdf'
set grid
set multiplot layout 1,2
set title '{/Symbol D}T=T_L-T_R,T_0=300K,{/Symbol e}_{/Symbol \255}=10,{/Symbol e}_{/Symbol \257}=15,{/Symbol G}_L=6'
# set xtics 0.2
set ylabel 'I_s'
set xlabel '{/Symbol D}T'
set xrange [-600:600]
set k right top
plot 'I-dT.txt' u 1:2 w lp t '{/Symbol D}{/Symbol m}=0,{/Symbol m}_0=0'
set title '{/Symbol D}{/Symbol m}_s={/Symbol m}_{down}-{/Symbol m}_{up}, {/Symbol m}_0={/Symbol m}_R=0'#,{/Symbol m}_{down}={/Symbol m}_0+{/Symbol D}{/Symbol m}_s/2,{/Symbol m}_{up}={/Symbol m}_0-{/Symbol D}{/Symbol m}_s/2'
set xlabel '{/Symbol m}_{down}-{/Symbol m}_{up}'
set xrange [-30:30]
plot 'I-mu.txt' u 1:2 w lp t '{/Symbol D}T=0, T_0=300K'
# set yrange [0:10]
unset multiplot
set output

set terminal postscript font "Times,16" ps 2
set output 'Is-2D.eps'
set pm3d map
set size square
set contour base
unset k
set xrange [-600:600]
set yrange [-80:80]
set cbrange [-150:150]
set cntrparam levels incremental 0,1,0
# set palette defined (0 "black",\
#                      0.2 "black",\
#                      0.2 "red",\
#                      0.4 "red",\
#                      0.4 "orange-red",\
#                      0.6 "orange-red",\
#                      0.6 "orange",\
#                      0.8 "orange",\
#                      0.8 "yellow",\
#                      1.0 "yellow",\
#                      1.0 "light-green",\
#                      1.2 "light-green")

set palette defined (-1 "blue", 0 "white", 1 "red") #thermometer palette
# set palette rgbformulae 33,13,10
# set palette defined
set title '{/Symbol e}_{/Symbol \255}=10,{/Symbol e}_{/Symbol \257}=15,{/Symbol G}_L=6,T_0=300K'
set xlabel '{/Symbol D}T'
set ylabel '{/Symbol m}_{/Symbol \257}-{/Symbol m}_{/Symbol \255}'
set zlabel 'I_s'
splot 'current1.txt' t 'I_s'
set output
q