set terminal postscript font "Times,16" ps 2
set output 'Is-2D_2.eps'
set grid
set pm3d map
set size square
set contour base
unset k
set xrange [-600:600]
set yrange [-80:80]
set cbrange [-200:200]
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
set title '{/Symbol e}_{/Symbol \255}=-5,{/Symbol e}_{/Symbol \257}=5,{/Symbol G}_L=1,{/Symbol m}=0,T_0=300K'
set xlabel '{/Symbol D}T'
set ylabel '{/Symbol m}_{/Symbol \257}-{/Symbol m}_{/Symbol \255}'
set zlabel 'I_s'
splot 'current2.txt' t 'I_s'
set output
q