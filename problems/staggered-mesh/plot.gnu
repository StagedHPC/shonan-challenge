set term png  size 1024,768 


set pm3d
set pm3d map
set grid

set cbrange [-0.01:0.01]

set palette model RGB defined (-0.01 "blue", 0 "white", 0.01 "red")


set size ratio -1

set out "debug.png"

splot "debug.txt" t ''

