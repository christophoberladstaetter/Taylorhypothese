#gnuplot test
reset

# Ausgabeformat:
set term x11
#set term postscript
#set term png

#set output "ygauss-nz.ps"
##set view equal xy
#unset colorbox
set title " "

#filename = "Eky.dat"  ## flow potential phi (x,y)
filename = "w2d.dat"  ## vorticity omega (x,y)
#filename = "v2d.dat"  ## y velocity (x,y)


# einfacher 2D Farbplot
#set xrange [0 : 1023] 
#set yrange [0 : 1023] 
#set xrange [0 : 255] 
#set yrange [0 : 255] 
#set cbrange [-1.1 : 1.1] 
#set size square

#set palette defined (-3 "black", -2 "blue", -1 "cyan", 0 "green", 1 "yellow", 2 "red", 3 "black")
set palette defined (-0.5 "black", -0.025 "blue", -0.01 "cyan", -0.005 "green",-0.0000025 "light-green",-0.0000000005 "gray",0 "white", 0.000000005 "honeydew", 0.0000025 "lemonchiffon", 0.0005 "yellow", 0.01 "orange", 0.025 "red", 0.5 "black")

#set palette defined (-1 "blue", 0 "white", 1 "red")
#set palette defined (-1.2 "black", -1 "blue", 0 "white", 1 "red", 1.2
#"black")
#set palette defined (-2 "black", -1 "blue", 0 "white", 1 "red", 2 "black")
#set palette defined (-1.25 "black", -0.5 "blue", 0 "white", 0.5 "orange", 1.25 "red", 3.5 "black")

set pm3d map
splot filename notitle
#splot filename every 4 notitle

# 2D Farbplot mit Konturen
#
#set contour base; set cntrparam level 30
#unset surface
#set table 'contour.dat'
#splot filename
#unset table
#!awk "NF<2{printf\"\n\"}{print}" <contour.dat >contour1.dat
#reset
#
#set ticslevel 0
#set palette defined (-1.25 "black", -1 "blue", 0 "white", 1 "red", 1.25 "black")
##set xrange [0 : 127] 
##set yrange [0 : 127] 
##set size square
#set cbrange [0 : 2] 
#set pm3d map
#splot filename notitle with pm3d, 'contour1.dat' notitle with line lt -1
#!rm contour.dat contour1.dat

pause 1
replot
reread



