# Final for MATH3013 - PDE solver
Operation system: Win10

Compiler: MinGW

# Features
* Various numerical schemes for 1D heat equation, 2D heat equation and 1D wave equation

# Visualization
The numerical solutions are visualized by gnuplot.

Specific scheme for 1D heat equation:
```
filename = 'Forward Euler 1D'
set terminal gif animate
set output filename.'.gif'
stats filename.'.dat' nooutput
set xrange[0:1]
set yrange[0:120]
do for [i=0:int(STATS_blocks)-2] {
	plot filename.'.dat' using 2:3 index i
}
```
All schemes for 1D heat equation:
```
set terminal gif animate
set output '1D heat.gif'
stats 'Forward Euler 1D.dat' nooutput
set xrange[0:1]
set yrange[0:120]
do for [i=0:int(STATS_blocks)-2] {
	plot 'Forward Euler 1D.dat' using 2:3 index i, \
	'Backward Euler 1D.dat' using 2:3 index i, \
	'Richardson.dat' using 2:3 index i, \
	'Dufort-Frankel.dat' using 2:3 index i, \
	'Crank-Nicolson 1D.dat' using 2:3 index i
}
```
Specific scheme for 2D heat equation:
```
filename = 'Forward Euler 2D'
set terminal gif animate
set output filename.'.gif'
stats filename.'.dat' nooutput
set xrange[0:1]
set yrange[0:1]
set zrange[0:120]
do for [i=0:int(STATS_blocks)-2] {
	splot filename.'.dat' using 2:3:4 index i
}
```
Specific scheme for 1D wave equation:
```
filename = 'Explicit Finite Difference'
set terminal gif animate
set output filename.'.gif'
stats filename.'.dat' nooutput
set xrange[0:1]
set yrange[-2:2]
do for [i=0:int(STATS_blocks)-2] {
	plot filename.'.dat' using 2:3 index i
}
```
All schemes for 1D wave equation:
```
set terminal gif animate
set output '1D wave.gif'
stats 'Explicit Finite Difference.dat' nooutput
set xrange[0:1]
set yrange[-2:2]
do for [i=0:int(STATS_blocks)-2] {
	plot 'Explicit Finite Difference.dat' using 2:3 index i, \
	'Implicit Finite Difference.dat' using 2:3 index i, \
	'Exact solution.dat' using 2:3 index i
}
```