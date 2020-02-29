## This file is for generating plots from the resulted data of the simulation data in GNUPLOT
## Created by Jonathan Henning

## Set the output to a PNG file (comment to see directly in screen)
set terminal pngcairo dashed size 800,800

# Remove the subtitle
unset key

set size square

## Number of particles in the simulation
part = 1340

## Do a loop over the data (99 lines in the case) saving the plots with the current number in the "pict/" folder (that have to be created first)
do for [k=0:99:1] {
	## Comment the output to file if you want to see directly in the screen
	set output 'pict/'.k.'.png'
	set label 1 "".k at -0.01,0.0335 tc lt 8
	plot [-0.01373:0.01373] [-0.001:0.02646] "data2.dat" every ::(k*part)::((k+1)*part-1) with points pt 6 ps 1.5
	## If you want to see directly in the screen is better to pause or it will be very quick
	#pause 0.5
}