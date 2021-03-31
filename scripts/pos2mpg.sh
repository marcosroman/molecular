#!/bin/bash

folder=$1
prefix=$2

posfile=$folder/$prefix"pos.data"
infofile=$folder/$prefix"info.data"
# check that data files are available
if [[ ! -e $posfile ]]; then
	echo $posfile not found"
	exit
fi
if [[ ! -e $infofile ]]; then
	echo $infofile not found"
	exit
fi

# we need to reformat the <prefix>pos.data file:
# each line has the format: t x1[0] x1[1] x1[2] ... xn[0] xn[1] xn[3]
# (where xi[0], xi[1], xi[2] are the x, y, z coordinates for particle i)
# and to plot them in 3D with gnuplot we need them in the format
# x1[0] x1[1] x1[2]
# ...
# xn[0] xn[1] xn[2]
gpltposfile=$folder/$prefix"gpltpos.data"
./convpos2gpltpos.sh $posfile > $gpltposfile

# make dir for storing the intermediate image files
pngoutputfolder=$folder/$prefix
if [[ ! -e $pngoutputfolder ]]; then mkdir $pngoutputfolder; fi

# plot
gnuplot <(echo "gpltposdatafile=\""$gpltposfile"\"; outputfolder=\""$pngoutputfolder"\";"; echo $(sed 's/# //' $infofile | tr '\n' ';'); cat animate.gnuplot)

ffmpeg -i $pngoutputfolder/%d.png -c:v libx264 -framerate 40 -crf 20 $prefix.mpg

