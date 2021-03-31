#!/bin/bash

# converts positions data file to gnuplot format
# for plotting via the 'splot' command

if [[ -e $1 ]]; then
	# remove first line and count columns
	columns=$(sed '1d' $1|head -n1|wc -w)
	# remove first line and print rearranged data to stdout
	sed '1d' $1 | awk -v cols=$columns '{for (i=0;i<(cols-1)/3;i++) {print $(3*i+1+1), $(3*i+2+1), $(3*i+3+1);}; print ""; print ""}'
fi
