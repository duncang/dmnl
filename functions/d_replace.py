#!/usr/bin/env python

# D replace function written by Duncan Greer 12 April 2007
#
# $Id: d_replace.py 807 2007-04-12 23:03:11Z greerd $
# $Log$
#

import os,sys,re

usage = "usage: %s infile outfile" %         os.path.basename(sys.argv[0])

if len(sys.argv) < 3:
	print usage
else:

	infilename = sys.argv[1]
	outfilename = sys.argv[2]

	infile = open(infilename,'r')
	outfile = open(outfilename,'w')

	# compile pattern
	#pattern = re.compile('D')

#	for line in infile:
#		i = pattern.finditer(line)
#		for match in i

	currentline = 0

	for line in infile:
		currentline = currentline + 1
		if currentline == 3:
			outfile.write(line.replace('D','e'))
		elif currentline == 4:
			outfile.write(line.replace('D','e'))
		elif currentline == 5:
			outfile.write(line.replace('D','e'))
		elif currentline > 7:
			outfile.write(line.replace('D', 'e'))
		else:
			outfile.write(line)
	
	infile.close()
	outfile.close()
