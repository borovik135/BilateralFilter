#!/usr/bin/python

# python executable
import string
import getopt
import sys
import math
import fileinput
from random import random
from math import sqrt
import numpy as np

def main():

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o:v:",
                       ["help","input=","output=","p="])
    except getopt.GetoptError as err:
        print str(err) 
        usage()
        sys.exit(2)

    for o, a in opts:
       if o in ("-h", "--help"):
           usage()
       if o in ("-i", "--input"):
           input = str(a)
       if o in ("-o", "--output"):
           output = str(a)
       if o in ("-v", "--variance"):
           percentage = str(a)

    if len(sys.argv) == 1 :
        usage()
        sys.exit(2)

    v = float(percentage);
    if v <= 0:
        print "input variance should be a percentage (>0); exiting\n"
        sys.exit()
    elif v >= 100:
        print "input variance should be a percentage (<100); exiting\n"
        sys.exit()

    print "Determining size of point cloud in file %s"%input

    i =0;
    for line in fileinput.input(input):
	p = [float(r) for r in line.split()]
	if i == 0 :
	    xmin = p[0]
	    xmax = p[0]
	    ymin = p[1]
	    ymax = p[1]
	    zmin = p[2]
	    zmax = p[2]
	else :	
	    if p[0] < xmin :
                xmin = p[0]
	    elif p[0]>xmax :
	        xmax = p[0];
	    if p[1] < ymin :
                ymin = p[1]
	    elif p[1]>ymax :
	        ymax = p[1];
	    if p[2] < zmin :
                zmin = p[2]
	    elif p[2]>zmax :
	        zmax = p[2];
	i = i+1;

    sx = xmax - xmin
    sy = ymax - ymin
    sz = zmax - zmin
    diag = sqrt(sx*sx + sy*sy + sz*sz)

    print "%d points in the file"%i
    print "min : %f %f %f "%(xmin,ymin,zmin)
    print "max : %f %f %f "%(xmax,ymax,zmax)
    print "size: %f %f %f "%(sx,sy,sz)
    print "diagonal: %f"%(diag)


    variance = v*diag/100;
    print "adding noise with variance %f\n"%variance

    f = open(output,'w')
    for line in fileinput.input(input):
	p = [float(r) for r in line.split()]
        noise = np.random.normal(0.0, variance, 3);
        p[0] += noise[0]
        p[1] += noise[1]
        p[2] += noise[2]
        lineout='\t'.join(str(r) for r in p)
        f.write("%s\n"%lineout)

    f.close()


def usage():
    print 'USAGE:\n'
    print '-i <input file>'
    print '-o <output file>'
    print '-v <variance expressed as a percentage of the shape diagonal>'
    sys.exit(' ')


if __name__ == "__main__":
   main()

