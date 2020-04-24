"""Adapted from Lorenza Pacini's code 
https://github.com/lorpac

"""

from CreateNetwork import AANetwork
import os 
from os.path import join as jn
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Create the aa network of a signle proteins')
parser.add_argument('f',  type=str,
                    help='Input file/folder')
parser.add_argument('o',  type=str,
                    help='Output path')
parser.add_argument('-avg',  type=bool, default=False,
                    help='Does the average of the input on folders')
parser.add_argument('-cutoffs', type=int, nargs='+', default=[5],
		    help='Set a list of interaction cutoffs (in Angstrom) (can be only one number)')
parser.add_argument('-i', type=bool, default=False,
                     help='1 means to select only interface, 0 (default) the full graph')    

args = parser.parse_args()

for cutoff in args.cutoffs:
    if args.i:
        cnet = AANetwork(cutoff=cutoff)
        if args.avg:
            cnet.create_average(args.f)
        else:
            cnet.create(args.f)
        cnet.get_interface()
        cnet.threshold_loop(args.o, args.f)

    else:
        if args.avg:
            AANetwork(cutoff=cutoff).save_avg(args.f, args.o) 
        else:
            AANetwork(cutoff=cutoff).save(args.f, args.o)

