"""Adapted from Lorenza Pacini's code 
https://github.com/lorpac

"""

from CreateNetwork import CreateNetwork
import os 
from os.path import join as jn
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Create the perturbation network between two proteins')
parser.add_argument('f',  type=str,
                    help='Input file/folder of the unperturbed state')
parser.add_argument('o',  type=str,
                    help='Name of the output')
parser.add_argument('-avg',  type=str,
                    help='Does the average of the input on folders')
parser.add_argument('-cutoffs', type=int, nargs='+', default=[5],
		    help='Set a list of interaction cutoffs (in Angstrom) (can be only one number)')    

args = parser.parse_args()
    
def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

for cutoff in args.cutoffs:
    if args.avg:
        CreateNetwork(cutoff=cutoff).save_avg(args.f, args.o) 
    else:
        CreateNetwork(cutoff=cutoff).save(args.f, args.o)

