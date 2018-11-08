from CreatePerturbationNetwork import CreatePerturbationNetwork
import os 
from os.path import join
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Create the perturbation network of two proteins')
parser.add_argument('path1',  type=str,
                    help='First input file/folder')
parser.add_argument('path2',  type=str,
                    help='Second input file/folder')
parser.add_argument('output',  type=str,
                    help='Output folder for the results')  
parser.add_argument('-avg',  type=str,
                    help='does the average of the input on folders')
parser.add_argument('-range',  type=float, nargs=3,
                    help='create the perturbation network for a range of thresholds')
parser.add_argument('-threshold',  type=float,
                    help='create the perturbation network for one threshold')  
parser.add_argument('-cutoff', type=float, default=5,
		    help='set the interaction cutoff (in Angstrom)')
parser.add_argument('-rearrange',  type=str, nargs=2,
                    help='display the full network so that two chains are separated')
parser.add_argument('-save',  type=str, default=True,
                    help='pickles the network')                      


args = parser.parse_args()
if args.range:
    threshold = np.arange(args.range[0], args.range[1], args.range[2]).tolist()
if args.threshold:
    threshold = args.threshold
if args.rearrange:
    rearrange = tuple(args.rearrange)
    
if args.avg:
    avg = True
else:
    avg = False

CreatePerturbationNetwork(path1=args.path1, path2=args.path2, avg=avg, cutoff=args.cutoff).draw_perturbation(threshold=threshold, output=args.output, rearrange=args.rearrange, save=args.save)
