"""Adapted from Lorenza Pacini's code 
https://github.com/lorpac

Main file to compute perturbation network
"""

from CreatePerturbationNetwork import CreatePerturbationNetwork
import os 
from os.path import join as jn
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Create the perturbation network between two proteins')
parser.add_argument('path1',  type=str,
                    help='Input file/folder of the unperturbed state')
parser.add_argument('path2',  type=str,
                    help='Input file/folder of the perturbed state')
parser.add_argument('output',  type=str,
                    help='Folder where to put the results')  
parser.add_argument('-avg',  type=str,
                    help='Does the average of the input on folders')
parser.add_argument('-cutoffs', type=int, nargs='+', default=[5],
		    help='Set a list of interaction cutoffs (in Angstrom) (can be only one number)')    
parser.add_argument('-drawing_method',  type=str,
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-pdb_path',  type=str,
                    help='PDB structure file to help draw the network (works only with the IGPS drawing method as of now)')
parser.add_argument('-drawing_colors',  type=list, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')


args = parser.parse_args()
    
if args.avg:
    avg = True
else:
    avg = False

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

for cutoff in args.cutoffs:
    out_dir = jn(args.output, 'cutoff'+str(cutoff))
    mkdir(out_dir)
    CreatePerturbationNetwork(path1=args.path1, path2=args.path2, avg=avg, cutoff=cutoff).draw_perturbation(out_dir, pdb_path=args.pdb_path, method=args.drawing_method, colors=args.drawing_colors)
