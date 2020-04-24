"""Adapted from Lorenza Pacini's code 
https://github.com/lorpac

Main file to compute perturbation network
"""

from CreatePerturbationNetwork import PerturbationNetwork
import os 
from os.path import join as jn, exists
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Create the perturbation network between two proteins')
parser.add_argument('path1',  type=str,
                    help='Input file/folder of the unperturbed state')
parser.add_argument('path2',  type=str,
                    help='Input file/folder of the perturbed state')
parser.add_argument('output',  type=str,
                    help='Folder where to put the results')  
# parser.add_argument('-avg',  type=bool, default=False,
#                     help='Does the average of the input on folders')
# parser.add_argument('-std',  type=bool, default=False,
#                     help='Does the std modulated by avg of the input on folders')
parser.add_argument('-c', type=int, nargs='+', default=[5],
		    help='Set a list of interaction cutoffs (in Angstrom) (can be only one number)')    
parser.add_argument('-dm',  type=str, default='default',
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-p',  type=str,
                    help='PDB structure file to help draw the network (works only with the IGPS drawing method as of now)')
parser.add_argument('-dc',  type=str, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')
parser.add_argument('-af',  type=str, nargs='+', default=None,
                    help='Toggles the Alignment mode and align between flags')


args = parser.parse_args()

if args.dm != 'default':
    assert exists(args.p), '%s not found' %args.p

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

if len(args.c) > 1:
    for cutoff in args.c:
        out_dir = jn(args.output, 'cutoff'+str(cutoff))
        mkdir(out_dir)
        pn = PerturbationNetwork(path1=args.path1, path2=args.path2, out_dir=out_dir, cutoff=cutoff)
else:
    mkdir(args.output)
    pn = PerturbationNetwork(path1=args.path1, path2=args.path2, out_dir=args.output, cutoff=args.c[0])

if args.af:
    pn.align(args.af)

pn.draw_perturbation(pdb_path=args.p, method=args.dm, colors=args.dc)    
