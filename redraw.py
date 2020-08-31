from CreatePerturbationNetwork import PerturbationNetwork
import os 
from os.path import dirname 
import argparse
import numpy as np
from DrawNetwork import DrawNetwork
import networkx as nx

parser = argparse.ArgumentParser(description='Create the perturbation network between two proteins')
parser.add_argument('f',  type=str, 
                    help='Input pickled network (should be 0.p in most cases) ')
parser.add_argument('-dm',  type=str, default='default',
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-p',  type=str,
                    help='PDB structure file to help draw the network (works only with the IGPS drawing method as of now)')
parser.add_argument('-dc',  type=str, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')
parser.add_argument('-i',  type=float, default=1,
                    help='Increment step to redraw')
parser.add_argument('-oxy',  type=float, default=None, nargs=9,
                    help='the 9 coordinates needed for OXY method (3 coordinates for each point of the frame)')
parser.add_argument('-c',  type=str, default=None,
                    help='Name of the top chain to separate the interface')

args = parser.parse_args()
    
net = nx.read_gpickle(args.f)
output = dirname(args.f) 
DrawNetwork(net, output, pdb_path=args.p, method=args.dm, colors=args.dc, increment=args.i, L_OXY=args.oxy, chaintop=args.c)    
