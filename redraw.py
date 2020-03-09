from CreatePerturbationNetwork import CreatePerturbationNetwork
import os 
from os.path import dirname 
import argparse
import numpy as np
from DrawNetwork import DrawNetwork
import networkx as nx

parser = argparse.ArgumentParser(description='Create the perturbation network between two proteins')
parser.add_argument('f',  type=str, 
                    help='Input pickled network (should be 0.p in most cases) ')
parser.add_argument('-drawing_method',  type=str, default='default',
                    help='Method used to draw the graphs. Default = Networkx default. IGPS = IGPS splitting.')
parser.add_argument('-pdb_path',  type=str,
                    help='PDB structure file to help draw the network (works only with the IGPS drawing method as of now)')
parser.add_argument('-drawing_colors',  type=str, nargs=2, default=['red', 'dodgerblue'],
                    help='Color used to draw the edges')
parser.add_argument('-increment',  type=float, default=1,
                    help='Increment step to redraw')

args = parser.parse_args()
    
net = nx.read_gpickle(args.f)
output = dirname(args.f) 
DrawNetwork(net, output, pdb_path=args.pdb_path, method=args.drawing_method, colors=args.drawing_colors, increment=args.increment)    
