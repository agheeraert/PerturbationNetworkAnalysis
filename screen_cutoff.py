import argparse
import numpy as np
from CreateNetwork import CreateNetwork

parser = argparse.ArgumentParser(description='Screen the number of degree and weight of each node and compare it within the different states')
parser.add_argument('path1',  type=str,
                    help='Static input file path of reference (apo)')
parser.add_argument('path2',  type=str,
                    help='Static input file path to compare (complex)')
parser.add_argument('range',  type=float, nargs=3,
                    help='Range within which to screen the cutoff')

args = parser.parse_args()

L_nets = []
L_cutoffs = np.arange(args.range[0], args.range[1], args.range[2]).tolist()
for elt in L_cutoffs:
    L_nets.append(((CreateNetwork(cutoff=elt).create(args.path1)), CreateNetwork(cutoff=elt).create(args.path2)))

for node in L_nets[0][0].nodes():
    L_degrees_apo, L_weights_apo, L_degrees_bound, L_weights_bound  = [], [], [], []
    for elt in L_cutoffs:
        L_degrees_apo.append(L_nets[elt][0].degree(node))
        L_degrees_bound.append(L_nets[elt][1].degree(node))
        L_degrees_apo.append(L_nets[elt][0].degree(node, weight='weight'))
        L_degrees_bound.append(L_nets[elt][1].degree(node, weight='weight'))

