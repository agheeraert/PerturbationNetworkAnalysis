import argparse
import numpy as np
from CreateNetwork import CreateNetwork
import matplotlib.pyplot as plt
import os.path
from scipy.stats import linregress

# parser = argparse.ArgumentParser(description='Screen the number of degree and weight of each node and compare it within the different states')
# parser.add_argument('path1',  type=str,
#                     help='Static input file path of reference (apo)')
# parser.add_argument('path2',  type=str,
#                     help='Static input file path to compare (complex)')
# parser.add_argument('output',  type=str,
#                     help='Output directory for the results')
# parser.add_argument('range',  type=float, nargs=3,
#                     help='Range within which to screen the cutoff')
# args = parser.parse_args()
def r2(x, y):
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    return str(r_value)[:5]

class args:
    path1 = "/home/agheerae/PDB/Apo_frames/Sim1/frame_1.pdb"
    path2 = "/home/agheerae/PDB/Prfar_frames/Sim1/frame_1.pdb"
    output = "/home/agheerae/Python/PerturbationNetworkAnalysis/results/screen/"
    range = [3, 10, 1]

L_nets_apo, L_nets_bound = [], []
L_cutoffs = np.arange(args.range[0], args.range[1], args.range[2]).tolist()
for elt in L_cutoffs:
    L_nets_apo.append(CreateNetwork(cutoff=elt).create(args.path1))
    L_nets_bound.append(CreateNetwork(cutoff=elt).create(args.path2))

for node in L_nets_apo[0].nodes():
    L_degrees_apo, L_weights_apo, L_degrees_bound, L_weights_bound  = [], [], [], []
    for i, elt in enumerate(L_cutoffs):
        L_degrees_apo.append(L_nets_apo[i].degree(node))
        L_degrees_bound.append(L_nets_bound[i].degree(node))
        L_weights_apo.append(L_nets_apo[i].degree(node, weight='weight'))
        L_weights_bound.append(L_nets_bound[i].degree(node, weight='weight'))
    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax1.scatter(L_cutoffs, L_degrees_apo, marker='x', c='b', label='degree in apo, R^2= '+r2(L_cutoffs, L_degrees_apo))
    ax1.scatter(L_cutoffs, L_degrees_bound, c='r', marker='+', label='degree in PRFAR-bound, R^2= '+r2(L_cutoffs, L_degrees_bound))
    f.suptitle("Degree distribution for node " + node)
    ax1.set_xlabel('Cutoff (in A)')
    ax1.set_ylabel('Degree of the node')
    ax1.legend(loc='upper left')
    plt.savefig(os.path.join(args.output, node+'_degree.pdf'))
    f = plt.figure()
    ax2 = f.add_subplot(111)
    ax2.scatter(L_cutoffs, L_weights_apo, c='b', marker='x', label='weight in apo, R^2= '+r2(L_cutoffs, L_weights_apo))
    ax2.scatter(L_cutoffs, L_weights_bound, c='r', marker='+', label='weight in PRFAR-bound, R^2= '+r2(L_cutoffs, L_weights_bound))
    ax2.set_xlabel('Cutoff (in A)')
    ax2.set_ylabel('Weight of the node')
    ax2.legend(loc='upper left')
    plt.savefig(os.path.join(args.output, node+'_weight.pdf'))
