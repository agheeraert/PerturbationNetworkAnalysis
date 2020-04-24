import networkx as nx
import numpy as np
from Bio.PDB.Polypeptide import aa1, aa3
from Bio.PDB import PDBParser
from tqdm import tqdm
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)
import argparse

from CreateNetwork import three2one, AANetwork

parser = argparse.ArgumentParser(description='Creates the node weight line plot')
parser.add_argument('f',  type=str, nargs='+',
                    help='List of network files')
parser.add_argument('m',  type=str, nargs=1,
                    help='Name of the method used (all) all, (sum) simple sum, (n1) norm1, (n2), norm2')                   

args = parser.parse_args()

if args.m == ['all']:
    m = ['sum', 'n1', 'n2']

for f in args.f:
    net = AANetwork()
    net.net = nx.read_gpickle(f)
    assert '.p' in f, print('Should use the .p extension')
    for method in m:
        net.node_weigths_line(method, f.replace('.p', '_node_weights%s.png' %method))