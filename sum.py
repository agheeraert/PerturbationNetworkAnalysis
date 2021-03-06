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

parser = argparse.ArgumentParser(description='Get a network and outputs only the interface(s)')
parser.add_argument('f',  type=str, nargs='+',
                    help='List of network files')

args = parser.parse_args()

for f in args.f:
    net = AANetwork()
    net.net = nx.read_gpickle(f)
    print(f, net.sum())


