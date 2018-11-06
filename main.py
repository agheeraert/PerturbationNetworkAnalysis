from CreatePerturbationNetwork import CreatePerturbationNetwork
import os 
from os.path import join
import argparse

# cwd = os.getcwd()
# data_dir = join(cwd, 'data')

# for directory in os.listdir(data_dir):
#     dir_path = join(data_dir, directory)
#     file1 = None
#     for filename in os.listdir(dir_path):
#         if filename[-4:] == '.pdb':
#             if file1 is None:
#                 file1=filename
#             else:
#                 file2=filename

#     for w in range(40, 41, 2):
#         CreatePerturbationNetwork(pdb1=join(dir_path, file1), pdb2=join(dir_path, file2)).draw_perturbation(threshold=w, output=join(dir_path, "w_"+str(w)+"_pnet.pdf"), rearrange=('H', 'F'))

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
parser.add_argument('-rearrange',  type=str, nargs=2,
                    help='display the full network so that two chains are separated')
parser.add_argument('-save',  type=str,
                    help='pickles the network')                      


args = parser.parse_args()
if args.range:
    threshold = range(args.range[0], args.range[1], args.range[2])
if args.threshold:
    threshold = args.threshold
if args.rearrange:
    rearrange = tuple(args.rearrange)
    
if args.avg:
    avg = True
else:
    avg = False

CreatePerturbationNetwork(path1=args.path1, path2=args.path2, avg=avg).draw_perturbation(threshold=threshold, output=args.output, rearrange=args.rearrange, save=args.save)
