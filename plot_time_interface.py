from TimePertNet import TimePertNet
import os 
from os.path import join as jn, exists, basename, isdir, dirname
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Creates a timed plot of the evolution of number of contacts at the interface')

parser.add_argument('-f',  type=str, nargs='+',
                    help='List of files to load')
parser.add_argument('-o',  type=str,
                    help='Output path for the figure')
parser.add_argument('-x',  type=float,
                    help='Time scale factor between the number of frames and the selected unit')
parser.add_argument('-u',  type=str,
                    help='Time unit')
parser.add_argument('-p',  type=str,
                    help='Operation to perform: (i) interfaces')
parser.add_argument('-l',  type=str, nargs='+',
                    help='List of labels')
parser.add_argument('-m',  type=int, default=1,
                    help='Moving average of specified time-window will be performed on the final plot')

args = parser.parse_args()

def smart_label_generator(files):
    labels = []
    for f in files:
        if len(basename(f))==0: 
            labels.append(basename(f[:-1]))
        else:
            labels.append(basename(f).split('.')[0])
    return labels

if not args.l:
    args.l = smart_label_generator(args.f)

tpn = TimePertNet(outdir=dirname(args.o))
tpn.create_comp(args.f, args.p)
if args.x and args.u:
    tpn.plot_interfaces(args.o, factor=args.x, unit=args.u, labels=args.l, moving=args.m)
else:
    tpn.plot_interfaces(args.o, labels=args.l, moving=args.m)
