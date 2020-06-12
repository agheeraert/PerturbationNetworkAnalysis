from CreateNetwork import AANetwork
from DrawNetwork import DrawNetwork
from multiprocessing import Pool
import networkx as nx
from tqdm import tqdm
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import itertools
import os
import numpy as np
from os.path import join as jn, dirname, isdir, basename
from os import listdir
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

plt.style.use('classic')

class TimePertNet(AANetwork):
    """Creates the Time Perturbation Network"""
    def __init__(self, cutoff=5, n_procs=4, outdir=None):
        super().__init__(cutoff=cutoff)
        self.n_procs = n_procs
        self.outdir = outdir

    def create_timed(self, folder):
        L_mat = []
        for filepath in tqdm(listdir(folder)):
            net = self.create(jn(folder, filepath))
            mat = np.array(nx.to_numpy_matrix(net))
            if self.operation == 'i':
                L_mat.append(self.sum_interfaces(mat, net))
            elif self.operation == 'nw':
                L_mat.append(self.nw(mat))
        return np.stack(L_mat, axis=-1), net
    
    def smart_loader(self, path):
        if len(basename(path)) == 0: #should get the name of a folder with /
            output = jn(self.outdir, basename(path[:-1])+'.npy')
        else: #should get all the rest
            output = jn(self.outdir, basename(path).split('.')[0]+'.npy')
        if isdir(path): 
            assert output, print('no output path specified !')
            #loads a folder and do the average aanet, saves
            mat, net = self.create_timed(path)
            np.save(output, mat)
            nx.write_gpickle(net, output.replace('.npy', '.net'))
        elif path[-4:] == '.npy': 
            #directly loads a precomputed network
            mat = np.load(path)
            net = nx.read_gpickle(path.replace('.npy', '.net'))
        else:
            print('Unknown extension for file %s' %path)
            return None
        return mat, net
    
    def create_comp(self, L_comp, operation):
        self.operation = operation
        pool = Pool(processes=self.n_procs)
        L_res = pool.map(self.smart_loader, L_comp)
        self.contact_diffs = [elt[0] for elt in L_res]
        self.nettop = L_res[0][1]

    # def sum_interfaces(self, limit=None):
    #     assert self.L_timemat and self.nettop, print('No comparison has been computed')
    #     if not limit:
    #         firstchain = next(iter(self.nettop.nodes()))[-1]
    #         limit = next(i for i, x in enumerate(self.nettop.nodes()) if x[-1] != firstchain)
    #     self.contact_diffs = []
    #     for timemat in self.L_timemat:
    #         self.contact_diffs.append((np.add(np.sum(timemat[limit:,:limit], axis=(0,1)),np.sum(timemat[:limit,limit:], axis=(0,1))))/2)

    def sum_interfaces(self, mat, net):
        firstchain = next(iter(net.nodes()))[-1]
        limit = next(i for i, x in enumerate(net.nodes()) if x[-1] != firstchain)
        contact_diff = (np.add(np.sum(mat[limit:,:limit], axis=(0,1)),np.sum(mat[:limit,limit:], axis=(0,1))))/2
        return contact_diff
    
    def nw(self, mat):
        degree = np.sum((mat!=0), axis=0) #count non zeros to get degree
        weight = np.sum(mat, axis=0)
        return np.divide(weight, degree)

    def moving_avg(self, table, n=10):
        convo = np.convolve(table, np.ones((n+1,))/(n+1), mode='valid')
        nanarray = np.empty(n)
        nanarray[:] = np.nan
        return np.concatenate([convo, nanarray], axis=-1)

    def plot_interfaces(self, output, factor=None, unit=None, labels=None, moving=1):
        colors = itertools.cycle(('black', 'g', 'b', 'r', 'magenta', 'orange'))
        labels = iter(labels)
        f = plt.figure()
        for contact_diff in self.contact_diffs:
            if moving != 1:
                contact_diff = self.moving_avg(contact_diff, moving)
            color = next(colors)
            label = next(labels)
            if factor:
                plt.plot(np.arange(0, len(contact_diff)*factor, factor), contact_diff, label=label, color=color)
            else:
                plt.plot(contact_diff, color=color)
        if unit:
            plt.xlabel('Time (%s)' %unit)
        else:
            plt.xlabel('Frame number')
        plt.ylabel('Number of contacts at interface')
        lims = plt.gca().get_ylim()
        plt.ylim((lims[0], lims[1]+20))
        plt.legend(ncol=2)
        plt.tight_layout()
        plt.savefig(output)
    
    def covmatrix(self, output=None):
        neighborhood_watches = np.concatenate(self.contact_diffs, axis=-1)
        covmat = np.cov(neighborhood_watches)
        if output:
            np.save(output, covmat)
        self.covmat = covmat
    
    def plot_covmatrix(self, output):
        f = plt.figure()
        firstchain = next(iter(self.nettop.nodes()))[-1]
        limit = next(i for i, x in enumerate(self.nettop.nodes()) if x[-1] != firstchain)
        n_aa = len(self.nettop.nodes())
        ticks = [list(range(50,limit+1,50))+list(range(limit+50,n_aa+1,50)),list(range(50,limit+1,50))+list(range(50,n_aa+1-limit,50))]
        plt.xticks(*ticks)
        plt.yticks(*ticks)
        plt.plot([limit+0.5, limit+0.5], [0,n_aa+1], linestyle=':', linewidth=1, color='k')
        plt.imshow(self.covmat, origin='lower',cmap='jet', interpolation='none')
        plt.colorbar()
        plt.savefig(output, dpi=1000)
        plt.close()

    def topk_cov_matrix(self, k, reverse=False):
        N, N = self.covmat.shape
        if reverse:
            _indices = np.argsort(self.covmat, axis=None)[:2*k:2] 
        else:
            _indices = np.argsort(self.covmat, axis=None)[-1:-2*k:-2]
        indices = np.stack(np.divmod(_indices, N)).transpose()
        resid2name = {i: aa for i, aa in enumerate(self.nettop.nodes())}
        print(np.vectorize(resid2name.get)(indices))
        for elt in self.covmat.reshape(-1)[_indices]:
            print(elt)




        



