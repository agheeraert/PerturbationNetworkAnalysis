from CreatePerturbationNetwork import CreatePerturbationNetwork
import os 
from os.path import join

cwd = os.getcwd()
data_dir = join(cwd, 'data')

for directory in os.listdir(data_dir):
    dir_path = join(data_dir, directory)
    file1 = None
    for filename in os.listdir(dir_path):
        if filename[-4:] == '.pdb':
            if file1 is None:
                file1=filename
            else:
                file2=filename

    for w in range(40, 41, 2):
        CreatePerturbationNetwork(pdb1=join(dir_path, file1), pdb2=join(dir_path, file2)).draw_perturbation(threshold=w, output=join(dir_path, "w_"+str(w)+"_pnet.pdf"), rearrange=('H', 'F'))