from allopathway import allopathway_nodes
from os import path, listdir
from CreatePerturbationNetwork import CreatePerturbationNetwork

traj_path = 'data/sim1/'
apo_path = path.join(traj_path, 'apo')
prfar_path = path.join(traj_path. 'prfar')

for _frame in listdir(apo_path):
    apo_frame = path.join(apo_path, _frame)
    prfar_frame = path.join(prfar_frame, _frame)
    CreatePerturbationNetwork(apo_frame, prfar_frame).create()

