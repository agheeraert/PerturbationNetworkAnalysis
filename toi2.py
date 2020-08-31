import networkx as nx
import argparse
import matplotlib.pyplot as plt 
import numpy as np
import pickle as pkl

parser = argparse.ArgumentParser(description='Graph the type of interaction for type II')
parser.add_argument('f',  type=str, nargs=2,
                    help='Input networks [HYDROPHOBIC 0.p] [POLAR 0.p]')
parser.add_argument('o',  type=str,
                    help='Output path')                   
parser.add_argument('-s',  type=str, default=1,
                    help='Step to increase the threshold') 
args = parser.parse_args()


hydro_net = nx.read_gpickle(args.f[0])
polar_net = nx.read_gpickle(args.f[1])

def is_saltbridge(u, v):
    u, v = u[0], v[0]
    if u in ['R', 'K']:
        if v in ['D', 'E']:
            return True
    elif u in ['D', 'E']:
        if v in ['R', 'K']:
            return True
    return False

def is_samecharge(u, v):
    u, v = u[0], v[0]
    if u in ['R', 'K'] and v in ['R', 'K']:
            return True
    elif u in ['D', 'E'] and v in ['D', 'E']:
            return True

def count_polar(net, polar_cont, sb_cont, sc_cont):
    c_p, c_sb, c_sc = 0, 0, 0
    for u, v in net.edges():
        if is_saltbridge(u, v):
            c_sb+=1
        elif is_samecharge(u, v):
            c_sc+=1
        else:
            c_p+=1
    polar_cont.append(c_p), sb_cont.append(c_sb), sc_cont.append(c_sc)

def cut_net(net, threshold):
    to_remove_edges = []
    for u, v in net.edges():
        weight = net.get_edge_data(u, v)['weight']
        if weight < threshold:
            to_remove_edges.append((u, v))
    net.remove_edges_from(to_remove_edges)
    net.remove_nodes_from(list(nx.isolates(net)))
    return net


hydro_cont, polar_cont, sbp_cont, sbh_cont, scp_cont, sch_cont = [], [], [], [], [], []

rh, rp = True, True
count_polar(hydro_net, hydro_cont, sbh_cont, sch_cont)
count_polar(polar_net, polar_cont, sbp_cont, scp_cont)
threshold = args.s
_hnet = cut_net(hydro_net, threshold)
_pnet = cut_net(polar_net, threshold)

while rh or rp:
    if rh:
        count_polar(_hnet, hydro_cont, sbh_cont, sch_cont)
    else:
        hydro_cont.append(0)
        sbh_cont.append(0)  
        sch_cont.append(0)      
    if rp:
        count_polar(_pnet, polar_cont, sbp_cont, scp_cont)
    else:
        polar_cont.append(0)
        sbp_cont.append(0)
        scp_cont.append(0)
    threshold += args.s
    _hnet = cut_net(_hnet, threshold)
    _pnet = cut_net(_pnet, threshold)
    if len(_hnet.edges()) == 0:
        rh = False
    if len(_pnet.edges()) == 0:
        rp = False
f = plt.figure()

hydro_cont, polar_cont, sbp_cont, sbh_cont, sch_cont, scp_cont = np.array(hydro_cont), np.array(polar_cont), np.array(sbp_cont),np.array(sbh_cont), np.array(sch_cont), np.array(scp_cont)
total = hydro_cont + polar_cont + sbp_cont + sbh_cont + sch_cont + scp_cont

plt.plot(hydro_cont/total*100, c='b', label='Hydrophobic contacts')
plt.plot(polar_cont/total*100, c='purple', label='Polar contacts')
plt.plot(sbp_cont/total*100, c='yellow', label='Salt bridge (polar)')
plt.plot(sbh_cont/total*100, c='orange', label='Salt bridge (hydrophobic)')
plt.plot(scp_cont/total*100, c='g', label='Same charge (polar)')
plt.plot(sch_cont/total*100, c='cyan', label='Same charge (hydrophobic)')
plt.legend()
plt.savefig(args.o)







