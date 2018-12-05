import networkx as nx

L_cutoffs = list(range(6, 10))
L_thresh = list(range(30, 100))

def cut(net, thresh):
    to_remove_edges = []
    for u, v in net.edges():
        weight = net.get_edge_data(u, v)['weight']
        if weight < thresh:
            to_remove_edges.append((u, v))
    net.remove_edges_from(to_remove_edges)
    net.remove_nodes_from(list(nx.isolates(net)))
    return net



for i in L_cutoffs:
    load_apo = '/home/agheerae/results/sim1/nets/threshold/apo_'+str(i)+'_0.p'
    load_prfar = '/home/agheerae/results/sim1/nets/threshold/prfar_'+str(i)+'_0.p'
    net_apo = nx.read_gpickle(load_apo)
    net_prfar = nx.read_gpickle(load_prfar)
    for j in L_thresh:
        output_apo = '/home/agheerae/results/sim1/nets/threshold/apo_'+str(i)+'_'+str(j)+'.p'
        output_prfar = '/home/agheerae/results/sim1/nets/threshold/prfar_'+str(i)+'_'+str(j)+'.p'
        cut_apo = cut(net_apo, j)
        if len(cut_apo.nodes()) != 0:
                nx.write_gpickle(cut_apo, output_apo)
        cut_prfar = cut(net_prfar, j)
        if len(cut_prfar.nodes) !=0:
                nx.write_gpickle(cut_prfar, output_prfar)
        


