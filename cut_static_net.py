import networkx as nx

L_cutoffs = list(range(3, 10))
L_thresh = list(range(1, 30))

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
    load_apo = '/home/agheerae/results/static/apo_'+str(i)+'.p'
    load_prfar = '/home/agheerae/results/static/prfar_'+str(i)+'.p'
    net_apo = nx.read_gpickle(load_apo)
    net_prfar = nx.read_gpickle(load_prfar)
    for j in L_thresh:
        output_apo = '/home/agheerae/results/static/apo_'+str(i)+'_'+str(j)+'.p'
        output_prfar = '/home/agheerae/results/static/prfar_'+str(i)+'_'+str(j)+'.p'
        nx.write_gpickle(cut(net_apo, j), output_apo)
        nx.write_gpickle(cut(net_prfar, j), output_prfar)
        


