from CreateNetwork import CreateNetwork

L_cutoffs = list(range(4, 10))

for i in L_cutoffs:
    CreateNetwork(cutoff=i).save_avg(folder='/home/agheerae/Python/PerturbationNetworkAnalysis/data/sim1/apo/', output='/home/agheerae/results/sim1/apo_nets/'+str(i)+'.p')
    CreateNetwork(cutoff=i).save_avg(folder='/home/agheerae/Python/PerturbationNetworkAnalysis/data/sim1/prfar/', output='/home/agheerae/results/sim1/prfar_nets/'+str(i)+'.p')
    CreateNetwork(cutoff=i).save_avg(folder='/home/agheerae/Python/PerturbationNetworkAnalysis/data/apo_all/', output='/home/agheerae/results/all_simu/apo_nets/'+str(i)+'.p')
    CreateNetwork(cutoff=i).save_avg(folder='/home/agheerae/Python/PerturbationNetworkAnalysis/data/prfar_all/', output='/home/agheerae/results/all_simu/prfar_nets/'+str(i)+'.p')
