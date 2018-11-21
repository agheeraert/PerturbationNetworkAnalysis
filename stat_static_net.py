from CreateNetwork import CreateNetwork

L_cutoffs = list(range(3, 10))
for i in L_cutoffs:
    CreateNetwork(cutoff=i).save_avg(folder='/home/agheerae/Python/PerturbationNetworkAnalysis/data/sim1/apo/', output='/home/agheerae/results/static/apo_'+str(i))
    CreateNetwork(cutoff=i).save_avg(folder='/home/agheerae/Python/PerturbationNetworkAnalysis/data/sim1/prfar/', output='/home/agheerae/results/static/prfar_'+str(i))
