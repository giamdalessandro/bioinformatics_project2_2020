import numpy as np
import networkx as nx
from connectivity_graph import load_conn_graph



def mapTriadCodes(census, rand_census, triad_cfg):
	"""
	Maps nx.triadic_census() subgraph codes to explicit to triadic patterns names.
	"""
	real = {}
	for k,v in sorted(census.items()):
		if k in triad_cfg:
			real[triad_cfg[k]] = v

	random = {}
	for rc in rand_census:
		for k,v in sorted(rc.items()):
			if k in triad_cfg:
				if triad_cfg[k] not in random.keys():
					random[triad_cfg[k]] = []
				
				random[triad_cfg[k]].append(v)
	
	return (real, random)

def triadSignificanceProfile(G, triad_cfg):
    """
    Compute the significance profile of the patterns mapped in triad_cfg,
    inside directed graph G.
        - G          : directed graph representing the network; 
        - triads_cfg : dict mapping interesting triadic patterns codes, 
            as in nx.triadic_census(), with explicit names. 
            (e.g. triad_cfg = {'003' : 'Null', '012' : 'Single-edge'})
    """
    print('Computing Z-score...')
	census = nx.triadic_census(G)
    in_degree_sequence = [d for n, d in G.in_degree()]  # in degree sequence
    out_degree_sequence = [d for n, d in G.out_degree()]  # out degree sequence
    #print("In_degree sequence %s" % in_degree_sequence)
    #print("Out_degree sequence %s" % out_degree_sequence)

    random_nets_census = []
    for i in range(100):
        rand_G = nx.directed_configuration_model(in_degree_sequence, out_degree_sequence, create_using=nx.DiGraph, seed=i)
        random_nets_census.append(nx.triadic_census(rand_G))

    real_census, random_census = mapTriadCodes(census,random_nets_census,triad_cfg)
    #print(real_census)
    #print(random_census)

    z_score = []
    for p in real_census.keys():
        #print(p)
        N_real_p = real_census[p]
        N_rand_p = np.mean(random_census[p])
        std = np.std(random_census[p])

        z_p =  ((N_real_p - N_rand_p)/std if std != 0 else 0)
        z_score.append(z_p)

    sp = []
    for i in range(len(z_score)):
        z_norm = np.linalg.norm(z_score)
        norm_z_score = (z_score[i]/z_norm if z_norm != 0 else z_score[i])
        sp.append(round(norm_z_score,4))
    
    print('Significance profile is:')
    print(sp)
    return sp


triad_cfg = {
	'021D': 'type-1',
	'021C': 'Three-chain',
	'021U': 'type-3',
    '111D': 'type-4',
    '111U': 'type-5',
    '030T': 'Feed-forward',
	'030C': 'Feedback',
    '201' : 'type-8',
    '120D': 'type-9',
    '120U': 'type-10',
    '120C': 'type-11',
    '210' : 'type-12',
    '300' : 'type-13'
}


net_G = load_conn_graph()
census = nx.triadic_census(net_G)

sp = triadSignificanceProfile(net_G, triad_cfg)

f_census = {}
print('\ncomputing network triadic census...\n')
print('triadType  \tN')
print('--------------------')
for k,v in sorted(census.items()):
	if k in triad_cfg:
		f_census[triad_cfg[k]] = [v]
		print(triad_cfg[k] + ': \t' + str(v))