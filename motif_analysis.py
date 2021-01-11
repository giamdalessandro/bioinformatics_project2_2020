import networkx as nx
from connectivity_graph import load_conn_graph


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

#sp = triadSignificanceProfile(net_G, triad_cfg)

f_census = {}

print('\nNetwork Triadic Census:')
for k,v in sorted(census.items()):
	if k in triad_cfg:
		f_census[triad_cfg[k]] = [v]
		print('  ' + triad_cfg[k] + ': ' + str(v))