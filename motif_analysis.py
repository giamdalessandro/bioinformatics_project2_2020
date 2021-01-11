import networkx
from connectivity_graph import load_conn_graph


triad_cfg = {
	'003' : 'Null',
	'012' : 'Single-edge',
	'021C': 'Pass-along',
	'021D': 'Double-dominant',
	'021U': 'Double-subordinate',
	'030C': 'Cycle',
	'030T': 'Transitive'
}


net_G = load_conn_graph()
census = nx.triadic_census(net_G)

#sp = triadSignificanceProfile(net_G, triad_cfg)

f_census = {}
f_census['group-size'] = [N_IND]
f_census['flee-dist'] = [params['female.FleeDist']]
f_census['aggr-intensity'] = [('mild' if params['Rating.Dom.female.Intensity'] == 0.1 else 'fierce')]
f_census['steepness'] = round(steep,4)

print('\nNetwork Triadic Census:')
for k,v in sorted(census.items()):
	if k in triad_cfg:
		f_census[triad_cfg[k]] = [v]
		print('  ' + triad_cfg[k] + ': ' + str(v))