from commons import *


## point 1.1 func ##

def save_matrices(mat, n_channels=64, conn_meth='pdc', freq=8, run="R01"):
    """
    Save adjacency matrices obtained from DTF and PDC connectivity analysis to file
        - dtf_mat    : connectivity matrix obtained with DTF measure;
        - pdc_mat    : connectivity matrix obtained with PDC measure;
        - n_channels : number of channels in the data, i.e. the resulting matrices 
                dim (n_channels,n_channels);
        - freq       : frequecy value related to the matrix data;
        - run        : the related run of the experiment, one of {'R01','R02'}.
    """
    path = "data/{}_{}_{}hz_auto.txt".format(conn_meth, run,freq)

    print("\nSaving to {}".format(path))
    f = open(path, "w")
    for i in range(n_channels):
        for j in range(n_channels):
            if j == 63:
                f.write(str(mat[i][j]) + "\n")
            else:
                f.write(str(mat[i][j]) + " ")
    f.close()
    return


def load_matrix(conn_method="pdc", freq=10, run="R01", auto='auto', verbose=True):
    """
    Load the adjacency matrix from file
        - conn_method : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq        : the frequqncy value related to the matrix data;
        - run         : the related run of the experiment, one of {'R01','R02'}.
    """
    mat_file = "data/{}_{}_{}hz_{}.txt".format(conn_method, run, freq, auto) 
    mat_list = []
    if verbose:
        print("Loading matrix from '{}' ...".format(mat_file))

    with open(mat_file, "r") as f:
        for row in f.readlines():
            mat_list.append(row.strip().split(" "))
        f.close()

    conn_mat = np.array(mat_list, dtype=np.float32)
    np.fill_diagonal(conn_mat, 0.)
    return conn_mat


def compute_adjacency(conn_mat, threshold=0.1226):    
    """
    Compute binary adjacency matrix from the given connectivity matrix.
        - conn_mat  : the connectivity matrix to be binarified;
        - threshold : each cell of the resulting matrix will be considered 1 if its value
                is greater or equal than 'threshold', 0 otherwise. 
    """
    adj_mat = np.zeros(shape=conn_mat.shape)

    for i in range(conn_mat.shape[0]):
        for j in range(conn_mat.shape[1]):
            if conn_mat[i][j] >= threshold and i != j:
                adj_mat[i][j] = 1

    return adj_mat


def compute_mean_adjacency(conn_method='pdc', run='R01'):
    aux_mat = []
    for i in range(8, 14):
        conn_mat = load_matrix(conn_method=conn_method, freq=i, run=run, verbose=False)
        aux_mat.append(conn_mat)
    
    mean_mat = np.mean(aux_mat, axis=0)
    plt.matshow(mean_mat)
    plt.title("Mean {} adjacency matrix of run {} in the alpha-band".format(conn_method, run))
    plt.colorbar()
    plt.show()

    aux_mat = []
    for i in range(8, 14):
        conn_mat = load_matrix(conn_method=conn_method, freq=i, run=run,  verbose=False)
        aux_mat.append((conn_mat - mean_mat)**2)
    var_mat = np.mean(aux_mat, axis=0)
    plt.matshow(var_mat)
    plt.title("Variance of the {} adjacency matrix of run {} in the alpha-band".format(conn_method, run))
    plt.colorbar()
    plt.show()
    return var_mat


def check_mean_var_EEG_rithm(plot=True):
    if plot:
        for meth in ['pdc', 'dtf']:
            for run in ['R01', 'R02']:
                compute_mean_adjacency(conn_method=meth, run=run)


def map_index_to_channels():
    """
    Maps channel coordinates to indeces in the file
    """
    # relabel nodes cause there's no labels in the original EDF file
    # (well, there are, but they are not automatically loaded, so...)
    with open("data/channel_locations.txt") as f:
        pos = {}
        for line in f:
            # yes, there are 8 spaces in the file.
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                pos.update({int(l[0])-1: str(l[1])})
    return pos


def load_conn_graph(conn="pdc", freq=10, run="R01", auto="auto", threshold=0.1226):
    """
    Load the connectivity graph from the related connectivity matrix.
        - conn : the method used to compute the connectivity matrix, one of {'dtf','pdc'};
        - freq : the frequqncy value related to the matrix data;
        - run  : the related run of the experiment, one of {'R01','R02'}.
    """
    print("\nInitializing graph from {}-{}-{}hz matrix ...".format(conn,run,freq))
    adj_mat = compute_adjacency(load_matrix(conn_method=conn,freq=freq,run=run,auto=auto), threshold=threshold)

    G = nx.from_numpy_array(adj_mat,create_using=nx.DiGraph)
    mapping = map_index_to_channels()
    return nx.relabel_nodes(G, mapping)


def already_computed():
    """
    If we have already computed the adj mats for point 1.1 do not compute them again\n
    They are computed for both runs for all freq in the alpha band.
    """
    for run in ['R01', 'R02']:
        for freq in range(8, 14):
            dtf_path = "data/dtf_{}_{}hz_auto.txt".format(run, freq)
            pdc_path = "data/pdc_{}_{}hz_auto.txt".format(run, freq)
            if not os.path.isfile(dtf_path) or not os.path.isfile(pdc_path):
                return False
    return True


def print_adj(conn_method='pdc', freq=10, run='R01', threshold=None, auto='auto'):
    """
    Prints adjacency matrix    
    """

    mat = load_matrix(conn_method=conn_method, freq=freq, run=run, auto=auto)
    plt.matshow(mat)
    plt.title("{} adjacency matrix of run {} @{}Hz".format(conn_method, run, freq))
    plt.colorbar()
    plt.show()

    if threshold is not None:
        mat = compute_adjacency(mat, threshold=threshold)
        density = 100*np.sum(mat)/4032
        print("Density = {:.02f}%".format(density))
        plt.matshow(mat)
        plt.title("{} binary adjacency matrix of run {} @{}Hz with density = {:.02f}%".format(conn_method, run, freq, density))
        plt.show()

 
def p1_1(file_name=None, freq=10, run='R01', point='1'):
    
    #### Load EEG data from edf file
    if point == '4':        ### <<<<<<<<<<<<<<<<<<
        if not os.path.isfile(file_name + '_dropped.edf'):
            small_group = ['Fp1.', 'Fp2.', 'F7..', 'F3..', 'Fz..',
                       'F4..', 'F8..', 'T7..', 'C3..', 'Cz..',
                       'C4..', 'T8..', 'P7..', 'P3..', 'Pz..',
                       'P4..', 'P8..', 'O1..', 'O2..']
            pyedflib.highlevel.drop_channels(file_name+".edf", to_keep=small_group)
            # NOTE >>>> first time it computes the edf files it crashes, but the file is there... <<<<<<<<<<<<<<<<<<
        file_name = file_name + '_dropped'

    print("\n[1.{}] >> Analyzing file {}".format(point, file_name))
    f = pyedflib.EdfReader(file_name + ".edf")
    n = f.signals_in_file
    signal_labels = f.getSignalLabels()
    sigbufs = np.zeros((n, f.getNSamples()[0]))
    for i in np.arange(n):
        sigbufs[i, :] = f.readSignal(i)

    print("[1.{}] >> Loaded matrix with shape {}".format(point, sigbufs.shape))
    f.close()

    data = cp.Data(sigbufs, fs=160., chan_names=signal_labels, data_info=file_name)
    ''' not required
    with warnings.catch_warnings():             # stupid warning about the plot...
        warnings.simplefilter("ignore")
        fxn()
        data.plot_data(trial=3)
    '''

    print("[1.{}] >> Optimizing p...".format(point))
    mv = cp.Mvar
    best_p, crit = mv.order_akaike(sigbufs, p_max=30, method='yw')    
    if PLOTS:
        plt.plot(1+np.arange(len(crit)), crit, 'g')
        plt.title("Model order estimation")
        plt.xlabel("order(p)")
        plt.ylabel("AIC(p)")
        plt.grid()
        plt.show()

    print("[1.{}] >> Best p = {}".format(point, best_p))
    data.fit_mvar(p=best_p, method='yw')
    
    if point == '4':
        return data

    #### Compute connectivity matrices with DTF and PDC measures
    #if not already_computed() or COMPUTE_MATS:
    
    dtf_path = "data/dtf_{}_{}hz_auto.txt".format(run, freq)
    if not os.path.isfile(dtf_path):
        # investigate connectivity using DTF
        dtf_values = data.conn('dtf',resolution=80)
        dtf_significance = data.significance(Nrep=100, alpha=0.05)
        print("[1.1] >> dtf_shape:",dtf_values.shape)
        print("\n[1.1] >> DTF sign:", dtf_significance)
        #if PLOTS:
        #    data.plot_conn("DTF measure")
        save_matrices(mat = dtf_values[freq], n_channels=64, conn_meth='dtf', freq=freq, run=run)

    pdc_path = "data/pdc_{}_{}hz_auto.txt".format(run, freq)
    if not os.path.isfile(pdc_path):
        # investigate connectivity using PDC
        pdc_values = data.conn('pdc',resolution=80)
        pdc_significance = data.significance(Nrep=100, alpha=0.05)
        print("[1.1] >> pdc_shape:", pdc_values.shape)
        print("\n[1.1] >> PDC sign:", pdc_significance)
        #if PLOTS :
        #    data.plot_conn("PDC measure")
        save_matrices(mat=pdc_values[freq], n_channels=64, conn_meth='pdc', freq=freq, run=run)

    if PLOTS and point != '6':
        # NOTE: the threshold values here will yield ~20% density in the network.
        if run == 'R01':
            print_adj(conn_method='pdc', freq=freq, run=run, threshold=0.1226)
            print_adj(conn_method='dtf', freq=freq, run=run, threshold=0.1378)
        elif run == 'R02':
            print_adj(conn_method='pdc', freq=freq, run=run, threshold=0.1167)
            print_adj(conn_method='dtf', freq=freq, run=run, threshold=0.1322)


def p1_4(R='R01'):
    data = p1_1(point='4')
    pdc_path = "data/pdc_{}_10hz_19_channels.txt".format(R)

    if not os.path.isfile(pdc_path):
        pdc_values = data.conn('pdc', resolution=80)
        data.significance(Nrep=100, alpha=0.05)         # returns pdc_significance but what to do with it? However, it does side effect
        pdc_mat = pdc_values[10]
        f_pdc = open(pdc_path, "w")
        for i in range(19):
            for j in range(19):
                if j == 18:
                    f_pdc.write(str(pdc_mat[i][j]) + "\n")
                else:
                    f_pdc.write(str(pdc_mat[i][j]) + " ")
        f_pdc.close()
    
    if PLOTS:
        print("[1.4] >> Plotting connectivity...")
        print_adj(conn_method='pdc', freq=10, run='R01', threshold=None, auto='19_channels')
        # data.plot_conn("PDC measure")     # is it even useful?


def load_channel_coordinates(label=True, map_ch=False):
    """
    Loads channels coordinates in a disctionary and returns it
    """
    with open("data/channel_locations.txt") as f:
        pos = {}
        for line in f:
            # yes, there are 8 spaces in the file.
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                if label:
                    pos.update({str(l[1]): [float(l[2]), float(l[3])]})
                elif map_ch:
                    pos.update({int(l[0]) : str(l[1])})
                else:
                    pos.update({int(l[0])-1 : [float(l[2]), float(l[3])]})
    return pos


def p1_5(G, point='1.5', communities=None, nodelist=None, edgelist=None):
    """
    Prints a topological representation of the networks
    Node colors depend on their degree
    """
    if nodelist is None:
        nodelist = G.nodes()
    if edgelist is None:
        edgelist = G.edges()
    
    pos = load_channel_coordinates()

    def p1_5_helper(G, pos, degree, node_color, point='1.5'):
        """
        Helper function to now write two times the same plt stuff
        """
        cmap = 'viridis' if point == '1.5' else 'plasma'
        vmin = min(node_color)
        vmax = max(node_color)

        nc = nx.draw_networkx_nodes(G, pos=pos, vmin=vmin, vmax=vmax, edgecolors='black', node_size=700, node_color=node_color, cmap=cmap)

        if point == '1.5':
            _  = nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='black', arrows=True, node_size=700)

        elif point == '3.2':
            edge_color = [G[u][v]['color'] for u, v in G.edges()]
            count = 0
            for i in edge_color:
                if i == 'r':
                    count += 1
            print("[3.2] >> FOUND {} on {} EDGES INVOLVED IN MOTIF 1".format(count, len(edge_color)))
            print("[3.2] >> {:.2f}% of edges are involved in motif 1".format(100*count/len(edge_color)))
            _  = nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color=edge_color, arrows=True, node_size=700)

        _  = nx.draw_networkx_labels(G, pos)
        if point == '1.5':
            plt.title("Topological representation of the network - {} degree".format(degree))
        elif point == '3.2':
            plt.title("Topological representation of the network's edges involved in motif 1 - {} degree".format(degree))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm._A = []
        plt.colorbar(sm)
        plt.show()


    def p2_5_helper(G, pos):
        node_color_in  = []
        node_color_out = []
        node_color_sum = []
        for node in G.nodes():
            node_color_in.append(G.in_degree(node))
            node_color_out.append(G.out_degree(node))
            node_color_sum.append(G.out_degree(node)-G.in_degree(node))
        p1_5_helper(G, pos, degree='in',  node_color=node_color_in,  point='1.5')
        p1_5_helper(G, pos, degree='out', node_color=node_color_out, point='1.5')
        p1_5_helper(G, pos, degree='sum', node_color=node_color_sum, point='1.5')


    def p3_2_helper(G, pos):
        node_color_in  = []
        node_color_out = []
        node_color_sum = []
        for node in G.nodes():
            node_color_in.append( G.in_degree(node))
            node_color_out.append(G.out_degree(node))
            node_color_sum.append(G.out_degree(node) - G.in_degree(node))
        p1_5_helper(G, pos, degree='in' , node_color=node_color_in,  point='3.2')
        p1_5_helper(G, pos, degree='out', node_color=node_color_out, point='3.2')
        p1_5_helper(G, pos, degree='sum', node_color=node_color_sum, point='3.2')


    def p4_2_helper(G, pos, communities):
        
        colorbar_labels = []
        S = set(communities.values())
        for i in S:
            colorbar_labels.append('Community {}'.format(i+1))

        cmap = 'Spectral'
        vmin = min(communities.values())
        vmax = max(communities.values())
        _  = nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='black', arrows=True, node_size=700)
        nc = nx.draw_networkx_nodes(G, pos=pos, vmin=vmin, vmax=vmax, edgecolors='black', node_size=700,
                                    nodelist=communities.keys(), node_color=list(communities.values()), cmap=cmap)
        _  = nx.draw_networkx_labels(G, pos)
        cbar = plt.colorbar(nc, ticks=np.arange(len(communities)), spacing='proportional')
        cbar.ax.set_yticklabels(colorbar_labels)
        plt.title("Topological representation of the network's communities found with Louvain")
        plt.show()


    def p4_3_helper(G, pos, communities):

        colorbar_labels = []
        S = set(communities)
        for i in S:
            colorbar_labels.append('Community {}'.format(i+1))

        cmap = 'Spectral'
        vmin = min(communities)
        vmax = max(communities)
        _  = nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color='black', arrows=True, node_size=700)
        nc = nx.draw_networkx_nodes(G, pos=pos, vmin=vmin, vmax=vmax, edgecolors='black', node_size=700, node_color=communities, cmap=cmap)
        _  = nx.draw_networkx_labels(G, pos)
        cbar = plt.colorbar(nc, ticks=np.arange(len(communities)), spacing='proportional')
        cbar.ax.set_yticklabels(colorbar_labels)
        plt.title("Topological representation of the network's communities found with Infomap")
        plt.show()

    
    if point == '1.5':
        node_color_in  = []
        node_color_out = []
        node_color_sum = []
        for node in G.nodes():
            node_color_in.append(G.in_degree(node))
            node_color_out.append(G.out_degree(node))
            node_color_sum.append(G.out_degree(node) + G.in_degree(node))
        p1_5_helper(G, pos, degree='in' , node_color=node_color_in)
        p1_5_helper(G, pos, degree='out', node_color=node_color_out)
        p1_5_helper(G, pos, degree='sum', node_color=node_color_sum)
    elif point == '2.5':
        p2_5_helper(G, pos)
    elif point == '3.2':
        p3_2_helper(G, pos)
    elif point == '4.2':
        p4_2_helper(G, pos, communities)
    elif point == '4.3':
        p4_3_helper(G, pos, communities)


def find_threshold(mat, target, start, point='3'):
    print('[1.{}] >> Optimizing thresold to reach {}% density..'.format(point, target))
    old = 100
    best_t = -1
    eps = 0.05
    # it's faster if we decrease our threshold
    vv = np.arange(start, 0.06, -0.00005)
    for t in vv:
        new_mat = compute_adjacency(mat, threshold=t)
        curr_d = 100*np.sum(new_mat)/4032
        if abs(curr_d - target) <= eps:
            if abs(curr_d - target) <= old:
                best_t = t
            else:
                break
        old = abs(curr_d - target)
    return best_t


def p1_3(conn_method, freq, run, point='3'):
    print("------------------------------------------------------")
    print("[1.{}] >> Optimizing ".format(point), conn_method, "matrix, run", run)
    mat = load_matrix(conn_method=conn_method, freq=freq, run=run, verbose=False)
    
    # threshold for minimum freq with pdc
    best_t = 0.27   
    for target in list([1, 5, 10, 20, 30, 50]):
        best_t = find_threshold(mat, target, best_t)
        new_mat = compute_adjacency(mat, threshold=best_t)
        density = 100*np.sum(new_mat)/4032
        print("Density ~ {:.02f}% with threshold = {}\n".format(density, best_t))
        plt.matshow(new_mat)
        plt.title("{} binary adjacency matrix of run {} @{}Hz with density = {:.02f}%".format(conn_method, run, freq, density))
        plt.show()


def p1_6(file_name="data/S003R01_fixed", freq=25, run='R01'):
    p1_1(file_name=file_name, freq=freq, run=run, point='6')
    mat = load_matrix(conn_method='pdc', freq=freq, run=run, verbose=True)
    threshold = find_threshold(mat, target=20, start=0.4, point='6')
    print("[1.6] >> best threshold = {}\n".format(threshold))
    print_adj(conn_method='pdc', freq=freq, run=run, threshold=threshold)


if __name__ == "__main__":
    '''
    mat = load_matrix(conn_method='pdc', freq=10, run='R01')#, threshold=THRES_PDC_10HZ_R01_20percent)
    #print(mat)
    pos = {}
    with open("data/channel_locations.txt") as f:
        for line in f:
            # yes, there are 8 spaces in the file.
            l = line.split(sep='        ')
            if l[0] != '\ufeff#':
                pos.update({str(l[1]): int(l[0])-1})
    #print(pos)
    left_indices   = []
    right_indices  = []
    center_indices = []
    for k in pos.keys():
        if '1' in k or '3' in k or '5' in k or '7' in k or '9' in k:
            left_indices.append(pos[k])
        elif '0' in k or '2' in k or '4' in k or '6'in k  or '8' in k:
            right_indices.append(pos[k])
        elif 'z' in k:
            center_indices.append(pos[k]) 
    print(left_indices,   '\n')
    print(right_indices,  '\n')
    print(center_indices, '\n')
    print(len(left_indices)+len(right_indices)+len(center_indices)) # if 64 ok

    shuffled = np.zeros(shape=mat.shape)

    for idx in left_indices:
        for i in range(len(left_indices)):
            for j in range(len(left_indices)):
                shuffled[i][j] = mat[idx][idx]  #wrong
    plt.matshow(mat)
    plt.colorbar()
    plt.show()
    plt.matshow(shuffled)
    plt.colorbar()
    plt.show()
    '''
