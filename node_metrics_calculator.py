import graph_constructor as gc
import os, csv
import pandas as pd 

create_data = False # Set to true if first time running and there is no PSI_chls.csv data.

# Load dataframe to initiate Nodes
ps1_df=pd.read_csv('chls_1JB0.csv')
ps1_df.drop(columns=['Index'], inplace=True)
ps1_df['Energy (cm^-1)'] = float('nan')

# Load site energy directory
ps1_site_energy_dir = os.getcwd() + '\site energies\PSI'

for i in range(len(os.listdir(ps1_site_energy_dir))):
    df_energy = pd.read_csv(ps1_site_energy_dir + '\\' + os.listdir(ps1_site_energy_dir)[i])
    # merged_df = pd.merge(merged_df, df_energy, on='ID', how='left') #Merge two data frames 
    # df['Energy (cm^-1)']=df['Energy (cm^-1)'].fillna((df_energy.set_index)('ID')['Energy (cm^-1)'])
    ps1_df = ps1_df.set_index('ID').combine_first(df_energy.set_index('ID')).reset_index() 

if create_data:
    g = gc.Graph()
    N = 96
    #Initialize all nodes
    for i in range(N): 
        g.addNode(gc.Node(ps1_df, i))

    for i in range(N):
        for j in range(N):
            if i != j:
                node_i = g.nodes[i]
                node_j = g.nodes[j]
                edge_add = gc.Edge(node_i, node_j)
                edge_add.setWeight(gc.transfer_rate(node_i, node_j))
                g.addEdge(edge_add)

    g.init_Weights()
    g.init_Costs()
    g.P700lookup()
    fields = ['ID', 'inStrength', 'outStrength', 'farness', 'BC', 'transferToP700', 'euclidean_distance_to_P700']

    chls = [str(i) for i in g.nodes] 
    inStrengths = [g.nodeInStrength(g.nodes[i]) for i in range(len(chls))]
    outStrengths = [g.nodeOutStrength(g.nodes[i]) for i in range(len(chls))]
    transferToP700 = [g.transferToP700(g.nodes[i]) for i in range(len(chls))]
    p700_euclidean_distances =[g.euclidean_node_P700distance(g.nodes[i]) for i in range(len(chls))]
    distances = [min( g.distance(g.nodes[i], g.p700a), g.distance(g.nodes[i], g.p700b) ) for i in range(len(chls))]

    bc_dict = g.betweenessCentrality()

    BCs = [bc_dict[g.nodes[i]] for i in range(len(chls))]

    filename = os.path.join(os.getcwd(), "data", "PSI_chls.csv")
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(filename, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        
        # Write the header row
        writer.writeheader()
        
        for i in range(len(chls)):
            row = {'ID': chls[i], 'inStrength': inStrengths[i], 'outStrength': outStrengths[i], 'farness': distances[i], 'transferToP700': transferToP700[i], 'BC': BCs[i], 
                'euclidean_distance_to_P700': p700_euclidean_distances[i]}
            writer.writerow(row)

