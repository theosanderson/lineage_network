import streamlit as st
import json
import requests
import networkx as nx
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from itertools import combinations
from io import BytesIO
import numpy as np
from scipy.optimize import minimize


def get_graph_layout(dists):
    


    # If any pairwise distances are 0 raise an error and say what they are
    zero_distances = [pair for pair, distance in dists.items() if distance == 0]
    if zero_distances:
        st.error(f"The following pairs have a distance of 0: {zero_distances}. A graph cannot be plotted. Please remove one of the pair.")

    G = nx.Graph()


    for pair, distance in dists.items():
        G.add_edge(pair[0], pair[1], weight=distance)
    
    # create random initial positions for the nodes
    
    initial_pos = nx.random_layout(G)

    pos = get_custom_layout2(G, iterations=100)
    return pos, G



def get_custom_layout2(G, iterations=50):
    """
    Calculate a custom layout for the nodes in a graph G based on stress minimization.
    """

    def stress(flat_positions):
        positions = flat_positions.reshape((len(nodes), 2))
        stress_value = 0
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                distance = np.linalg.norm(positions[i] - positions[j])
                graph_distance = distances.get((nodes[i], nodes[j]), distances.get((nodes[j], nodes[i]), 0))
                stress_value += (graph_distance - distance) ** 2
        return stress_value

    nodes = list(G.nodes())
    distances = nx.get_edge_attributes(G, 'weight')

    # Initial random positions flattened into a one-dimensional array
    initial_positions = np.random.rand(len(nodes) * 2)

    # Minimize the stress function
    result = minimize(stress, initial_positions, method='CG', options={'maxiter': iterations})

    # Reshape the result and create a position dictionary
    optimized_positions = result.x.reshape((len(nodes), 2))
    pos = {nodes[i]: optimized_positions[i] for i in range(len(nodes))}

    return pos

def optimize_layout(dists, num_iterations=100):
    best_correlation = -1
    best_layout = None
    best_graph = None

    for _ in range(num_iterations):
        pos, G = get_graph_layout(dists)

        # Calculate the correlation between graph distances and layout distances
        correlation , _ , _= calculate_correlation(G, pos)
        
        if correlation > best_correlation:
            print(f"New best correlation: {correlation}")
            best_correlation = correlation
            best_layout = pos
            best_graph = G

    return best_layout, best_graph

def calculate_correlation(G, pos):
    graph_distances = []
    layout_distances = []

    for edge in G.edges():
        node1, node2 = edge
        graph_distance = G[node1][node2]['weight']
        layout_distance = np.linalg.norm(np.array(pos[node1]) - np.array(pos[node2]))

        graph_distances.append(graph_distance)
        layout_distances.append(layout_distance)

    correlation = np.corrcoef(graph_distances, layout_distances)[0, 1]
    return correlation, graph_distances, layout_distances

# Function to download and load JSON data from URL
def download_data(url):
    response = requests.get(url)
    if response.status_code == 200:
        return json.loads(response.text)
    else:
        return None

# Function to collect mutation paths
def collect_mutation_paths(node, current_path, paths, mutation_key):
    node_path = current_path.copy()
    if 'branch_attrs' in node and 'mutations' in node['branch_attrs']:
        node_mutations = node['branch_attrs']['mutations'].get(mutation_key, [])
    else:
        node_mutations = []
        
        
    # remove any sites that are already there
    new_node_path = []
    node_poses = [item[1:-1] for item in node_mutations]
    for item in node_path:
      
        item_number = item[1:-1]
        
        if item_number not in node_poses:

            new_node_path.append(item)
  
    node_path = new_node_path
    node_path.extend(node_mutations)
    paths[node['name']] = node_path

    if 'children' in node:
        for child in node['children']:
            collect_mutation_paths(child, node_path, paths, mutation_key)

# Function to calculate edit distance
def calculate_edit_distance(str1, str2, mutation_paths,debug_print):
    list1 = mutation_paths[str1]
    list2 = mutation_paths[str2]

    list1 = [item[1:] for item in list1]
    list2 = [item[1:] for item in list2]

    set1 = set(list1)
    set2 = set(list2)

    unique_set1 = set1 - set2
    unique_set2 = set2 - set1

    processed_set = set([item[:-1] for item in unique_set1.union(unique_set2)])
    items = sorted([int(item) for item in processed_set])
    result = len(processed_set)
    if(debug_print):
        st.write("Sites with differences between",str1," and ",str2, ":",result, ",".join([str(item) for item in items]))

    return result

# Function to calculate all pairwise distances
def calculate_all_pairwise_distances(keys, mutation_paths,debug_print):
    pairwise_distances = {}
    for pair in combinations(keys, 2):
        distance = calculate_edit_distance(pair[0], pair[1], mutation_paths, debug_print = debug_print)
        pairwise_distances[pair] = distance

    return pairwise_distances

# Function to draw network graph
def draw_graph(pos, G, keys, show_edge_labels, mutation_key, file_format='png'):

    plt.figure(figsize=(12, 12))
    plt.tight_layout()
    #plt.axis('off')
    plt.gca().set_aspect('equal')
    plt.gca().set_axis_off()
    plt.margins(0, 0)

    # tight layout
  

    # Draw nodes

    for key in keys:
        node_labels = {key: key}
        nx.draw_networkx_nodes(G, pos, nodelist=[key], node_size=500, node_color='skyblue', label=key)
        nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10)
    
    
    if show_edge_labels:
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

    # Draw edges
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)




    buffer = BytesIO()
    
    plt.savefig(buffer, format=file_format)
    buffer.seek(0)
    return buffer

# Streamlit app layout
st.title('SARS CoV-2 Variant Mutation Network')


data_url = "https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2"
json_data = download_data(data_url)
mutation_type = st.radio("Select Mutation Type", ('Nucleotide', 'Spike Protein'))
mutation_key = 'nuc' if mutation_type == 'Nucleotide' else 'S'

if json_data:
    tree_data = json_data['tree']
    mutation_paths = {'ancestral':[]}
    collect_mutation_paths(tree_data, [], mutation_paths, mutation_key)

    # Selectable keys for visualization
    keys = st.multiselect('Select Variants', list(mutation_paths.keys()), default=['ancestral', 'B.1.1.7','BA.1','B.1.617.2','BA.2','BA.5', 'XBB.1.5', 'BA.2.86']
                          
                          )

    # Toggle for mutation type
    
    

    # Toggle for edge labels
    show_edge_labels = st.checkbox("Show edge labels", value=True)

    enable_site_display = st.checkbox("List sites with variation", value=False)

    
    if keys:

        dists = calculate_all_pairwise_distances(keys, mutation_paths, debug_print = enable_site_display)

        pos, G=  optimize_layout(dists, num_iterations=10)
       
        svg_buffer  =draw_graph(pos, G, keys, show_edge_labels, mutation_key, file_format='svg')

        buffer = draw_graph(pos, G, keys, show_edge_labels, mutation_key, file_format='png')
        st.image(buffer)

        # Plot correlation
        correlation, graph_distances, layout_distances = calculate_correlation(G, pos)

        # scatter
        fig, ax = plt.subplots(figsize=(5,4), dpi=60)
        # lowdpi

       

        ax.scatter(graph_distances, layout_distances)
        ax.set_xlabel('True distance')
        ax.set_ylabel('Plotted distance')
        ax.set_title(f'Correlation of distances: {correlation:.2f}')
        # save to buffer 
        buffer = BytesIO()
        fig.savefig(buffer, format='png')
        buffer.seek(0)
        st.image(buffer)




    
        st.download_button(
            label="Download graph as SVG",
            data=svg_buffer,
            file_name=f"graph.svg",
            mime="image/svg+xml"
        )

        st.download_button(
            label="Download graph as PNG",
            data=buffer,
            file_name=f"graph.png",
            mime="image/png"
        )

        # give credit to Nextstrain
        st.markdown("Data from [Nextclade](https://nextstrain.org/staging/nextclade/sars-cov-2)")

else:
    st.error("Failed to load data from the URL.")
