import streamlit as st
import json
import requests
import networkx as nx
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'
from itertools import combinations
from io import BytesIO

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
        node_path.extend(node_mutations)

    paths[node['name']] = node_path

    if 'children' in node:
        for child in node['children']:
            collect_mutation_paths(child, node_path, paths, mutation_key)

# Function to calculate edit distance
def calculate_edit_distance(str1, str2, mutation_paths):
    list1 = mutation_paths[str1]
    list2 = mutation_paths[str2]

    list1 = [item[1:] for item in list1]
    list2 = [item[1:] for item in list2]

    set1 = set(list1)
    set2 = set(list2)

    unique_set1 = set1 - set2
    unique_set2 = set2 - set1

    processed_set = set([item[:-1] for item in unique_set1.union(unique_set2)])

    return len(processed_set)

# Function to calculate all pairwise distances
def calculate_all_pairwise_distances(keys, mutation_paths):
    pairwise_distances = {}
    for pair in combinations(keys, 2):
        distance = calculate_edit_distance(pair[0], pair[1], mutation_paths)
        pairwise_distances[pair] = distance

    return pairwise_distances

# Function to draw network graph
def draw_graph(mutation_paths, keys, show_edge_labels, mutation_key, file_format='png'):
    dists = calculate_all_pairwise_distances(keys, mutation_paths)

    # If any pairwise distances are 0 raise an error and say what they are
    zero_distances = [pair for pair, distance in dists.items() if distance == 0]
    if zero_distances:
        st.error(f"The following pairs have a distance of 0: {zero_distances}. A graph cannot be plotted. Please remove one of the pair.")


    G = nx.Graph()
    for pair, distance in dists.items():
        G.add_edge(pair[0], pair[1], weight=distance)

    plt.figure(figsize=(12, 12))
    pos = nx.kamada_kawai_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size = 2000)
    
    if show_edge_labels:
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)


    buffer = BytesIO()
    
    plt.savefig(buffer, format=file_format)
    buffer.seek(0)
    return buffer

# Streamlit app layout
st.title('SARS CoV-2 Variant Mutation Network')

# Toggle for mutation type
mutation_type = st.radio("Select Mutation Type", ('Nucleotide', 'Spike Protein'), index=1)
mutation_key = 'nuc' if mutation_type == 'Nucleotide' else 'S'

# Toggle for edge labels
show_edge_labels = st.checkbox("Show Edge Labels", value=True)

data_url = "https://nextstrain.org/charon/getDataset?prefix=staging/nextclade/sars-cov-2"
json_data = download_data(data_url)

if json_data:
    tree_data = json_data['tree']
    mutation_paths = {'ancestral':[]}
    collect_mutation_paths(tree_data, [], mutation_paths, mutation_key)

    # Selectable keys for visualization
    keys = st.multiselect('Select Variants', list(mutation_paths.keys()), default=['ancestral', 'B.1.1.7'])
    
    if keys:

        svg_buffer  =draw_graph(mutation_paths, keys, show_edge_labels, mutation_key, file_format='svg')

        buffer = draw_graph(mutation_paths, keys, show_edge_labels, mutation_key)
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
