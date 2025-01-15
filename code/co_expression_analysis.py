# load dependencies 
!pip install scipy statsmodels networkx seaborn
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# Import data
file_path = 'countlist_5xFAD.csv'  # Path to your CSV file
data = pd.read_csv(file_path, index_col=0)  # Ensure the first column has gene names

# Subset the data to the first 100 rows (genes) for testing
data = data.dropna(subset=["Gene_Name"]) # Drop Gene_Name column
data_subset = data.iloc[:100, :]  # Select the first 100 rows (genes)
data_subset = data_subset.set_index("Gene_Name") # Set Gene_Name as index

# Ensure the expression data is numeric
expression_data_numeric = data_subset.apply(pd.to_numeric, errors='coerce')

# Transpose the expression data to focus on gene co-expression (genes as columns)
expression_data_transposed = expression_data_numeric.T

# Function to calculate pairwise correlations
def calculate_correlation_matrix(expression_data):
    correlation_matrix = expression_data.corr(method="spearman")  # Spearman correlation
    return correlation_matrix

# Perform combined co-expression analysis and then display first few rows of results data frame
combined_corr_matrix = calculate_correlation_matrix(expression_data_transposed)
combined_corr_matrix.head()

# Create a Network from the Correlation Matrix
def create_network(corr_matrix, threshold=0.7):
    """
    Creates a co-expression network from a correlation matrix.
    Only edges with correlation >= threshold are added.
    """
    G = nx.Graph()
    for i, gene1 in enumerate(corr_matrix.index):
        for j, gene2 in enumerate(corr_matrix.columns):
            if i < j:  # Avoid duplicate pairs
                correlation = corr_matrix.iloc[i, j]
                if abs(correlation) >= threshold:  # Apply threshold
                    G.add_edge(gene1, gene2, weight=correlation)
    return G

# Define threshold for correlation
threshold = 0.7

# Create combined network for all samples and analyze the combined network
combined_network = create_network(combined_corr_matrix, threshold)
print(f"Combined Network: {combined_network.number_of_nodes()} nodes and {combined_network.number_of_edges()} edges")

# Visualize the network with all gene names as labels
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(combined_network, seed=42)
nx.draw(combined_network, pos, with_labels=False, node_size=20, edge_color="lightblue")
nx.draw_networkx_labels(combined_network, pos, labels={node: node for node in combined_network.nodes}, font_size=7)
plt.title("Gene Co-expression Network with All Gene Names")
plt.show()

# identify hub genes in network
degrees = dict(combined_network.degree())
sorted_genes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
print("Top 10 hub genes:")
print(sorted_genes[:10])

degree_centrality = nx.degree_centrality(combined_network)
sorted_centrality = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
print("Top 10 genes by degree centrality:")
print(sorted_centrality[:10])

clustering = nx.clustering(combined_network)
sorted_clustering = sorted(clustering.items(), key=lambda x: x[1], reverse=True)
print("Top 10 genes by clustering coefficient:")
print(sorted_clustering[:10])

betweenness = nx.betweenness_centrality(combined_network)
sorted_betweenness = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)
print("Top 10 genes by betweenness centrality:")
print(sorted_betweenness[:10])

eigenvector = nx.eigenvector_centrality(combined_network)
sorted_eigenvector = sorted(eigenvector.items(), key=lambda x: x[1], reverse=True)
print("Top 10 genes by eigenvector centrality:")
print(sorted_eigenvector[:10])

# Select top 10 genes for each metric
top_degree_centrality_genes = [gene for gene, _ in sorted_centrality[:10]]
top_clustering_genes = [gene for gene, _ in sorted_clustering[:10]]
top_betweenness_genes = [gene for gene, _ in sorted_betweenness[:10]]
top_eigenvector_genes = [gene for gene, _ in sorted_eigenvector[:10]]

# Create subgraphs
subgraph_degree_centrality = combined_network.subgraph(top_degree_centrality_genes)
subgraph_clustering = combined_network.subgraph(top_clustering_genes)
subgraph_betweenness = combined_network.subgraph(top_betweenness_genes)
subgraph_eigenvector = combined_network.subgraph(top_eigenvector_genes)

# Plot subgraphs
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# Degree Centrality Subgraph
nx.draw(subgraph_degree_centrality,ax=axes[0, 0],with_labels=True,node_size=500,node_color="skyblue",font_size=10,edge_color="gray",)
axes[0, 0].set_title("Top 10 by Degree Centrality")

# Clustering Coefficient Subgraph
nx.draw(subgraph_clustering,ax=axes[0, 1],with_labels=True,node_size=500,node_color="lightgreen",font_size=10,edge_color="gray",)
axes[0, 1].set_title("Top 10 by Clustering Coefficient")

# Betweenness Centrality Subgraph
nx.draw(subgraph_betweenness,ax=axes[1, 0],with_labels=True,node_size=500,node_color="salmon",font_size=10,edge_color="gray",)
axes[1, 0].set_title("Top 10 by Betweenness Centrality")

# Eigenvector Centrality Subgraph
nx.draw(subgraph_eigenvector, ax=axes[1, 1],with_labels=True,node_size=500,node_color="plum",font_size=10,edge_color="gray",)
axes[1, 1].set_title("Top 10 by Eigenvector Centrality")

plt.tight_layout()
plt.show()

# Subset the data to exclude rows with NaN values in the "Gene_Name" column
data = data.dropna(subset=["Gene_Name"])

# Subset data for 5xFAD (mutant) and BL6 (WT) groups separately
five_x_fad_data = data.filter(regex="5xFAD")  # Select columns with "5xFAD" in their names
wt_data = data.filter(regex="BL6")  # Select columns with "BL6" in their names

# Ensure the expression data is numeric
five_x_fad_data_numeric = five_x_fad_data.apply(pd.to_numeric, errors='coerce')
wt_data_numeric = wt_data.apply(pd.to_numeric, errors='coerce')

# Set Gene_Name as the index for easier gene identification
five_x_fad_data_numeric = five_x_fad_data_numeric.set_index(data["Gene_Name"])
wt_data_numeric = wt_data_numeric.set_index(data["Gene_Name"])

# Select the first 100 rows (genes) for each group
five_x_fad_data_numeric = five_x_fad_data_numeric.iloc[:100, :]
wt_data_numeric = wt_data_numeric.iloc[:100, :]

# Transpose the expression data to focus on gene co-expression (genes as rows)
five_x_fad_transposed = five_x_fad_data_numeric.T  # Transpose for co-expression (samples as rows)
wt_transposed = wt_data_numeric.T  # Transpose for WT group

# Function to calculate pairwise correlations
def calculate_correlation_matrix(expression_data):

    correlation_matrix = expression_data.corr(method="spearman")  # Spearman correlation
    return correlation_matrix

# Calculate the correlation matrices separately for 5xFAD and WT groups
five_x_fad_corr_matrix = calculate_correlation_matrix(five_x_fad_transposed)
wt_corr_matrix = calculate_correlation_matrix(wt_transposed)

# Create networks for both treatment and control groups
def create_network(corr_matrix, threshold=0.7):
    G = nx.Graph()
    for i, gene1 in enumerate(corr_matrix.index):
        for j, gene2 in enumerate(corr_matrix.columns):
            if i < j:  # Avoid duplicate pairs
                correlation = corr_matrix.iloc[i, j]
                if abs(correlation) >= threshold:  # Apply threshold
                    G.add_edge(gene1, gene2, weight=correlation)
    return G

# Create networks for 5xFAD and WT
five_x_fad_network = create_network(five_x_fad_corr_matrix)
wt_network = create_network(wt_corr_matrix)

# Plot side-by-side subgraphs
plt.figure(figsize=(14, 6))
plt.subplot(1, 2, 1)  # 1 row, 2 columns, first subplot
pos_five_x_fad = nx.spring_layout(five_x_fad_network, seed=42)
nx.draw(five_x_fad_network, pos_five_x_fad, with_labels=False, node_size=20, edge_color="lightblue")
plt.title("5xFAD Gene Co-expression Network")
plt.subplot(1, 2, 2)  # 1 row, 2 columns, second subplot
pos_wt = nx.spring_layout(wt_network, seed=42)
nx.draw(wt_network, pos_wt, with_labels=False, node_size=20, edge_color="lightgreen")
plt.title("WT Gene Co-expression Network")
plt.tight_layout()
plt.show()

# Identify and print top 10 hub genes for both networks
five_x_fad_degrees = dict(five_x_fad_network.degree())
sorted_five_x_fad_genes = sorted(five_x_fad_degrees.items(), key=lambda x: x[1], reverse=True)
print("Top 10 hub genes in 5xFAD network:")
print(sorted_five_x_fad_genes[:10])
wt_degrees = dict(wt_network.degree())
sorted_wt_genes = sorted(wt_degrees.items(), key=lambda x: x[1], reverse=True)
print("\nTop 10 hub genes in WT network:")
print(sorted_wt_genes[:10])

density_five_x_fad = nx.density(five_x_fad_network)
density_wt = nx.density(wt_network)
print(f"Density of 5xFAD network: {density_five_x_fad}")
print(f"Density of WT network: {density_wt}")

avg_degree_five_x_fad = sum(dict(five_x_fad_network.degree()).values()) / len(five_x_fad_network.nodes)
avg_degree_wt = sum(dict(wt_network.degree()).values()) / len(wt_network.nodes)
print(f"Average degree in 5xFAD network: {avg_degree_five_x_fad}")
print(f"Average degree in WT network: {avg_degree_wt}")

# Set of hub genes for 5xFAD and WT networks
hub_genes_five_x_fad = set([gene for gene, degree in sorted_five_x_fad_genes[:10]])
hub_genes_wt = set([gene for gene, degree in sorted_wt_genes[:10]])

# Find common and unique hub genes
common_hub_genes = hub_genes_five_x_fad.intersection(hub_genes_wt)
unique_five_x_fad_hub_genes = hub_genes_five_x_fad - hub_genes_wt
unique_wt_hub_genes = hub_genes_wt - hub_genes_five_x_fad

print(f"Common hub genes: {common_hub_genes}")
print(f"Unique 5xFAD hub genes: {unique_five_x_fad_hub_genes}")
print(f"Unique WT hub genes: {unique_wt_hub_genes}")

# Extract edge weights from 5xFAD and WT networks
edge_weights_five_x_fad = [data['weight'] for _, _, data in five_x_fad_network.edges(data=True)]
edge_weights_wt = [data['weight'] for _, _, data in wt_network.edges(data=True)]

# Plot the distribution of edge weights for both networks
plt.figure(figsize=(10, 6))
plt.hist(edge_weights_five_x_fad, bins=30, alpha=0.5, label='5xFAD')
plt.hist(edge_weights_wt, bins=30, alpha=0.5, label='WT')
plt.legend()
plt.title('Edge Weight Distribution Comparison')
plt.xlabel('Edge Weight (Correlation)')
plt.ylabel('Frequency')
plt.show()
