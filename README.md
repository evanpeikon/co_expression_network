# 🧬 Unraveling Gene Co-Expression Networks: Insights into Transcriptional Regulation
## 🧬 Introduction 
Co-expression network analysis is a powerful approach for identifying groups of genes, referred to as "modules," that exhibit coordinated expression patterns across a dataset. These modules often represent functionally related genes, or sets of genes, that work together in shared biological pathways or cellular processes. By capturing patterns of gene activity, co-expression networks enable researchers to infer relationships between gene expression and biological traits or conditions. This makes them invaluable tools for uncovering the mechanisms driving disease progression, developmental processes, or responses to environmental changes.

One of the most widely adopted tools for constructing co-expression networks is Weighted Gene Co-expression Network Analysis (WGCNA). This method leverages pairwise correlations between gene expression levels to group genes into clusters, or modules, of tightly interconnected genes. Additionally, WGCNA highlights key genes within each module—known as hub genes—that are highly connected and often play pivotal regulatory roles.

This guide will introduce you to the principles of co-expression network analysis, walk you through constructing a co-expression network using WGCNA, and demonstrate how to identify key modules and hub genes associated with specific biological traits. Let’s get started!

## 🧬 How Does Co-Expression Network Analysis Differ from PPI Network Analysis?
As we’ve established, co-expression networks are a valuable tool for exploring transcriptional regulation and identifying coordinated patterns of gene expression. By grouping co-expressed genes into modules, researchers can uncover functionally related gene sets, correlate these modules with external traits, and pinpoint key regulatory genes. However, you might wonder how this approach differs from protein-protein interaction (PPI) network analysis—a topic I previously discussed in a tutorial titled [From Genes to Networks: Uncovering Biological Insights through Protein-Protein Interactions](https://github.com/evanpeikon/PPI_Network_Analysis)/

While both approaches reveal relationships between biological entities, they do so at different levels of the biological hierarchy. Co-expression networks focus on transcriptional regulation, capturing the coordination of gene expression across samples. In contrast, PPI networks highlight the physical and functional interactions between proteins. Essentially, co-expression networks uncover the upstream regulatory relationships that drive gene expression, while PPI networks represent the downstream interactions between the proteins encoded by those genes.

Additionally, co-expression analysis provides a broader perspective, capturing genes that are functionally linked through shared expression patterns, even if the proteins they code for do not physically interact. For example, genes A, B, and C may form a co-expression module due to being regulated by a common transcription factor. While their proteins might eventually interact to form a functional complex, co-expression analysis focuses on their shared transcriptional regulation. In contrast, PPI networks explore the mechanistic relationships between the proteins themselve, offering a detailed view of how they collaborate physically to perform biological functions. Together, these two approaches provide a comprehensive understanding of biological systems, linking upstream regulation to downstream molecular mechanisms.

## 🧬 Constructing a Gene Co-Expression Network: Background and Approach
Now, in this section, I’m going to walk you through how to create a gene co-expression network. Before we dive into the practical steps, it’s important to first understand a key characteristic of this type of analysis. Unlike a protein-protein interaction (PPI) network, which requires a pre-generated list of differentially expressed genes (DEGs) as input, co-expression network analysis can be performed even if you only have sample data from a single condition. This makes co-expression analysis particularly versatile because it does not rely on prior differential expression analysis. Instead, it identifies patterns of shared expression among genes across your dataset, revealing modules of co-regulated genes that can offer unique insights into the regulatory framework of the biological condition under study.

For this tutorial, we’ll be using a file named ```countlist_5xFAD```. This dataset includes gene expression data from 20 samples: 10 from 5xFAD mutant mice and 10 from wild-type (WT) control mice. This file is the output of a [prior analysis I performed](https://github.com/evanpeikon/mouse_AD_models), exploring the similarities, differences, and overlapping features between three mutant mouse models of Alzheimer’s disease. While we are working with this specific dataset, the methods I’ll demonstrate can be applied to any dataset formatted as a CSV file where rows represent genes and columns represent samples, with only minor adjustments to the code.

> Note: You can use [this code](https://github.com/evanpeikon/co_expression_network/blob/main/code/get_5xFAD_data.py) to retrieve the 5xFAD mouse data (+ it's corresponding WT control data) from GEO, then clean and reformat it, so the resultant output file ```countlist_5xFAD``` is the same dataset I'm working with in this tutorial. 

Notably, when you have a dataset with two experimental conditions, such as treatment and control, you have a choice in how you approach the co-expression analysis. You can either (1) combine all the data into a single analysis or (2) split the dataset by condition. Each approach has its merits depending on the research question you are addressing. A combined co-expression network analysis, which pools the data from both conditions, is typically used when the goal is to identify gene modules that are conserved across the conditions. These conserved modules often represent core biological processes or pathways that are not significantly altered by the experimental treatment (or in this case, disease progression).

On the other hand, separate co-expression analyses for the treatment and control groups are particularly useful if you are interested in comparing network structures and identifying condition-specific modules or hub genes. By analyzing each group independently, you can uncover gene modules that are unique to a specific condition. For instance, modules found only in the 5xFAD group may represent biological pathways disrupted in Alzheimer’s disease, while those unique to the WT group may reflect normal regulatory networks. This approach also allows you to compare network properties such as connectivity, hub gene identification, and module composition between the two groups, which can reveal critical differences in gene regulation under disease versus normal conditions.

In practice, researchers often start with a combined analysis to identify shared and genotype-associated modules. Afterward, they may refine their findings by performing group-specific analyses to further explore differences in network topology or module composition, which is the approach I’ll take in this tutorial. This stepwise approach—starting with a combined network and following up with group-specific networks—not only provides a broader perspective on the biological system but also ensures that subtle, condition-specific effects are not overlooked.

Now, in the code block below, we'll start by loading the libraries used in this tutorial, importing our data, and then calculating pairwise correlations between gene expression:

> Note: I'm going to limit this analysis to the first 100 genes in my countlist_5xFAD file. This is to reduce computational load, as co-expression analysis can be very computationally expensive when using large datasets (and i'm using a personal machine to run the code in this tutorial).

```python
# Import all libraries for this tutorial
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# Load the count list CSV
file_path = 'countlist_5xFAD.csv'  # Path to your CSV file
data = pd.read_csv(file_path, index_col=0)

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
```
<img width="1198" alt="Screenshot 2025-01-15 at 10 53 54 AM" src="https://github.com/user-attachments/assets/bf6195bd-6d47-44f5-b7d2-3b699eb844e9" />

As you can see, the output above is a correlation matrix generated from our gene expression dataset using Spearman correlation. Note, that the matrix has the same set of genes in both the rows and columns (it's an 100x100 matrix). Each cell in the matrix represents the pairwise Spearman correlation between two genes based on their expression profiles across all samples. The diagonal values are all 1.0 because each gene is perfectly correlated with itself, and the off diagonal values represent the Spearman correlation coefficient between pairs of genes (the coefficient ranges from -1, which is a perfect negative correlation, to 1, which is a perfect positive correlation). Additionally, you'll note that some cells in the matrix show an NaN value. This occurs when at least one of the genes in a pair has constant expression across all samples (i.e., no variation) or when the gene has a very low (near zero) expression level across samples. 

Now that we've created our correlation matrix, we'll use it as input to create a co-expression network, as demonstrated below:

```python
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
```
Which produces the following output:
- ```Combined Network: 69 nodes and 583 edges```

The output above tells us that our of the 100 gene in our dataset, 69 genes are included in our co-expression network and there are 586 edges (connections) between these genes. The reason we have 69 nodes in the network, and not 100, is that some genes may not have have a significant correlation (above ```threshold = 0.7```) with any other genes, and as a result they are excluded from the network. The edges in our network represent pairwise relationships between genes based on their correlation. In other words,  an edge exists between two genes if the absolute value of their correlation is greater than or equal to the threshold of 0.7. 

Now, let's visualize the graph of our network along with the top 10 hub genes based on degree:

```python
# Visualize the network with all gene names as labels
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(combined_network, seed=42)
nx.draw(combined_network, pos, with_labels=False, node_size=20, edge_color="lightblue")
nx.draw_networkx_labels(combined_network, pos, labels={node: node for node in combined_network.nodes}, font_size=7)
plt.show()

# identify hub genes in network
degrees = dict(combined_network.degree())
sorted_genes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
print("Top 10 hub genes:")
print(sorted_genes[:10])
```
Which produces the following output:

<img width="600" alt="Screenshot 2025-01-15 at 3 05 22 PM" src="https://github.com/user-attachments/assets/0320b6e7-475b-4fc1-bc01-5c1f76425ca6" />

- ```Top 10 hub genes:[('Trappc10', 32), ('Gpr107', 32), ('Mkrn2', 31), ('Tfe3', 31), ('Brat1', 30), ('Dlat', 30), ('Hnrnpd', 30), ('Nalcn', 30), ('Dgke', 30), ('Cdh4', 29)]```

In the figure above you can see the gene co-expression network with nodes labeled by gene name. Additionally, the code outputs the top 10 hub genes in the co-expression network, ranked by the degree (i.e., the number of connections) each gene has. For example, Trappc10 has 32 connections (edges), meaning it is co-expressed with 33 other genes are or above the threshold correlation of 0.7. 

Using the degree to rank nodes is a simple and intuituve way to identify hub genes, as it reflects how many other genes a particular gene is co-expressed with above the threshold. However, we can also use more sophisticated metrics like degree centrality, clustering coefficient, betweenness centrality, and eigenvector centrality, which provide additional insights into the network structure and relationships between nodes, which we'll explore in the next section below. 

## 🧬 Exploring Gene Co-Expression Network Connectivity
In the last sub-section, I introduced the metrics degree, degree centrality, clustering coefficient, betweenness centrality, and eigenvector centrality. Each of these metrics has it's own strengths and advantages, which I'll briefly outline below:

- ```Degree```: The degree is the number of edges (connections) that a given node has, which reflects how many other genes a given gene is co-expressed with. This measurement provides a simple and fast measurement of network connectivity, despite it's simplicity. 
- ```Degree Centrality:``` Degree centrality normalizes the degree of a node by dividing it by the maximum possible degree in the network. As a result, degree centrality gives a relative measure of how connected a node is compared to others in the network, making it a useful metric when comparing nodes in networks of different sizes.
- ```Clustering Coefficient:``` The clusterng coefficient of a node measures the connectivity of a node's neighbors, indicating how close they are to forming a complete (i.e. fully connected) graph. This measure is useful for understanding localized network structure and helps identify genes that are part of tightly connected modules.
- ```Betweeness Centrality:``` Betweeness centrality measures the frequency at which a node appears on the shortest paths between other nodes in the network, and as a result it identifies nodes that act as bridges or bottlenecks in the network. Nodes that rank highly on betweenness centrality may serve as critical regulators or mediators in the network. 
- ```Eigenvector Centrality:``` This metric measures a node's influence in the network based on the influence of its neighbors. High eigenvector centrality indicates a gene is connected to many well-connected nodes, and as a result it identifies nodes with influence beyond their immediate connections (i.e. nodes with "global" influence). 

For co-expression networks, degree, degree centrality, and clustering coefficient are most commonly used. However, if you're analyzing regulatory hubs or genes involved in key pathways, eigenvector or betweenness centrality could provide deeper insights. In the code below, we'll identify the top 10 hub genes baseds on degree centrality, clustering coefficient, betweeness centrality, and eigenvector centrality, and then we'll visualize sub-plots for each of these sets of hub genes:

```python
# Identify hub genes for each metric
degree_centrality = nx.degree_centrality(combined_network)
sorted_centrality = sorted(degree_centrality.items(), key=lambda x: x[1], reverse=True)
clustering = nx.clustering(combined_network)
sorted_clustering = sorted(clustering.items(), key=lambda x: x[1], reverse=True)
betweenness = nx.betweenness_centrality(combined_network)
sorted_betweenness = sorted(betweenness.items(), key=lambda x: x[1], reverse=True)
eigenvector = nx.eigenvector_centrality(combined_network)
sorted_eigenvector = sorted(eigenvector.items(), key=lambda x: x[1], reverse=True)

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
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

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
```
<img width="600" alt="Screenshot 2025-01-15 at 3 24 27 PM" src="https://github.com/user-attachments/assets/098b2ba5-8aac-4605-b1c5-0cbeb603a551" />

As you can see in the figure above, analysis identified key genes in the co-expression network based on four centrality measures, each of which highlights different aspects of network topology, revealing the roles of these genes within the broader network structure.
- Top Genes By Degree Centrality: The genes with the highest degree centrality, such as Trappc10 and Gpr107, are highly connected within the network, indicating they may act as hubs coordinating interactions among other genes. These hubs are critical for the overall cohesion and robustness of the network.
- Top Genes by Clustering Coefficient: Several genes, including Cdc45, Fap, and Btbd17, had clustering coefficients of 1.0, suggesting they form tightly interconnected clusters. This implies that these genes might belong to highly specialized or functionally coherent modules within the network.
- Top Genes by Betweenness Centrality: Genes like Mkrn2, Comt, and Rtca exhibited high betweenness centrality, indicating their importance as bridges or intermediaries connecting different parts of the network. These genes could play a pivotal role in the flow of information or regulation between gene modules.
- Top Genes by Eigenvector Centrality: Genes such as Gpr107 and Trappc10 were ranked highest by eigenvector centrality, reflecting their influence in the network through their connections to other well-connected genes. These genes may be critical for maintaining the overall regulatory structure of the network.

Now that we've explored the combined network for our treatment (5xFAD mice) and control (corrsesponding WT mice) samples, we'll analyse the co-expression networks for each of these groups seperatley.

## 🧬 Exploring Gene Co-Expression Networks For Treatment and Control Samples (Seperatley)
As previously mentioned separate co-expression analyses for the treatment and control groups are particularly useful if you are interested in comparing network structures and identifying condition-specific modules or hub genes. By analyzing each group independently, you can uncover gene modules that are unique to a specific condition. For instance, modules found only in the 5xFAD group may represent biological pathways disrupted in Alzheimer’s disease, while those unique to the WT group may reflect normal regulatory networks. This approach also allows you to compare network properties such as connectivity, hub gene identification, and module composition between the two groups, which can reveal critical differences in gene regulation under disease versus normal conditions.

In the code block below, we'll create gene co-expression networks for our 5xFAD mouse samples and WT control samples, then identify the hub genes in each network:

```python
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
```
<img width="900" alt="Screenshot 2025-01-15 at 3 37 49 PM" src="https://github.com/user-attachments/assets/c62176ba-fca4-49ed-ab6b-6d5c0e48f9bf" />

Now that we've created two seperate gene co-expression networks, there are several interesting analyses you can perform. To start, let's compare the density of both networks, which will tell us how well connected the genes are in each group:

```python
density_five_x_fad = nx.density(five_x_fad_network)
density_wt = nx.density(wt_network)
print(f"Density of 5xFAD network: {density_five_x_fad}")
print(f"Density of WT network: {density_wt}")
```
Which produces the following output:
- ```Density of 5xFAD network: 0.18991640543364682```
- ```Density of WT network: 0.253140871424752```

Based on the results above, we can see that the density of the 5xFAD network is relatively low, indicating that the gene co-expression relationships in this network are sparser. This may suggest that the regulatory relationships between genes in the 5xFAD condition are more selective or disrupted compared to the WT condition, which may reflect altered or less cohesive interactions in the disease state (we should be careful to overinterpret these findings as we are not using the entire dataset for these analysis). Additionally, we can see that the WT network has a higher density, implying that the gene co-expression relationships in the control (wild-type) samples are more interconnected, which could indicate a more robust and tightly regulated gene network under normal physiological conditions. 

Next, we'll look at the average degree of nodes in each network, which can indicate the overall connectivity in each condition:

```python
avg_degree_five_x_fad = sum(dict(five_x_fad_network.degree()).values()) / len(five_x_fad_network.nodes)
avg_degree_wt = sum(dict(wt_network.degree()).values()) / len(wt_network.nodes)
print(f"Average degree in 5xFAD network: {avg_degree_five_x_fad}")
print(f"Average degree in WT network: {avg_degree_wt}")
```
- ```Average degree in 5xFAD network: 16.522727272727273```
- ```Average degree in WT network: 21.770114942528735```

As expected bases on the previous analysis, the 5xFAD network has a lower average degree than the WT network. Now that we've highlighted some differences between the two networks, we'll perform a cross-network hub gene analysis to see if there are any overlapping hub genes between the two networks:

```python
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
```
- ```Common hub genes: {'Trappc10', 'Hnrnpd', 'Gpr107'}```
- ```Unique 5xFAD hub genes: {'Tfe3', 'Itga5', 'Nalcn', 'Acvr1b', 'Xpo6', 'Dgke', 'Scmh1'}```
- ```Unique WT hub genes: {'Mkrn2', 'Raf1', 'Sox9', 'Brat1', 'Fer', 'Gna12', 'Gnai3'}```

The analysis above reveals three common hub genes between the 5xFAD and WT networks including Hnrnpd, Trappc10, and Gpr107, which suggests these genes may play a central role in maintaining core co-expression relationships in both conditions. However, the networks also exhibit distinct unique hub genes, with seven genes exclusive to the 5xFAD network  and seven genes unique to the WT network. These unique hub genes likely reflect condition-specific regulatory roles, with 5xFAD-specific hubs potentially driving processes related to Alzheimer’s pathology, while WT-specific hubs maintain normal regulatory functions. 

Finally, we'll look at edge weight distributions, which allows us to compare the distribution of correlation strengths between the 5xFAD and WT networks, telling us if certain gene pairs have stronger or weaker co-expression in one group compared to another:

```python
# Extract edge weights from 5xFAD and WT networks
edge_weights_five_x_fad = [data['weight'] for _, _, data in five_x_fad_network.edges(data=True)]
edge_weights_wt = [data['weight'] for _, _, data in wt_network.edges(data=True)]
plt.figure(figsize=(10, 6))
plt.hist(edge_weights_five_x_fad, bins=30, alpha=0.5, label='5xFAD')
plt.hist(edge_weights_wt, bins=30, alpha=0.5, label='WT')
plt.legend()
plt.title('Edge Weight Distribution Comparison')
plt.xlabel('Edge Weight (Correlation)')
plt.ylabel('Frequency')
plt.show()
```
<img width="644" alt="Screenshot 2025-01-15 at 3 50 35 PM" src="https://github.com/user-attachments/assets/8672fdda-3b1e-4da7-9bf7-5a8cc20a22f6" />

The edge weight distributions reveal notable differences in co-expression patterns between the 5xFAD and WT networks. The WT network exhibits a higher frequency of strong positive correlations (edge weights in the 0.75 to 1.0 range), suggesting more robust and coordinated gene expression relationships under normal conditions. In contrast, the 5xFAD network shows a stronger presence of negative correlations (edge weights around -0.75. These findings suggest that Alzheimer's pathology in the 5xFAD model may disrupt the normal regulatory architecture, reducing positive co-expression while introducing dysregulated, oppositional interactions between genes, which may contribute to disease-specific processes.

# 🧬 If You Enjoy This Guide...
## You May Also Enjoy the Following DIY Guides...
- [Python Fundamentals For Biologists](https://github.com/evanpeikon/Python_Fundamentals_Biology)
- [Bash Fundamentals For Bioinformatics](https://github.com/evanpeikon/bash_fundamentals)
- [Assorted Functions In Numerical Methods & Systems Biology](https://github.com/evanpeikon/systems_biology)
- [Introduction to Functional Enrichment Analysis](https://github.com/evanpeikon/functional_enrichment_analysis)
- [From Genes to Networks: Uncovering Biological Insights through Protein-Protein Interactions](https://github.com/evanpeikon/PPI_Network_Analysis)





