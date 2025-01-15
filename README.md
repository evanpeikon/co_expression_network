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
# Visualize the network 
plt.figure(figsize=(10, 8))
pos = nx.spring_layout(combined_network, seed=42)
nx.draw(combined_network, pos,with_labels=False, node_size=20, edge_color="lightblue")
plt.title("Gene Co-expression Network")
plt.show()

# identify hub genes in network
degrees = dict(combined_network.degree())
sorted_genes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
print("Top 10 hub genes:")
print(sorted_genes[:10])
```
Which produces the following output:

<img width="600" alt="Screenshot 2025-01-15 at 1 37 42 PM" src="https://github.com/user-attachments/assets/ee66bf64-d852-4b3f-a0de-dc37f43878c8" />

- ```Top 10 hub genes:[('Trappc10', 32), ('Gpr107', 32), ('Mkrn2', 31), ('Tfe3', 31), ('Brat1', 30), ('Dlat', 30), ('Hnrnpd', 30), ('Nalcn', 30), ('Dgke', 30), ('Cdh4', 29)]```
