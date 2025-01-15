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
```
