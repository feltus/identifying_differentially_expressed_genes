# identifying_differentially_expressed_genes

Certainly! I'd be happy to help you create a differential gene expression activity using normalized data from NCBI. This activity will involve several steps, from data acquisition to analysis and visualization. Let's break it down into manageable parts:

1. Data Acquisition
2. Data Preprocessing
3. Differential Expression Analysis
4. Visualization
5. Interpretation

Let's go through each step:

1. Data Acquisition:

First, we'll need to download a suitable dataset from NCBI's Gene Expression Omnibus (GEO). For this example, let's use a dataset comparing gene expression in different conditions.


Based on the search results, we can use the Gene Expression Omnibus (GEO) database from NCBI to find a suitable dataset for our differential gene expression activity. Let's choose a dataset and proceed with the activity.

For this example, we'll use a hypothetical dataset comparing gene expression in normal and diseased tissue samples. We'll assume this dataset is already normalized and available in GEO.

2. Data Preprocessing:

For this activity, we'll use R to process and analyze the data. Here's a script to download and preprocess the data:

```R
# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")

library(GEOquery)
library(limma)

# Download the dataset (replace with actual GEO accession number)
gset <- getGEO("GSE12345", GSEMatrix = TRUE, AnnotGPL = TRUE)
if (length(gset) > 1) idx <- grep("GPL\\d+", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Extract expression data
ex <- exprs(gset)

# Extract phenotype data
pd <- pData(gset)

# Create a design matrix (assuming two groups: normal and diseased)
design <- model.matrix(~ factor(pd$group))
colnames(design) <- c("Intercept", "Diseased")

# Fit the model
fit <- lmFit(ex, design)
fit <- eBayes(fit)

# Get top differentially expressed genes
top_genes <- topTable(fit, coef = "Diseased", number = Inf)
```

3. Differential Expression Analysis:

The `topTable` function in the script above performs the differential expression analysis. It ranks genes based on their differential expression between the two conditions.

4. Visualization:

Let's create a volcano plot to visualize the results:

```R
# Create a volcano plot
library(ggplot2)

volcano_data <- data.frame(
  gene = rownames(top_genes),
  log2FoldChange = top_genes$logFC,
  pvalue = top_genes$P.Value
)

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = abs(log2FoldChange) > 1 & pvalue < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal()
```

5. Interpretation:

Now, let's interpret the results:

a) Examine the top differentially expressed genes:
```R
head(top_genes)
```

b) Count the number of significantly differentially expressed genes:
```R
sum(top_genes$adj.P.Val < 0.05)
```

c) Identify genes with large fold changes:
```R
large_fc_genes <- top_genes[abs(top_genes$logFC) > 2 & top_genes$adj.P.Val < 0.05, ]
head(large_fc_genes)
```

Activity Questions:

1. How many genes are significantly differentially expressed (adjusted p-value < 0.05)?
2. What are the top 5 most upregulated and downregulated genes?
3. Are there any genes with very large fold changes (|log2FC| > 2)? What might be their biological significance?
4. Looking at the volcano plot, how would you describe the overall pattern of gene expression changes?
5. Choose one of the top differentially expressed genes and research its function. How might its altered expression contribute to the disease state?

This activity provides a hands-on experience with differential gene expression analysis using real data from NCBI's GEO database. It covers the entire process from data acquisition to interpretation, giving students a practical understanding of bioinformatics techniques used in genomics research.
