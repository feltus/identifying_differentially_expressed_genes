# lab-identifying_differentially_expressed_genes

# Lab Overview
In this lab, we will obatin RNAseq data from NCBI and perform a differential gene expression analysis.

# Getting help

Review how to use Palmetto here: https://docs.rcd.clemson.edu/palmetto/

Useful Prometheus prompts:

```
* What are the molecular biology steps to construct and sequence an RNA sample for RNAseq analysis?
* What are the differences between count, FPKM, and TPM units for RNAseq?
* Please explain how DESeq2 works including and explanation of the negative binomial distribution.
```

# Lab Objectives
* Open R-studio on the Palmetto Clsuter
* Perform Differentially Expressed Gene (DEG) Analysis on Palmetto.
* Visualization and interpret the results.

# Task A. Experimental setup.

***Step A.  Launch R studio, create a working directory, and install R packages***
First, we'll need to download a suitable dataset from NCBI's Gene Expression Omnibus (GEO). For this example, let's use a dataset comparing gene expression in different conditions. 

* Log into Palmetto and create a working directory.
Access Palmetto using ondemand: https://ondemand.rcd.clemson.edu/. 

* Start an R-studio server and open a seperate terminal on Palmetto.
You will find this in the interactive sessions link.  You can try a default server or modify to 16GB of RAM, two CPU cores, and 12 hours of walltime.

* Open a terminal and create a working directory in */scratch*. Dont forget that the directory will eed to be nested inside a directory with your user name.

* Go to the R-studio console and clear it (CTRL-L).

* Set your working directory to the one you created using the *session* drop-down menu. 

**Step B. Install R packages***

* Execute these R commands in the R console.

```R
# Install and load required packages
BiocManager::install(version = '3.19')
BiocManager::install("GEOquery")
BiocManager::install("limma")

# Load libraries
library(GEOquery)
library(limma)
```

**Step C. Perform the DEG analysis***

```
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
**Step C. Visualize the results***

#Let's create a volcano plot to visualize the results. Execute these commands in R.

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

**Step C.Interpret the results***
Examine the top differentially expressed genes:
```R
head(top_genes)
```

Count the number of significantly differentially expressed genes:
```R
sum(top_genes$adj.P.Val < 0.05)
```

Identify genes with large fold changes:
```R
large_fc_genes <- top_genes[abs(top_genes$logFC) > 2 & top_genes$adj.P.Val < 0.05, ]
head(large_fc_genes)
```

# Task A. Experimental setup.
Upload your volcao plot and a text file with the top DEGs into the Praxis LXP VM.
