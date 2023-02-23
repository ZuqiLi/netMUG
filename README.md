# netMUG
A network-guided multi-view clustering framework

netMUG (**net**work-guided **MU**lti-view clusterin**G**) is a novel pipeline that clusters samples on the basis of these sample-specific feature interactions or wirings. In particular, in the presence of 2-view data, multi-view features are jointly selected via SmCCNet [1], based on canonical correlations and additional extraneous variable. ISNs (individual-specific networks) [2] are constructed from the selected features, taking as a measure of edge strength the overall correlation between every pair of features and the extraneous data. The Euclidean distance metric representing dissimilarities between ISNs is fed into Ward’s hierarchical clustering, using the Dynamic Tree Cut R library to automatically derive the number of clusters [3].

### Installation
To install netMUG in your R environment, you can either download the code files or use the following code in R:
```
# Install package
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ZuqiLi/netMUG")

# Load package
library(netMUG)
```
The current netMUG is developed under R version 4.2.1 with the following packages:
- parallel (4.2.1)
- devtools (2.4.3)
- dynamicTreeCut (1.63.1)
- SmCCNet (0.99.0)

### Usage
There're 2 ways to use netMUG: call the all-in-one-go function or break it down to steps.
#### Strategy 1: All-in-one-go function
```
# X is the first data view of shape [n x p1]
# Y is the second data view of shape [n x p2]
# Z is the extraneous variable of shape [n]
# l1, l2 are the sparsity parameters (more explanation can be found in SmCCNet), which can be tuned via cross validation
# s1, s2 are the subsampling parameters (more explanation can be found in SmCCNet), which can be tuned via cross validation
library(netMUG)
res = netMUG(X, Y, Z, l1, l2, s1, s2)
# netMUG returns a list: the selected features from X, the selected features from Y, ISNs, and the final clustering.
```
#### Strategy 2: step-by-step pipeline
```
library(netMUG)
# Step 1: select multi-view features informed by an extraneous variable
smccnet <- selectFeatures(X, Y, Z, l1, l2, s1, s2)
Xsub <- X[, smccnet$featureX]
Ysub <- Y[, smccnet$featureY]

# Step 2: build ISNs from the selected features
V <- cbind(Xsub, Ysub)
ISNs <- buildInfISNs(V, Z, nCores = 1)

# Step 3: compute distances between ISNs
dis <- computeDist(ISNs)

# Step 4: Ward's hierarchical clustering with Dynamic Tree Cut
dendro <- hclust(as.dist(dis_m), method = "ward.D2")
clust <- cutreeDynamic(dendro, minClusterSize = 1, distM = dis, 
                      deepSplit = 0)
clust <- as.factor(clust)
```

### References
> [1] Shi WJ, Zhuang Y, Russell PH et al. Unsupervised discovery of phenotype-specific multi-omics networks. Bioinformatics 2019;35:4336–43.\
> [2] Kuijjer ML, Tung MG, Yuan G et al. Estimating Sample-Specific Regulatory Networks. iScience 2019;14:226–40.\
> [3] Langfelder P, Zhang B, Horvath S. Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 2008;24:719–20.
