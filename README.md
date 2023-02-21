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

### References:
> [1] Shi WJ, Zhuang Y, Russell PH et al. Unsupervised discovery of phenotype-specific multi-omics networks. Bioinformatics 2019;35:4336–43.
> [2] Kuijjer ML, Tung MG, Yuan G et al. Estimating Sample-Specific Regulatory Networks. iScience 2019;14:226–40.
> [3] Langfelder P, Zhang B, Horvath S. Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. Bioinformatics 2008;24:719–20.
