# Skmeans
The datasets are uploaded in the mydata folder. The last column denotes the ground truth. The R code for the S-k-means algorithm and other peer algorithms are uploaded in the functions.R file. An example run is given in Example.pdf and Example.Rmd.

## Details
The details for S-k-means are as follows:

#### Inputs:
1. X       : An n * p matrix to be clustered, whose rows denote the data points.
2. M       : A k * p matrix, whose rows denote the inital cluster centroids.
3. itermax : Maximum number of itearions to run. Default is 30.
#### Outputs:

1. label     : An n-length vector denoting the class labels. 
2. centroids : A k * p matrix, whose rows denote the cluster centroids.


## Paper
S. Chakraborty, S. Das, k−Means clustering with a new divergence-based distance metric: Convergence and
performance analysis, Pattern Recognition Letters (2017), https://doi.org/10.1016/j.patrec.2017.09.025
