# Skmeans
The datasets are uploaded in the mydata folder. The last column denotes the ground truth. The R code for the S-k-means algorithm and other peer algorithms are uploaded in the functions.R file. An example run is given in Example.pdf and Example.Rmd.

The details for S-k-means are as follows:

#### Inputs:
X       : An n * p matrix to be clustered, whose rows denote the data points.
M       : A k * p matrix, whose rows denote the inital cluster centroids.
itermax : Maximum number of itearions to run. Default is 30.
#### Outputs:

label     : An n-length vector denoting the class labels. 
centroids : A k * p matrix, whose rows denote the cluster centroids.


If you use the code or the datasets, please acknowledge so by citing the following article.
S. Chakraborty, S. Das, k−Means clustering with a new divergence-based distance metric: Convergence and
performance analysis, Pattern Recognition Letters (2017), https://doi.org/10.1016/j.patrec.2017.09.025
