# VisClust

This is visClust, a novel clustering algorithm based on theory from mathematical imaging and linear projections.

This repository contains the implementation of [visClust: A visual clustering algorithm based on orthogonal projections](https://arxiv.org/abs/2211.03894), Anna Breger, Clemens Karner, Martin Ehler, 2023.

**Written by Anna Breger and Clemens Karner**\
University of Vienna, Faculty of Mathematics\
Vienna, Austria\
https://homepage.univie.ac.at/anna.breger/
https://homepage.univie.ac.at/clemens.karner/

This software is licensed under the MIT license, see the file 'LICENSE'.

For all questions, bugs and suggestions please email
clemens.karner@univie.ac.at or anna.breger@univie.ac.at

# Usage instructions

The file 'example.m' contains all possible parameter configurations applied to the Iris dataset.

Example:
```
X = iris_dataset'
[prediction, projector] = visclust(X,3)
```

**VisClust** Unsupervised classification via visual clustering
   - IDX = VISCLUST(X, OPTIONAL) returns a suggested number of clusters for N data points of dimension D, stored in an N-by-D matrix X.

   - IDX = **visclust**(X, NUMCLUSTERS, OPTIONAL) partitions N data points of dimension D, stored in an N-by-D matrix X, into K clusters. VISCLUST returns a vector IDX of length N containing the cluster indices of each data point.

   - [IDX, P] = **visclust**(X, NUMCLUSTERS, OPTIONAL) returns a cell P containing the optimal projection(s) identified by visClust.

   - [ ... ] = **visclust**(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies optional parameter name/value pairs to control visClust. Parameters are:

      - 'division'... Wanted division of clusters. Choices are:
         - '0'... Uniform division (default).
         - K-by-1 matrix... Wanted relative division, e.g. division = [0.4,0.4,0.2].

      - 'thresh... Allowed deviation of wanted division, ideally between 0.01 and 0.3. Default: 0.1 (corresponds to 10%)

      - 'scaling'... Parameter to control the standard deviation in the Gaussian filters, ideally between 0.5 and 2. Default: 1.25

      - 'subsampling'... Number of data points clustered without subsampling. Default: 10000
    
      - 'projections'... Projections to be used for clustering. Choices are:
         - 'random'... Random orthogonal projections (default).
         - 'pca'... Projector obtained by Principal Component Analysis.
         - 'tsne'... Dimension reduction using the t-Distributed Stochastic Neighbor Embedding method.
         - M-by-2 cell... M projectors stored in a M-by-2 cell. The first cell dimension stores the matrices of projections to R^2 (shape Dx2), if available the projections to R^3 are stored in the second dimension (shape Dx3).

      - 'method'... Clustering method. Choices are:
         - 'all'... Runs through all methods (default).
         - 'all2'... 2D clustering.
         - 'all3'... 3D clustering.
         - 'vis2'... Simultaneous clustering in 2D.
         - 'vis3'... Simultaneous clustering in 3D.
         - 'pvis2'... Partial clustering in 2D.
         - 'pvis3'... Partial clustering in 3D.

      - 'showimg'... Generates a plot of projected clustering. Choices are:
         - 'false'... No plot (default).
         - 'true'... Generates a plot.
         - 'PATH'... Path to store the generated plot.
