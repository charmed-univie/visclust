%% Default parameters
% load iris dataset
[X, y] = iris_dataset;
X=X';
[y,~]=find(y~=0); % convert vector containing correct clustering to suitable format for evaluation
NumClusters=3;
% cluster using default parameters and optional 'showimg' parameter
[prediction, projector]=visclust(X,NumClusters,'showimg','true');

% if the target vector is given, you can compute the accuracy, Rand index or adjusted Rand index as follows
disp("Clustering accuracy: "+evaluation(prediction, y,"ACC"))
disp("Clustering Rand index: "+evaluation(prediction, y,"RI"))
disp("Clustering adjusted Rand index: "+evaluation(prediction, y,"ARI"))

% % Optional parameter configurations
% % Clustering without division or number of clusters
% [prediction, projector]=visclust(X);
%
% % Modifying 'division' parameter
% % uniform division, default
% [prediction, projector]=visclust(X,NumClusters,'division','0');
% % predefined division
% %[prediction, projector]=visclust(X,NumClusters,'division',[0.5,0.5]);
% 
% % Modifying 'thresh' parameter
% [prediction, projector]=visclust(X,NumClusters,'thresh',0.2);
% 
% % Modifying 'scaling' parameter
% [prediction, projector]=visclust(X,NumClusters,'scaling',1.2);
% 
% % Modifying 'subsampling' parameter
% [prediction, projector]=visclust(X,NumClusters,'subsampling',50);
% 
% % Modifying 'projectors' parameter
% % random orthogonal projections, default
% [prediction, projector]=visclust(X,NumClusters,'projections','random');
% % PCA projection
% [prediction, projector]=visclust(X,NumClusters,'projections','pca');
% % t-SNE dimension reduction
% [prediction, projector]=visclust(X,NumClusters,'projections','tsne');
% % custom projectors
% P2=rand(4,2);
% P3=rand(4,3);
% P=cell(1,2);
% P{1,1}=P2;
% P{1,2}=P3;
% [prediction, projector]=visclust(X,NumClusters,'projections',P);
% 
% % Modifying 'method' parameter
% % all methods, default
% [prediction, projector]=visclust(X,NumClusters,'method','all');
% %'all2' - Simultanious clustering in 2D and Partial clustering in 2D
% [prediction, projector]=visclust(X,NumClusters,'method','all2');
% %'all3' - Simultanious clustering in 3D and Partial clustering in 3D
% [prediction, projector]=visclust(X,NumClusters,'method','all3');
% %'vis2' - Simultanious clustering in 2D
% [prediction, projector]=visclust(X,NumClusters,'method','vis2');
% %'vis3' - Simultanious clustering in 3D
% [prediction, projector]=visclust(X,NumClusters,'method','vis3');
% %'pvis2' - Partial clustering in 2D
% [prediction, projector]=visclust(X,NumClusters,'method','pvis2');
% %'pvis3' - Partial clustering in 3D
% [prediction, projector]=visclust(X,NumClusters,'method','pvis3');
% 
% % Modifying 'showimg' parameter
% [prediction, projector]=visclust(X,NumClusters,'showimg','true');
% [prediction, projector]=visclust(X,NumClusters,'showimg','true','method','vis3');
% 
% % Modifying 'imgpath' parameter
% [prediction, projector]=visclust(X,NumClusters,'showimg','./test');
% [prediction, projector]=visclust(X,NumClusters,'showimg','./test','method','vis3');