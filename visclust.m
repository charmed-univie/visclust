function [idx,P]= visclust(X,numClusters,varargin)
% VISCLUST Unsupervised classification via visual clustering
%
%   IDX = VISCLUST(X, NUMCLUSTERS, OPTIONAL) partitions N data points of
%   dimension D, stored in an N-by-D matrix X, into K clusters. VISCLUST
%   returns a vector IDX of length N containing the cluster indices of each
%   data point.
%
%   [IDX, P] = VISCLUST(X, NUMCLUSTERS, OPTIONAL) returns a cell P containing
%   the optimal projection(s) identified by visClust.
%
%   [ ... ] = VISCLUST(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%      optional parameter name/value pairs to control visClust.
%      Parameters are:
%
%      'division' - Wanted division of clusters. Choices are:
%         '0' - Uniform division (default).
%          K-by-1 matrix - Wanted relative division, e.g. division = [0.4,0.4,0.2].
%
%      'thresh' - Allowed deviation of wanted division, ideally between 0.01 and 0.3.
%          Default: 0.1 (corresponds to 10%)
%
%      'scaling' - Parameter to control the standard deviation in the Gaussian
%          filters, ideally between 0.5 and 2.
%          Default: 1.25
%
%      'subsampling' - Number of data points clustered without subsampling.
%          Default: 10000
%
%      'projectors' - Projectors to be used for clustering. Choices are:
%         'random' - Random orthogonal projections (default).
%         'pca' - Projector obtained by Principal Component Analysis.
%          M-by-2 cell - M projectors stored in a M-by-2 cell. The first cell
%          dimension stores the matrices of projections to R^2 (shape Dx2), if available
%          the projections to R^3 are stored in the second dimension (shape Dx3).
%
%      'method' - Clustering method. Choices are:
%         'all' - Runs through all methods (default).
%         'all2' - 2D clustering.
%         'all3' - 3D clustering.
%         'vis2' - Simultaneous clustering in 2D.
%         'vis3' - Simultaneous clustering in 3D.
%         'pvis2' - Partial clustering in 2D.
%         'pvis3' - Partial clustering in 3D.
%
%      'showimg' - Generates a plot of projected clustering. Choices are:
%         'false' - No plot (default).
%         'true' - Generates a plot.
%         'PATH' - Path to store the generated plot.
%
%   Example:
%       X = iris_dataset'
%       [prediction, projector] = visClust(X,3)
%
%
% This is VISCLUST, written by Anna Breger and Clemens Karner
% University of Vienna, Faculty of Mathematics
% Vienna, Austria
% Copyright (c) 2022
% https://homepage.univie.ac.at/anna.breger/
% https://homepage.univie.ac.at/clemens.karner/
%
% For all questions, bugs and suggestions please email
% clemens.karner@univie.ac.at or anna.breger@univie.ac.at
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add visClust files to path
addpath("./core_functions/");

%% Initialize variables
defaultScaling = 1.25;
defaultThresh=0.1;
defaultSubsampling=10000;
iternum_k2=5000; % number of iterations (=projections) for 2-dim clustering
iternum_k3=2000; % number of iterations (=projections) for 3-dim clustering

%% Process input arguments
iP = inputParser;
% default input values
defaultMethod = 'all';
validMethods = {'all','all2','all3','vis2','vis3','pvis2','pvis3'};
checkMethod = @(x) any(validatestring(x,validMethods));
defaultProjector= 'random';
validProjectors = {'random','pca'};
checkProjectors = @(x) iscell(x)||any(validatestring(x,validProjectors));
defaultDivision= '0';
defaultImage='false';
checkSubsampling = @(x) assert(isnumeric(x) && (x >= 0)&& (x <= size(iP.Results.X,1)));
% check input values and set default
addRequired(iP,'X',@ismatrix);
addRequired(iP,'numclusters',@isnumeric);
addParameter(iP,'division',defaultDivision,@ismatrix)
addParameter(iP,'projectors',defaultProjector,checkProjectors)
addParameter(iP,'scaling',defaultScaling,@isnumeric)
addParameter(iP,'subsampling',defaultSubsampling,checkSubsampling)
addParameter(iP,'thresh',defaultThresh,@isnumeric)
addParameter(iP,'method',defaultMethod,checkMethod)
addParameter(iP,'showimg',defaultImage,@(x) (ischar(x)||isstring(x)))
parse(iP,X,numClusters,varargin{:})
% create uniform division if requested
div=iP.Results.division;
divset=true;
if (ischar(div)||isstring(div)) && div=="0"
    divset=false;
    div=(ones(iP.Results.numclusters,1)/iP.Results.numclusters)';
end
cl=length(div);

%% Initialize projectors
if iscell(iP.Results.projectors)
    try
        P2=iP.Results.projectors(:,1); % cell array with projections to R^2
        if size(iP.Results.projectors,2)==2
            P3=iP.Results.projectors(:,2); % cell array with projections to R^3
        end
    catch
        error("invalid projector format")
    end
elseif iP.Results.projectors=="random"
    P2 = getProj(iP.Results.X,2,iternum_k2); % cell array with random projections to R^2
    P3 = getProj(iP.Results.X,3,iternum_k3); % cell array with random projections to R^3
elseif iP.Results.projectors=="pca"
    P2 = getProj(iP.Results.X,2,"PCA"); % cell array with PCA to R^2
    P3 = getProj(iP.Results.X,3,"PCA"); % cell array with PCA to R^3
end

%% Split dataset (in case of more than 10 000 data points)
l=size(iP.Results.X,1); % number of data points
if iP.Results.subsampling < l
    idx=ismember(1:l,randsample(1:l,iP.Results.subsampling)); % randomly pick activeData indices and convert to logical array
    DATAA = iP.Results.X(idx,:) ;
    DATAI= iP.Results.X(~idx,:);
    dicta=1:length(iP.Results.X);
    dicta(~idx)=[];
    dicti=1:l;
    dicti(idx)=[];
else
    DATAA=iP.Results.X;
end

%% Start clustering
if iP.Results.method=="all"
    [~,~,~,itertime_k2]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting simultaneous 2D Clustering, worst case runtime: "+num2str(itertime_k2*iternum_k2)+"sec")
    [idx,~,Pi,~,~,tsd,success]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"vis");
    if Pi>0
        P={P2{Pi}};
    else
        P=cell(0,1);
    end
    if  ((~divset&&~success)||divset)&&(cl>2)&&(iscell(iP.Results.projectors)||iP.Results.projectors~="pca")
        disp("Starting partial 2D Clustering, worst case runtime: "+num2str(itertime_k2*iternum_k2*(cl-1))+"sec")
        [idxp,Pp]=visClustPartial_(DATAA,div,iP.Results.scaling,iP.Results.thresh,'pvisp0',iP.Results.projectors,tsd);
        outputVals=unique(idx);
        distr = zeros(cl,1);
        for iout=1:length(outputVals)
            distr(iout)=sum(idx==outputVals(iout))/l;
        end
        distr=sort(distr);
        distrp = zeros(cl,1);
        outputVals=unique(idxp);
        for iout=1:length(outputVals)
            distrp(iout)=sum(idxp==outputVals(iout))/l;
        end
        distrp=sort(distrp);
        div=sort(div);
        
        if (~divset&& max(distrp)<1)||(sum(abs(distrp(end-cl+1:end)-div'))<sum(abs(distr(end-cl+1:end)-div')))
            idx=idxp;
            P=Pp;
        end
    end
elseif iP.Results.method=="all2"
    [~,~,~,itertime_k2]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting simultaneous 2D Clustering, worst case runtime: "+num2str(itertime_k2*iternum_k2)+"sec")
    [idx,~,Pi,~,~,tsd,success]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"vis");
    if Pi>0
        P={P2{Pi}};
    else
        P=cell(0,1);
    end
    if  ((~divset&&~success)||divset)&&(cl>2)&&(iscell(iP.Results.projectors)||iP.Results.projectors~="pca")
        disp("Starting partial 2D Clustering, worst case runtime: "+num2str(itertime_k2*iternum_k2*(cl-1))+"sec")
        [idxp,Pp]=visClustPartial_(DATAA,div,iP.Results.scaling,iP.Results.thresh,'pvisp2',iP.Results.projectors,tsd);
        outputVals=unique(idx);
        distr = zeros(cl,1);
        for iout=1:length(outputVals)
            distr(iout)=sum(idx==outputVals(iout))/l;
        end
        distr=sort(distr);
        distrp = zeros(cl,1);
        outputVals=unique(idxp);
        for iout=1:length(outputVals)
            distrp(iout)=sum(idxp==outputVals(iout))/l;
        end
        distrp=sort(distrp);
        div=sort(div);
        if (~divset&& max(distrp)<1)||(sum(abs(distrp(end-cl+1:end)-div'))<sum(abs(distr(end-cl+1:end)-div')))
            idx=idxp;
            P=Pp;
        end
    end
elseif iP.Results.method=="all3"
    [~,~,~,itertime_k3]=visClust_(DATAA,div,P3,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting simultaneous 3D Clustering, worst case runtime: "+num2str(itertime_k3*iternum_k3)+"sec")
    [idx,~,Pi,~,~,tsd,success]=visClust_(DATAA,div,P3,iP.Results.scaling,iP.Results.thresh,"vis");
    if Pi>0
        P={P3{Pi}};
    else
        P=cell(0,1);
    end
    if  ((~divset&&~success)||divset)&&(cl>2)&&(iscell(iP.Results.projectors)||iP.Results.projectors~="pca")
        disp("Starting partial 3D Clustering, worst case runtime: "+num2str(itertime_k3*iternum_k3*(cl-1))+"sec")
        [idxp,Pp]=visClustPartial_(DATAA,div,iP.Results.scaling,iP.Results.thresh,'pvis3',iP.Results.projectors,tsd);
        outputVals=unique(idx);
        distr = zeros(cl,1);
        for iout=1:length(outputVals)
            distr(iout)=sum(idx==outputVals(iout))/l;
        end
        distr=sort(distr);
        distrp = zeros(cl,1);
        outputVals=unique(idxp);
        for iout=1:length(outputVals)
            distrp(iout)=sum(idxp==outputVals(iout))/l;
        end
        distrp=sort(distrp);
        div=sort(div);
        if (~divset&& max(distrp)<1)||(sum(abs(distrp(end-cl+1:end)-div'))<sum(abs(distr(end-cl+1:end)-div')))
            idx=idxp;
            P=Pp;
        end
    end
elseif iP.Results.method=="vis2"
    [~,~,~,itertime_k2]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting simultaneous 2D Clustering, worst case runtime: "+num2str(itertime_k2*iternum_k2)+"sec")
    [idx,~,Pi]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"visonly");
    if Pi>0
        P={P2{Pi}};
    else
        P=cell(0,1);
    end
elseif iP.Results.method=="vis3"
    [~,~,~,itertime_k3]=visClust_(DATAA,div,P3,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting simultaneous 3D Clustering, worst case runtime: "+num2str(itertime_k3*iternum_k3)+"sec")
    [idx,~,Pi]=visClust_(DATAA,div,P3,iP.Results.scaling,iP.Results.thresh,"visonly");
    if Pi>0
        P={P3{Pi}};
    else
        P=cell(0,1);
    end
elseif iP.Results.method=="pvis2"
    [~,~,~,itertime_k2]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting partial 2D Clustering, worst case runtime: "+num2str(itertime_k2*iternum_k2*(cl-1))+"sec")
    [idx,P]=visClustPartial_(DATAA,div,iP.Results.scaling,iP.Results.thresh,'pvisp2',iP.Results.projectors);
elseif iP.Results.method=="pvis3"
    [~,~,~,itertime_k3]=visClust_(DATAA,div,P2,iP.Results.scaling,iP.Results.thresh,"time");
    disp("Starting partial 3D Clustering, worst case runtime: "+num2str(itertime_k3*iternum_k3*(cl-1))+"sec")
    [idx,P]=visClustPartial_(DATAA,div,iP.Results.scaling,iP.Results.thresh,'pvis3',iP.Results.projectors);
end

%% Unify previously split dataset (in case of more than 10 000 data points)
if iP.Results.subsampling < l
    idxComplete=zeros(l,1);
    idxComplete(dicta)=idx;
    idxComplete(dicti)=idx(knnsearch(DATAA,DATAI));
    
    idx=idxComplete;
    targetvectorVals=unique(idx);
    for iout=1:length(targetvectorVals)
        div(iout)=sum(idx==targetvectorVals(iout))/length(iP.Results.X);
    end
end

%% Generate images
if size(P,1)>0 && iP.Results.showimg~="false"
    f=figure;
    numclust=unique(idx);
    for i=1:length(numclust)
        tdata=iP.Results.X(idx==numclust(i),:);
        if size(P{1,1},2)==2
            scatter(tdata*P2{Pi}(:,1),tdata*P2{Pi}(:,2),"*");
        elseif size(P{1,1},2)==3
            scatter3(tdata*P3{Pi}(:,1),tdata*P3{Pi}(:,2),tdata*P3{Pi}(:,3),"*");
        end
        hold on;
    end
    if iP.Results.showimg~="true"
        savefig(iP.Results.showimg);
        saveas(f,iP.Results.showimg+".png")
    end
end
end
