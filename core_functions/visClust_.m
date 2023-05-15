function [idx,numClustersFound,PIndOut,time,PIndSingleOut,targetvectorSingle,success]= visClust_(DATA,x,P,s,thresh,method,PInd)
% visual clustering with projections from R^d to R^k with k=2,3
%% input
% DATA ... normalized data stored in n1xn2 Matrix, n1 = #data points, n2 = original dimension
% x... wanted division of the dataset,
% P... cell with orthogonal projection(s)
% s... scaling factor for sigma
% thresh... threshold for wanted division, default = 0.1
% method... clustering method "vis"... visClust, "pvis"... partial visClust clustering, "time"... measures average time per iteration based on three iterations
% PInd... projector index to begin iteration with
%% output
% idx... class assignments stored in n1x1 Matrix, n1 = #data points
% numClustersFound... number of clusters found
% PIndOut... index of chosen projector
% time... runtime of iteration
% PIndSingleOut... index of projection for single cluster
% targetvectorSingle... clustering for single cluster
% success... logical value indicating successfull clustering
%
%
% This is the main part of the clustering algorithm VISCLUST, written by
% Anna Breger and Clemens Karner.
% University of Vienna, Faculty of Mathematics
% Vienna, Austria
% Copyright (c) 2023
% https://homepage.univie.ac.at/anna.breger/
% https://homepage.univie.ac.at/clemens.karner/
%
% For all questions, bugs and suggestions please email
% clemens.karner@univie.ac.at or anna.breger@univie.ac.at
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input arguments
if nargin < 7
    PInd=1;
end

%% Normalize data and filter outliers
if isstring(P)
    if P=="tsne2" || P=="tsne3"
        outIndices = [];
    end
else
    outIndices = find(any([DATA < mean(DATA)-4*std(DATA),DATA > mean(DATA)+4*std(DATA)],2))';
end
regularPtsIndices=1:size(DATA,1);
regularPtsIndices(outIndices)=[];
OUTLIERS = DATA(outIndices,:);
DATA(outIndices,:) = [];
scale=max(abs(max(DATA,[],1)),abs(min(DATA,[],1)));
scale(scale==0)=1;
NDATA = DATA./scale; % normalize data
OUTLIERS= OUTLIERS./scale; % normalize data

%% Initialize variables
success=false;
x=sort(x);
time=0;
assignTemporaryClusters=false;
PIndSingleOut=0;
PIndOut = 0;
targetvectorSingle=-1;
l=length(NDATA); % number of instances without outliers
if isstring(P)
    if P=="tsne2"
        m=10;
        k = 2; % projected dimension
    elseif P=="tsne3"
        m=10;
        k = 3; % projected dimension
    end
else
    m = length(P); % number of projections
    k = size(P{1},2); % projected dimension
end
cl = length(x); % number of clusters
res = 100; % resolution identifies the number of digits taken into account
numClustersFound=0;
wrongClust=0;
wrongClustIter=0;
threshOpt=1;

if method=="time" % measuring expected time of visClust
    timeArray=zeros(10,1);
    m=min(10,size(P,1));
elseif method =="vis" % actual clustering
    assignTemporaryClusters=true; % recycle clusters from visClust
end

%% Iterate over all projections until clusters are found
for i = PInd:m
    if method=="time"
        tic;
    end
    %% Iterative adaption of scaling factor
    wrongClustIter=wrongClustIter+1;
    if wrongClustIter == 250
        if wrongClust >200
            s=s*1.25;
        elseif wrongClust < -200
            s=s*0.75;
        end
        wrongClustIter=0;
        wrongClust=0;
    end
    %% Project and scale data
    if isstring(P)
        if P=="tsne2" || P=="tsne3"
            CDATA=tsne(NDATA,'NumDimensions',k);
            CDATA=CDATA/max(CDATA(:));
        end
    else
        CDATA = NDATA*P{i};
    end
    BDATA = round(CDATA*res);
    %% Compute appropriate sigma for Gaussian filter from distances in the projected data set
    if l > 500
        dist = sort(pdist(CDATA(randi(l,[500,1]),:)))./(2*(k-1)-1);
        dist(dist==0) = [];
        sigma = s*median(dist(1:1000));
    else
        dist = sort(pdist(CDATA))./(2*(k-1)-1);
        dist(dist==0) = [];
        sigma = s*median(dist(1:l))*l/250;
    end
    n = 2*ceil(res*3*sigma)+1; % odd filter size
    TDATA = BDATA - min(BDATA,[],1) + n; %scale projected data for filtering
    %% Compute and apply Gaussian filter
    GAUSSFILTER = mygaussfilt_multi(sigma,n,k);
    M1 = zeros(max(TDATA)+n);
    switch k
        case 1
            for j = 1:l
                M1(TDATA(j,1)-(n-1)/2:TDATA(j,1)+(n-1)/2) = M1(TDATA(j,1)-(n-1)/2:TDATA(j,1)+(n-1)/2) + (1/l)*GAUSSFILTER;
            end
        case 2
            for j = 1:l
                M1(TDATA(j,1)-(n-1)/2:TDATA(j,1)+(n-1)/2,TDATA(j,2)-(n-1)/2:TDATA(j,2)+(n-1)/2) = M1(TDATA(j,1)-(n-1)/2:TDATA(j,1)+(n-1)/2,TDATA(j,2)-(n-1)/2:TDATA(j,2)+(n-1)/2) + (1/l)*GAUSSFILTER;
            end
        case 3
            for j = 1:l
                M1(TDATA(j,1)-(n-1)/2:TDATA(j,1)+(n-1)/2,TDATA(j,2)-(n-1)/2:TDATA(j,2)+(n-1)/2,TDATA(j,3)-(n-1)/2:TDATA(j,3)+(n-1)/2) = M1(TDATA(j,1)-(n-1)/2:TDATA(j,1)+(n-1)/2,TDATA(j,2)-(n-1)/2:TDATA(j,2)+(n-1)/2,TDATA(j,3)-(n-1)/2:TDATA(j,3)+(n-1)/2) + (1/l)*GAUSSFILTER;
            end
        otherwise
            error("Projection dimension invalid, only integers between 1 and 3 are supported.")
    end
    [M2,n3] = bwlabeln(M1 > mean(M1(:))); % identify connected components

    %% Discard regions evolving from outliers
    nc = histcounts(M2,n3+1) >= n^k;
    v = (0:n3);
    cy = nonzeros(v(nc)); % clusters without outliers
    cn = [0,v(~nc)]; % outlier clusters & background

    %% Count number of wrong clusters
    if length(cy) > cl
        wrongClust = wrongClust+1;
    elseif length(cy)< cl
        wrongClust = wrongClust-1;
    end

    %% Measure time
    if method=="time"
        timeArray(i)=toc;
        %% Partially assign clusters
    elseif method=="pvis" || assignTemporaryClusters
        if length(cy) >= 2 || m == 1 % enter when at least two clusters are found or just one projection is provided
            CDATATEMP = num2cell(TDATA,1);
            outputTemp = M2(sub2ind(size(M2),CDATATEMP{:}));
            % combine data with outliers
            output=zeros(l+size(OUTLIERS,1),1);
            output(regularPtsIndices)=outputTemp;
            CDATA=zeros(size(output,1),k);
            CDATA(regularPtsIndices,:)=TDATA;
            if ~isstring(P)
                CDATA(outIndices,:)=OUTLIERS*P{i};
            end
            if nnz(output) > l/10
                distr = zeros(cl,1);
                % nearest neighbor for data points without assignment
                if sum(ismember(output,cn)) > 0 % check if there are data points without assignment
                    if l > 10000
                        M3=M2;
                        M3(ismember(M3,v(~nc)))=0;
                        if size(P{1},2)==2
                            edgedata = bwboundaries(M3,'noholes');
                            lenedge=sum(cellfun(@(c) size(c,1), edgedata));
                            edgepoints=zeros(lenedge,2);
                            edgelabels=zeros(lenedge,1);
                            lentemp=1;
                            for ie=1:length(edgedata)
                                clusttemp=M3(sub2ind(size(M3),edgedata{ie}(1,1),edgedata{ie}(1,2)));
                                edgepoints(lentemp:lentemp+size(edgedata{ie},1)-1,:)=edgedata{ie};
                                edgelabels(lentemp:lentemp+size(edgedata{ie},1)-1,:)=clusttemp*ones(size(edgedata{ie},1),1);
                                lentemp=lentemp+size(edgedata{ie},1);
                            end
                        elseif size(P{1},2)==3
                            E1=edge3(M3,"Sobel",0.1);
                            [E2,nE2]=bwlabeln(E1);
                            edgepoints=zeros(sum(E1(E1==1)),3);
                            edgelabels=zeros(sum(E1(E1==1)),1);
                            lentemp=1;
                            for ni=1:nE2
                                [i1,i2,i3] = ind2sub(size(E2),find(E2 == ni));
                                edgepoints(lentemp:lentemp+size(i1,1)-1,:)=[i1,i2,i3];
                                E3=M3(E2==ni);
                                edgelabels(lentemp:lentemp+size(i1,1)-1,:)=mode(E3(E3>0))*ones(size(i1,1),1);
                                lentemp=lentemp+size(i1,1);
                            end
                        end
                        while sum(ismember(output,cn)) > 0
                            index =dsearchn(edgepoints,CDATA(ismember(output,cn),:));
                            output((ismember(output,cn))) = edgelabels(index);
                        end
                    else
                        hdata = CDATA;
                        houtput = output;
                        while sum(ismember(output,cn)) > 0
                            hdata(ismember(output,cn),:) = [];
                            houtput(ismember(output,cn)) = [];
                            index =dsearchn(hdata,CDATA(ismember(output,cn),:));
                            output(ismember(output,cn)) = houtput(index);
                        end
                    end
                end
                outputTemp=output;
                outputVals=unique(outputTemp);
                for iout=1:length(outputVals)
                    distr(iout)=sum(outputTemp==outputVals(iout))/l;
                    output(outputTemp==outputVals(iout))=iout;
                end
                distrt=distr;
                if cl==2 % if only two clusters are searched, find the larger cluster first
                    numClustersFound=1;
                    [optmax,optid]=min(abs(distrt-x(2)));
                    if m ~= 1 && optmax < thresh
                        PIndOut=i;
                        if assignTemporaryClusters % assign temporary clusters
                            targetvectorSingle=zeros(length(output),1);
                            targetvectorSingle(output==optid)=1;
                            assignTemporaryClusters=false;
                            PIndSingleOut=i;
                        else % assign clusters
                            idx=ones(length(output),1);
                            idx(output==optid)=2;
                            return;
                        end
                    end
                else % if more than two clusters are searched, find the small cluster first
                    finished=false;
                    targetvectortemp=zeros(length(output),1);
                    cf=0;
                    for xtempindex=1:cl-1 % if more clusters than the smallest one are found, store them
                        [opt,optid]=min(abs(distrt-x(xtempindex)));
                        if m ~= 1 && opt < thresh
                            targetvectortemp(output==optid)=xtempindex;
                            distrt(optid)=100;
                            finished=true;
                            cf=cf+1;
                        else
                            break;
                        end
                    end
                    if finished
                        if assignTemporaryClusters % assign temporary clusters
                            assignTemporaryClusters=false;
                            PIndSingleOut=i;
                            targetvectorSingle=targetvectortemp;
                        else % assign clusters
                            idx=targetvectortemp;
                            idx(idx==0)=cl;
                            numClustersFound=cf;
                            PIndOut = i;
                            return;
                        end
                    end
                end
            end
        end
    end
    %% Simultaneously assign clusters
    if method == "vis" || method =="visonly"
        if length(cy) == cl || m == 1 % enter when correct number of clusters is found or just one projection is provided
            CDATATEMP = num2cell(TDATA,1);
            outputTemp = M2(sub2ind(size(M2),CDATATEMP{:}));
            % combine data with outliers
            output=zeros(l+size(OUTLIERS,1),1);
            output(regularPtsIndices)=outputTemp;
            CDATA=zeros(size(output,1),k);
            CDATA(regularPtsIndices,:)=TDATA;
            if ~isstring(P)
                CDATA(outIndices,:)=OUTLIERS*P{i};
            end

            if nnz(output) > l/10 || m==1
                distr = zeros(cl,1);
                % nearest neighbor for data points without assignment
                if sum(ismember(output,cn)) > 0 % check if there are data points without assignment
                    if l > 10000
                        M3=M2;
                        M3(ismember(M3,v(~nc)))=0;
                        if size(P{1},2)==2
                            edgedata = bwboundaries(M3,'noholes');
                            lenedge=sum(cellfun(@(c) size(c,1), edgedata));
                            edgepoints=zeros(lenedge,2);
                            edgelabels=zeros(lenedge,1);
                            lentemp=1;
                            for ie=1:length(edgedata)
                                clusttemp=M3(sub2ind(size(M3),edgedata{ie}(1,1),edgedata{ie}(1,2)));
                                edgepoints(lentemp:lentemp+size(edgedata{ie},1)-1,:)=edgedata{ie};
                                edgelabels(lentemp:lentemp+size(edgedata{ie},1)-1,:)=clusttemp*ones(size(edgedata{ie},1),1);
                                lentemp=lentemp+size(edgedata{ie},1);
                            end
                        elseif size(P{1},2)==3
                            E1=edge3(M3,"Sobel",0.1);
                            [E2,nE2]=bwlabeln(E1);
                            edgepoints=zeros(sum(E1(E1==1)),3);
                            edgelabels=zeros(sum(E1(E1==1)),1);
                            lentemp=1;
                            for ni=1:nE2
                                [i1,i2,i3] = ind2sub(size(E2),find(E2 == ni));
                                edgepoints(lentemp:lentemp+size(i1,1)-1,:)=[i1,i2,i3];
                                E3=M3(E2==ni);
                                edgelabels(lentemp:lentemp+size(i1,1)-1,:)=mode(E3(E3>0))*ones(size(i1,1),1);
                                lentemp=lentemp+size(i1,1);
                            end
                        end
                        while sum(ismember(output,cn)) > 0
                            index =dsearchn(edgepoints,CDATA(ismember(output,cn),:));
                            output((ismember(output,cn))) = edgelabels(index);
                        end
                    else
                        hdata = CDATA;
                        houtput = output;
                        while sum(ismember(output,cn)) > 0
                            hdata(ismember(output,cn),:) = [];
                            houtput(ismember(output,cn)) = [];
                            index =dsearchn(hdata,CDATA(ismember(output,cn),:));
                            output(ismember(output,cn)) = houtput(index);
                        end
                    end
                end
                outputTemp=output;
                outputVals=unique(outputTemp);
                for iout=1:length(outputVals)
                    distr(iout)=sum(outputTemp==outputVals(iout))/l;
                    output(outputTemp==outputVals(iout))=iout;
                end
                if m==1 % if only one projector is provided assign classes and terminate
                    numClustersFound=1;
                    PIndOut=i;
                    idx=output;
                    return;
                else
                    distrsort=sort(distr);
                    tempthresh=sum(abs(x-distrsort'));
                    if tempthresh<threshOpt % if a new optimal clustering is found assign classes
                        numClustersFound=1;
                        PIndOut=i;
                        idx=output;
                        threshOpt=tempthresh;
                        if tempthresh<thresh % if an optimal clustering within the searched threshold is found terminate
                            success=true;
                            return;
                        end
                    end
                end
            end
        end
    end
end
if method=="time"
    idx=zeros(l+size(OUTLIERS,1),1);
    if m>=3
        time=mean(timeArray(3:m));
    else
        time=timeArray(1:m);
    end
end
if PIndOut == 0
    idx=zeros(l+size(OUTLIERS,1),1);
end
end
