function [idx,PLIST] = visClustPartial_(DATA,x,s,thresh,method,Pinp,targetVectorInp)
% Partial visClust for 2D and 3D
%% input
% DATA ...stored in n1xn2 matrix, n1 = #data points, n2 = original dimension
% x... wanted division of the clusters
% s... scaling factor for sigma
% thresh... threshold for division
% method... clustering method, either "pvis" or "pvisp"
% targetVectorInp... matrix containing clusters from visClust
% Pinp... cell with projection(s)
%% output
% idx... class assignments stored in n1x1 Matrix, n1 = #data points
% PLIST... list of used projection(s)
%
%
% This is part of the clustering algorithm VISCLUST, written by Anna Breger
% and Clemens Karner.
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
%% Initialize variables
if isstring(Pinp) && (Pinp=="tsne2" ||Pinp=="tsne3")
    iternum_k2=30;
    iternum_k3=10;
else
    iternum_k2=3000;
    iternum_k3=1000;
end

PLIST=cell(0,1);
l=length(DATA);
cl=length(x);
x=sort(x);
outIndex=1;
idx=zeros(l,1);
targetMulti=zeros(l,0);
div=zeros(cl,1);
thresho=max([0.02,thresh/cl]); % adapt threshold to the number of clusters
numClustersFound=0;
projectorIndex=0;
smode=0;
scaling=false;

if iscell(Pinp)
    try
        P=Pinp(:,1); % cell array with projections to R^2
        if size(Pinp,2)==2
            P3=Pinp(:,2); % cell array with projections to R^3
            scaling=true;
        end
    catch
        error("invalid projector format")
    end
elseif Pinp=="random"
    if method(end)=="0"
        scaling=true;
        P3 = getProj(DATA,3,iternum_k3);
        smode=2;
    elseif method(end)=="2"
        smode=2;
    elseif method(end)=="3"
        smode=3;
    end
elseif Pinp=="tsne"
    if method(end)=="0"
        scaling=true;
        smode=12;
    elseif method(end)=="2"
        smode=12;
    elseif method(end)=="3"
        smode=13;
    end
end
method=method(1:end-1);

%% Store recycled clusters from visClust
if  exist("targetVectorInp","var")
    if method == "pvis"  % check if clusters are found once
        numClustersFound=length(unique(targetVectorInp(targetVectorInp~=0)));
        idx(targetVectorInp~=0)=targetVectorInp(targetVectorInp~=0);
    elseif method =="pvisp" % check if clusters are found twice in 2-dim
        for i=1:unique(targetVectorInp(targetVectorInp~=0))
            targetMulti(:,end+1)=(targetVectorInp==i);
        end
    end
end

%% Start partial clustering until all clusters are found
while numClustersFound<cl-1
    if smode==2
        P = getProj(DATA,2,iternum_k2); % cell array with random projections to R^2
        projectorIndex=0;
    elseif smode==3
        P = getProj(DATA,3,iternum_k3); % cell array with random projections to R^3
        projectorIndex=0;
    elseif smode==12
        P = "tsne2";
        P3 = "tsne3";
        projectorIndex=0;
    elseif smode==13
        P = "tsne3";
        projectorIndex=0;
    end
    foundOneEntry=false;
    numClustersFoundTemp=0;
    thresh=thresho/sum(x(numClustersFound+1:cl)); % adapt threshold to the number of clusters
    xscaled=x(numClustersFound+1:end)/sum(x(numClustersFound+1:cl)); % adapt distribution to the number of clusters
    % find cluster
    if method=="pvisp" % partial clustering with two projections per cluster
        dataIndices=(idx==0);
        indicesDict=1:l;
        indicesDict(~dataIndices)=[];
        while ~foundOneEntry
            [target,numClustersTemp,projectorIndex]=visClust_(DATA(dataIndices,:),xscaled,P,s,thresh,"pvis",projectorIndex+1);
            if projectorIndex ==0 && scaling % change projection dim to 3 and switch to pvis
                method="pvis";
                P=P3;
                smode=smode+1;
                scaling = false;
                projectorIndex=0;
                [~,~,~,itertime_k3]=visClust_(DATA(dataIndices,:),x,P,s,thresh,"time");
                disp("Starting partial 3D Clustering, worst case runtime: "+num2str(itertime_k3*length(P)*(cl-1))+"sec");
                break
            elseif projectorIndex==0
                return
            end
            targetTemp=zeros(l,numClustersTemp);
            for tempIndex=1:numClustersTemp % iterate through all detected clusters and store each one independently
                targetTemp(indicesDict(target==tempIndex),tempIndex)=1;
            end
            for iii=1:numClustersTemp
                foundEntry=false;
                % check if the iii-th found cluster corresponds to any of the previously found clusters up to accuracy 0.95
                for iiii=1:size(targetMulti,2)
                    if adjustedRandIndex(targetTemp(:,iii),targetMulti(:,iiii)) > 0.85
                        idx(targetTemp(:,iii)==1)=outIndex;
                        if ~isstring(P)
                            PLIST{end+1,1}=P{projectorIndex};
                        end
                        outIndex=outIndex+1;
                        foundEntry=true;
                        foundOneEntry=true;
                        numClustersFoundTemp=numClustersFoundTemp+1;
                        break
                    end
                end
                if ~foundEntry
                    targetMulti(:,end+1)=targetTemp(:,iii);
                end
            end
        end
    end
    if method=="pvis" % partial clustering with single projection per cluster
        dataIndices=(idx==0);
        indicesDict=1:l;
        indicesDict(~dataIndices)=[];
        [target, numClustersFoundTemp,projectorIndex] =visClust_(DATA(dataIndices,:),xscaled,P,s,thresh,"pvis",projectorIndex+1);
        if projectorIndex==0 && scaling % change projection dim to 3
            P=P3;
            scaling=false;
            smode=smode+1;
            projectorIndex=0;
            [~,~,~,itertime_k3]=visClust_(DATA,x,P,s,thresh,"time");
            disp("Starting partial 3D Clustering, worst case runtime: "+num2str(itertime_k3*length(P)*(cl-1))+"sec");
            [target,numClustersFoundTemp,projectorIndex] =visClust_(DATA(dataIndices,:),xscaled,P,s,thresh,"pvis",projectorIndex+1);
        elseif projectorIndex==0
            break
        end
        for tempIndex=1:numClustersFoundTemp % iterate through all found clusters and assign classes
            idx(indicesDict(target==tempIndex))=numClustersFound+tempIndex;
        end
        if ~isstring(P)
            PLIST{end+1,1}=P{projectorIndex};
        end
    end
    % increase offset by number of clusters found
    numClustersFound=numClustersFound+numClustersFoundTemp;
end

%% Assign remaining cluster
idx(idx==0)=cl;
for i=1:cl
    div(i)=sum(idx==i)/l;
end
end