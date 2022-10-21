function P = getProj(DATA,lowdim,n)
%% INPUT
% DATA... stored in n1xn2 Matrix, n1 = #data points, n2 = original dimension
% lowdim... lower dimension of projection
% n... either string "PCA" or number of random projections
%% OUTPUT
% P... projections
%
%
% This is part of the clustering algorithm VISCLUST, written by Anna Breger 
% and Clemens Karner.
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
if nargin ==3
    if (isstring(n))&& (isnumeric(lowdim))
        if upper(n)=="PCA" % compute PCA for given data
            P = cell(1,1);
            [coeff,~] = pca(DATA);
            P{1} = coeff(:,1:lowdim);
        else
            error("value of n invalid")
        end
    elseif (isnumeric(n)) &&  (isnumeric(lowdim)) % create n random projections
        d=size(DATA,2);
        P = cell(n,1);
        for i=1:n
            [Qi,~] = qr(randn(d,lowdim),0);
            P{i} = Qi;
        end
    else
        error("arguments of wrong data type")
    end
else
    error("3 arguments required")
end
end