function [GAUSSFILTER] = mygaussfilt_multi(sigma,n,k)
%% input:
% sigma... standard deviation
% n... filtersize
% k... dimension
%% output:
% gaussFilter... n-dim Gaussian filter of size n^k discretized in [-2sigma,2sigma]
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
x = linspace(-3*sigma,3*sigma,n);
% Create 1D gaussian distribution
GAUSSFILTER = 1/(sigma*sqrt(2*pi))*exp(-0.5*x.^2/sigma^2);
GAUSSFILTER = GAUSSFILTER';
fil_temp = GAUSSFILTER;
% Extend it to N-D
for i = 2:k
    fil_temp2 = GAUSSFILTER;
    GAUSSFILTER = GAUSSFILTER.*fil_temp(1);
    for j = 2:n
        GAUSSFILTER = cat(i,GAUSSFILTER,fil_temp2.*fil_temp(j));
    end
end
end

